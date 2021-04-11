;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writepds,outfile,chan=chan,help=help,_extra=_ext

; Write the loaded pair to disk as an pds file

common loadedBlock,loadedi,loaded1,loaded
common pdsBlock,pds

if n_params() ne 1 or keyword_set(help) then begin
  print,'writepds,outfile[,/overwrite][,/append][,chan=1 or 2][,/help]'
  print,"         '.csv' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writepds,outfile[,/overwrite][,/append][,chan=1 or 2][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loaded).ndata le 2 then begin
  print,'ERROR in writepds: No pair is loaded, so nothing can be written to disk'
  return
endif

if size(pds) le 0 then begin
  print, "pds structure not set"
  return
end

pairpointer = loaded

writeapdsfile,pairpointer,outfile,chan=chan

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setpds,imfile=infile,target=target,pname=pname,ttype=ttype,author=author,level=level,bookmark=bookmark, arecibo=arecibo,dss14=dss14

common pdsBlock,pds

;default
if size(pds) eq 0 then pds = {infile: '', target:'', pname:'',ttype:'Asteroid', author:'Planetary Radar Team', level:'Calibrated', bookmark:'AO TX;AO RX'}
if keyword_set(arecibo) then begin
  pds.bookmark = 'AO TX;AO RX;AO RI'
end
if keyword_set(dss14) then begin
  pds.bookmark = 'DSS14 TX;DSS14 RX'
endif
if keyword_set(infile) then pds.infile = infile
if keyword_set(target) then pds.target = target
if keyword_set(pname) then pds.pname = pname
if keyword_set(ttype) then pds.ttype = ttype
if keyword_set(author) then pds.author = author
if keyword_set(level) then pds.level = level

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getextral,extratags,name,help=help

; Returns the value of an extra tag directly from the structure
;
; /comment returns the comment rather than the tag value
;
; Returns a null string if the tag isn't found.
;
; 2006 Jul 7: Add 'comment' keyword

if size(name, /type) ne 7 then name = ''  ;  just so it's defined as a string

if keyword_set(help) or isnull(name) or n_params() ne 1 then begin
  print,' '
  print,'value = getextral(extratags,name[,stack=n][,/comment][,/help])'
  print,' '
  print,'Return the value of an extra tag from the pointer;'
  print,'    if the tag is not present, the null string is returned'
  print,''
  print,'name must be a quoted string'
  print,' '
  print,'/comment returns the tag comment rather than the tag value'
  print,' '
  return, ''
endif

; Search and retrieve the value

extranum = where(strlowcase(extratags.name) eq strlowcase(name), count)

if count gt 0 then begin

  ; Tag present

  extranum = extranum[0]
  if keyword_set(comment) then begin
    value = extratags[extranum].comment
  endif else begin
    case strlowcase(extratags[extranum].format) of
      's' : value = extratags[extranum].value
      'c' : value = extratags[extranum].value
      't' : value = long(extratags[extranum].value)
      'i' : value = long(extratags[extranum].value)
      'f' : value = float(extratags[extranum].value)
      'd' : value = double(extratags[extranum].value)
      'v' : value = extratags[extranum].value  ;  just treat vectors as strings
      else: begin
              print,"ERROR: getextral isn't set up to handle extra tag format '", $
                    extratags[extranum].format,"'",format='(3a)'
              return,''
            endcase
    endcase
  endelse
endif else begin

  ; Tag not present

  value = ''

endelse

return, value
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writeapdsfile,pairpointer,outfile,pds,tablelun,chan=chan
; This will actually do the writing for one file.
; First implement the writing of the data, TODO add a table

;
; PDS csv files. One pair per file. Starts with PDS keywords, then the column descriptions, then the data as three csv columns: frequency,pol1,pol2
; Should we allow one channel? I don't think it's needed, but should 
; almost work automagically.
; append will only mean anything if we make a table. Good idea?
; Need a way to specify bistatic rx station, nonstandard receiver or datataking setup.

common pdsBlock,pds

if n_params() lt 3 or keyword_set(help) then begin
  print, 'writepdsfile,pair,outfile,pds[,tablelun][,/overwrite][,chan=1 or 2][,/help]'
  return
endif

; pairpointter is a pointer to either the loaded pair or a stack entry. I *think* they look the same.
;pds is a structure containg the PDS keywords that aren't natively in the data

err = openOutfile(lun,outfile,'csv',/get_lun,_extra=_ext)
if err ne 0 then return

; Check which channel(s) to write

if n_elements(chan) eq 0 then begin
  writepol = [1,1]
  threecol = 1
  nchan = 1
endif else if chan eq 1 or chan eq 2 then begin
  writepol = (chan eq 1) ? [1,0] : [0,1]
  threecol = 0
  nchan = 2
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse

; Get some elements of the loaded spectrum
pair = (*pairpointer).spec
extratags = (*pairpointer).extratags
tags = (*pairpointer).tags
ndata = (*pairpointer).ndata
ntags = (*pairpointer).ntags
tname = (*pairpointer).tname
nextra = (*pairpointer).nextra
if nchan eq 2 then begin
  extratags_write = extratags
  nextra_write = nextra
  addacomma = ','
endif else begin
  splitExtra,chan,extratags,nextra,extratags_write,nextra_write
  addacomma=''
endelse

; Get the times Assuming jdstart and end are correct for now.
; convert to ISO. Do start last so we have start times in vars
jdend = getextral(extratags,'jdend')
caldat, jdend, mon,day,year,hh,mm,ss
ss = int(ss + 0.5) ; round to nearest to aboid jd rounding
endstring = timestamp(day=day,hour=hh,min=mm,month=mon,second=ss,year=year,/utc)
jdstart = getextral(extratags,'jdstart')
caldat, jdstart, mon,day,year,hh,mm,ss
ss = int(ss + 0.5) ; round to nearest to aboid jd rounding
startstring = timestamp(day=day,hour=hh,min=mm,month=mon,second=ss,year=year,/utc)

; Construct and print some PDS header strings
printf, lun, '#Keywords,',addacomma
printf, lun, 'Product Name,',pds.pname,addcomma
printf, lun, 'Product Description,CW spectrum converted from RDF format',addcomma
printf, lun, 'Start Time,', startstring,addcomma
printf, lun, 'Stop Time,', endstring, addcomma
printf, lun, 'Target Name,', pds.target, addcomma
printf, lun, 'Target Type,', pds.targettype, addcomma
printf, lun, 'Author List,', pds.author, addcomma
printf, lun, 'Product Processing Level,', pds.proclevel, addcomma
printf, lun, 'Science Search Facet,', pds.facet, addcomma
printf, lun, 'Product Wavelength Ranges,Microwave', addcomma
; next should be replaced with istruments and telescopes
printf, lun, 'Observing System Bookmark,', pds.bookmark, addcomma 
printf, lun, 'CW data file,', pds.infile, addcomma
;
; Now tags
;
; Should allow override
skiptags=['iyy','imm','idd','rchour','rcmin','rcsec']
for i = 0, ntags-1 do begin
dummy=where(strlowcase(tname[i]) eq skiptags, count)
if count gt 0 then continue
if threecol then begin
  tn = string(tname[i], format='(a,"_1,")')
  printf, lun, tn, tags[0].(i), addcomma
  tn = string(tname[i], format='(a,"_2,")')
  printf, lun, tn, tags[1].(i), addcomma
endif else begin
  tn = string(tname[i], format='(a,",")')
  printf, lun, tn, tags[chan].(i)
endelse

;
; and extra tags
;

for i = 0, nextra-1 do begin
  if (extratags[i].format eq 't') then continue ; skip tag names
  printf, lun, extratags[i].name, ',', extratags[i].value, addcomma
endfor

;Column definitions

printf, lun, "# Column Definitions,", addcomma
if threcol then begin
  printf, lun, 'frequency,pol1,pol2'
  printf, lun, 'real,real,real'
  printf, lun, 'Hz,,'
  printf, lun, 'Offset from Ephemeris frequency,Polarization 1 normalized to unit standard deviation and zero mean in baseline,Polarization 2 normalized to unit standard deviation and zero mean in baseline'
endif else begin
  printf, lun, 'frequency,pol1'
  printf, lun, 'real,real'
  printf, lun, 'Hz,'
  printf, lun, 'Offset from Ephemeris frequency,Polarization 1 normalized to unit standard deviation and zero mean in baseline'
endelse

;
; and now the data
;
printf, lun, '# Data,', addcomma

dfreq = tags[0].dfreq
posfr = tags[0].posfr
xjcen = tags[0].xjcen
freq = posfr*dfreq*(dindgen(ndata) - xjcen)
pair = (*stack[n-1]).spec


if threecol then begin
  for i = 0, ndata-1 do begin
    printf, lun, dfreq,pair[0,i],pair[1,i], format='(e0,e0,e0)'
  endfor
endif else begin
  for i = 0, ndata-1 do begin
    printf, lun, dfreq,pair[0,i],format='(e0,e0)'
  endfor
endelse

close, lun
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writepds1,outfile,tkplay=tkplay,help=help,_extra=_ext

; Write the loaded single-channel spectrum to disk as an pds file

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 1 or keyword_set(help) then begin
  print,'writepds1,outfile[,/overwrite][,/append][,/tkplay][,/help]'
  print,"          '.csv' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writepds1,outfile[,/overwrite][,/append][,/tkplay][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loaded1).ndata le 2 then begin
  print,'ERROR in writepds1: No single-channel spectrum is loaded,', $
        ' so nothing can be written to disk',format='(2a)'
  return
endif

; Get some elements of the loaded spectrum

spec = (*loaded1).spec
tags1 = (*loaded1).tags
extratags = (*loaded1).extratags
ndata = (*loaded1).ndata
ntags = (*loaded1).ntags
nextra = (*loaded1).nextra
tname = (*loaded1).tname
pol = (*loaded1).pol

; If requested, fix three tag values (xjcen, jsnr1, jsnr2) so that the output
; is consistent with tkplay, which counts these tags from 1 rather than from 0
;
; Also fix the jcp tag, which tkplay counts (sometimes) from 0 rather than from 1

if keyword_set(tkplay) then begin
  tags1.xjcen = tags1.xjcen + 1L
  tags1.jsnr1 = tags1.jsnr1 + 1L
  tags1.jsnr2 = tags1.jsnr2 + 1L
  tags1.jcp = tags1.jcp - 1L
  if tags1.jcp lt 0 then begin
    print,'ERROR in writepds1: jcp tag for OC is already 0 rather than 1'
    return
  endif
  extratags = [extratags, extratags[0]]
  extratags[nextra].format = ''
  extratags[nextra].name = ''
  extratags[nextra].value = ''
  extratags[nextra].comment = $
               '# Conv to   tkplay format ' + systime(/utc) + ' UT'
  nextra = nextra + 1L
endif

; Open the output file

err = openOutfile(lun,outfile,'csv',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the pds header

sformat = '(a,5x,a16,3x,a)'
iformat = '(a,5x,a16,3x,i0)'
printf,lun,'s','type','float',format=sformat
printf,lun,'i','ndata',ndata,format=iformat
printf,lun,'i','height',1,format=iformat
printf,lun,'i','width',ndata+ntags,format=iformat
printf,lun,'i','size',4,format=iformat
printf,lun,'s','machine','SPARC',format=sformat
printf,lun,'i','nchan',1,format=iformat
printf,lun,'i','ntags',ntags,format=iformat
printf,lun,'i','nspec',1,format=iformat
printf,lun,'s','format','RDF_spectrum_1',format=sformat
printf,lun,'.'

; Write the spectrum and tags (as binary)
; -- write all tags as floating point to be compatible with tkplay

writeu,lun,swap_endian(spec, /swap_if_little_endian)
for n=0L,ntags-1 do begin
  writeu,lun,swap_endian(float(tags1.(n)), /swap_if_little_endian)
endfor

; Write the pds footer

for n=0L,nextra-1 do begin
  printf,lun,strtrim(string(extratags[n],format='(a,5x,a16,3x,a,3x,a)'), 2)
endfor
printf,lun,'s','target',tname,format=sformat
printf,lun,'s','polarization',pol,format=sformat
printf,lun,'.'

print,'Wrote 1 pds frame, containing 1 ',pol, $
      ' spectrum, to file ',outfile,format='(4a)'

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writestackpds,outfile,group=g,ming=ming,maxg=maxg,arrayg=ag, $
                  chan=chan,tkplay=tkplay,help=help,_extra=_ext

; Take stack pairs in one or more groups and write them to disk as an pds file,
; one pds frame per pair.  The keywords g, ming and maxg, or ag specify the groups;
; omitting all of these causes all stack pairs to be written.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check if a valid set of group-related keywords is set

if n_elements(g) gt 0 then begin
  gKeywordsOK = n_elements(ming) eq 0 and n_elements(maxg) eq 0 and n_elements(ag) eq 0
endif else if n_elements(ming) gt 0 then begin
  gKeywordsOK = n_elements(maxg) gt 0 and n_elements(ag) eq 0
endif else if n_elements(maxg) gt 0 then begin
  gKeywordsOK = 0
endif else begin
  gKeywordsOK = 1
endelse

if n_params() ne 1 or (not gKeywordsOK) or keyword_set(help) then begin
  print,' '
  print,'writestackpds,outfile[,/overwrite][,/append]'
  print,'              [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'              [,chan=1 or 2][,/tkplay][,/help]'
  print,' '
  print,"    '.csv' output file extension is added if not already present"
  print,' '
  print,'    Setting keywords g, ming and maxg, or ag writes all pairs in ', $
        'the specified group(s) to an pds file',format='(2a)'
  print,'    -- these three keyword choices are mutually exclusive'
  print,'    Calling writestackpds with none of those keywords set writes ', $
        'all stack pairs to an pds file',format='(2a)'
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writestackpds,outfile[,/overwrite][,/append]'
  print,'              [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'              [,chan=1 or 2][,/tkplay][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if nstack eq 0 then begin
  print,'ERROR in writestackpds: The pair stack is empty'
  return
endif

; Create a sorted array of the group numbers which will be written

allgroups = stackgroups()
n_all = n_elements(allgroups)

if n_elements(g) gt 0 then begin
  ngroups = 1L
  groupvec = [long(g)]
endif else if n_elements(ming) gt 0 then begin
  if ming gt maxg then begin
    print,'ERROR in writestackpds: Must specify range with ming <= maxg'
    return
  endif
  ngroups = long(maxg) - long(ming) + 1L
  groupvec = lonarr(n_all)
  ngroups = 0L
  for k=0L,n_all-1 do begin
    if allgroups[k] ge ming and allgroups[k] le maxg then begin
      groupvec[ngroups] = allgroups[k]
      ngroups = ngroups + 1L
    endif
  endfor
  if ngroups eq 0 then begin
    print,'writestackpds: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
endif else if n_elements(ag) gt 0 then begin
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroups then begin
    print,"ERROR in writestackpds: Don't include duplicate group numbers ", $
          "in array ag",format='(2a)'
    return
  endif
endif else begin
  ngroups = n_all
  groupvec = allgroups
endelse

; Check that all specified groups actually exist

for j=0L,ngroups-1 do begin
  groupnum = where(allgroups eq groupvec[j], count)
  if count eq 0 then begin
    print,'ERROR in writestackpds: Group #',groupvec[j], $
          ' is not present in the pair stack',format='(a,i0,a)'
    return
  endif
endfor

; Check which channel(s) to write

if n_elements(chan) eq 0 then begin
  writepol = [1,1]
endif else if chan eq 1 or chan eq 2 then begin
  writepol = (chan eq 1) ? [1,0] : [0,1]
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse
nchan = total(writepol)

; Open the output file

err = openOutfile(lun,outfile,'csv',/get_lun,_extra=_ext)
if err ne 0 then return

; Go through the stack pair by pair and write out pairs which
; are included in one of the specified groups

sformat = '(a,5x,a16,3x,a)'
iformat = '(a,5x,a16,3x,i0)'
nwrite = 0L

for n=0L,nstack-1 do begin

  ; See if this pair is included in one of the groups which will be written

  group_el = where(groupvec eq (*stack[n]).group, count)

  if count gt 0 then begin

    ; This pair should be written to disk

    pair = (*stack[n]).spec
    tags = (*stack[n]).tags
    extratags = (*stack[n]).extratags
    ndata = (*stack[n]).ndata
    ntags = (*stack[n]).ntags
    nextra = (*stack[n]).nextra
    tname = (*stack[n]).tname
    if nchan eq 2 then begin
      extratags_write = extratags
      nextra_write = nextra
    endif else begin
      splitExtra,chan,extratags,nextra,extratags_write,nextra_write
    endelse

    ; If requested, fix three tag values (xjcen, jsnr1, jsnr2) so that the output
    ; is consistent with tkplay, which counts these tags from 1 rather than from 0
    ;
    ; Also fix the jcp tag, which tkplay counts (sometimes) from 0 rather than from 1

    if keyword_set(tkplay) then begin
      tags.xjcen = tags.xjcen + 1L
      tags.jsnr1 = tags.jsnr1 + 1L
      tags.jsnr2 = tags.jsnr2 + 1L
      tags.jcp = tags.jcp - 1L
      if (tags[0].jcp lt 0 or tags[1].jcp lt 0) then begin
        print,'ERROR in writestackpds: jcp tags are already [0,1] rather than [1,2]'
        return
      endif
      extratags_write = [extratags_write, extratags_write[0]]
      extratags_write[nextra_write].format = ''
      extratags_write[nextra_write].name = ''
      extratags_write[nextra_write].value = ''
      extratags[nextra_write].comment = $
                   '# Conv to   tkplay format ' + systime(/utc) + ' UT'
      nextra_write = nextra_write + 1L
    endif

    ; Write the pds header

    printf,lun,'s','type','float',format=sformat
    printf,lun,'i','ndata',ndata,format=iformat
    printf,lun,'i','height',nchan,format=iformat
    printf,lun,'i','width',ndata+ntags,format=iformat
    printf,lun,'i','size',4,format=iformat
    printf,lun,'s','machine','SPARC',format=sformat
    printf,lun,'i','nchan',nchan,format=iformat
    printf,lun,'i','ntags',ntags,format=iformat
    printf,lun,'i','nspec',nchan,format=iformat
    printf,lun,'s','format','RDF_spectrum_1',format=sformat
    printf,lun,'.'

    ; Write the spectrum and tags (as binary)
    ; -- write all tags as floating point to be compatible with tkplay

    for ch=1,2 do begin
      if writepol[ch-1] then begin
        writeu,lun,swap_endian(pair[ch-1,*], /swap_if_little_endian)
        for k=0L,ntags-1 do begin
          writeu,lun,swap_endian(float(tags[ch-1].(k)), /swap_if_little_endian)
        endfor
      endif
    endfor

    ; Write the pds footer

    for k=0L,nextra_write-1 do begin
      printf,lun,strtrim(string(extratags_write[k],format='(a,5x,a16,3x,a,3x,a)'), 2)
    endfor
    printf,lun,'s','target',tname,format=sformat
    chanstring = ['OC','SC']
    if nchan eq 1 then printf,lun,'s','polarization',chanstring[chan-1],format=sformat
    printf,lun,'.'

    nwrite = nwrite + 1L

  endif
endfor

if nchan eq 2 then begin
  print,'Wrote ',nwrite,' pds frames, each containing 1 OC and 1 SC spectrum, to file ', $
        outfile,format='(a,i0,2a)'
endif else begin
  print,'Wrote ',nwrite,' pds frames, each containing 1 ',chanstring[chan-1], $
        ' spectrum, to file ',outfile,format='(a,i0,4a)'
endelse

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writestackpds1,outfile,group=g,ming=ming,maxg=maxg,arrayg=ag, $
                   chan=chan,tkplay=tkplay,help=help,_extra=_ext

; Take single-channel stack spectra in one or more groups and write them to disk
; as an pds file, one pds frame per spectrum.  The keywords g, ming and maxg,
; or ag specify the groups; omitting all of these causes all stack1 spectra to 
; be written.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check if a valid set of group-related keywords is set

if n_elements(g) gt 0 then begin
  gKeywordsOK = n_elements(ming) eq 0 and n_elements(maxg) eq 0 and n_elements(ag) eq 0
endif else if n_elements(ming) gt 0 then begin
  gKeywordsOK = n_elements(maxg) gt 0 and n_elements(ag) eq 0
endif else if n_elements(maxg) gt 0 then begin
  gKeywordsOK = 0
endif else begin
  gKeywordsOK = 1
endelse

if n_params() ne 1 or (not gKeywordsOK) or keyword_set(help) then begin
  print,' '
  print,'writestackpds1,outfile[,/overwrite][,/append]'
  print,'               [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'               [,chan=1 or 2][,/tkplay][,/help]'
  print,' '
  print,"    '.csv' output file extension is added if not already present"
  print,' '
  print,'    Setting keywords g, ming and maxg, or ag writes all single-channel ', $
        'spectra in the specified group(s) to an pds file',format='(2a)'
  print,'    -- these three keyword choices are mutually exclusive'
  print,'    Calling writestackpds1 with none of those keywords set writes ', $
        'all stack1 spectra to an pds file',format='(2a)'
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writestackpds1,outfile[,/overwrite][,/append]'
  print,'               [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'               [,chan=1 or 2][,/tkplay][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if nstack1 eq 0 then begin
  print,'ERROR in writestackpds1: The single-channel stack is empty'
  return
endif

; Create a sorted array of the group numbers which will be written

allgroups = stackgroups1()
n_all = n_elements(allgroups)

if n_elements(g) gt 0 then begin
  ngroups = 1L
  groupvec = [long(g)]
endif else if n_elements(ming) gt 0 then begin
  if ming gt maxg then begin
    print,'ERROR in writestackpds1: Must specify range with ming <= maxg'
    return
  endif
  ngroups = long(maxg) - long(ming) + 1L
  groupvec = lonarr(n_all)
  ngroups = 0L
  for k=0L,n_all-1 do begin
    if allgroups[k] ge ming and allgroups[k] le maxg then begin
      groupvec[ngroups] = allgroups[k]
      ngroups = ngroups + 1L
    endif
  endfor
  if ngroups eq 0 then begin
    print,'writestackpds1: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
endif else if n_elements(ag) gt 0 then begin
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroups then begin
    print,"ERROR in writestackpds1: Don't include duplicate group numbers ", $
          "in array ag",format='(2a)'
    return
  endif
endif else begin
  ngroups = n_all
  groupvec = allgroups
endelse

; Check that all specified groups actually exist

for j=0L,ngroups-1 do begin
  groupnum = where(allgroups eq groupvec[j], count)
  if count eq 0 then begin
    print,'ERROR in writestackpds1: Group #',groupvec[j], $
          ' is not present in the pair stack',format='(a,i0,a)'
    return
  endif
endfor

; Check which channel(s) to write

if n_elements(chan) eq 0 then begin
  writepol = [1,1]
endif else if chan eq 1 or chan eq 2 then begin
  writepol = (chan eq 1) ? [1,0] : [0,1]
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse

; Before opening the output file, make sure that at least
; one spectrum will be written (not a foregone conclusion
; if the chan keyword is set)

n = 0L
count = 0L

while n lt nstack1 and count eq 0 do begin

  ; See if stack1 spectrum #(n+1) will be written

  group_el = where( (groupvec eq (*stack1[n]).group) and                 $
                    ( (writepol[0] and (*stack1[n]).pol eq 'OC') or      $
                      (writepol[1] and (*stack1[n]).pol eq 'SC')     ),  $
                    count)
  n = n + 1
endwhile

if count eq 0 then begin
  print,'writestackpds1: No spectra meet these group/polarization criteria'
  return
endif

; Open the output file

err = openOutfile(lun,outfile,'csv',/get_lun,_extra=_ext)
if err ne 0 then return

; Go through the single-channel stack spectrum by spectrum and write out 
; spectra which are included in one of the specified groups

nwrite = 0L

for n=0L,nstack1-1 do begin

  ; See if this spectrum is included in one of the groups which will be written

  group_el = where(groupvec eq (*stack1[n]).group, count)

  if count gt 0 then begin

    ; This spectrum should be written to disk

    spec = (*stack1[n]).spec
    tags1 = (*stack1[n]).tags
    extratags = (*stack1[n]).extratags
    ndata = (*stack1[n]).ndata
    ntags = (*stack1[n]).ntags
    nextra = (*stack1[n]).nextra
    tname = (*stack1[n]).tname
    pol = (*stack1[n]).pol

    ; If requested, fix three tag values (xjcen, jsnr1, jsnr2) so that the output
    ; is consistent with tkplay, which counts these tags from 1 rather than from 0
    ;
    ; Also fix the jcp tag, which tkplay counts (sometimes) from 0 rather than from 1

    if keyword_set(tkplay) then begin
      tags1.xjcen = tags1.xjcen + 1L
      tags1.jsnr1 = tags1.jsnr1 + 1L
      tags1.jsnr2 = tags1.jsnr2 + 1L
      tags1.jcp = tags1.jcp - 1L
      if tags1.jcp lt 0 then begin
        print,'ERROR in writestackpds1: jcp tag for OC is already 0 rather than 1'
        return
      endif
      extratags = [extratags, extratags[0]]
      extratags[nextra].format = ''
      extratags[nextra].name = ''
      extratags[nextra].value = ''
      extratags[nextra].comment = $
                   '# Conv to   tkplay format ' + systime(/utc) + ' UT'
      nextra = nextra + 1L
    endif

    ; Write the pds header

    sformat = '(a,5x,a16,3x,a)'
    iformat = '(a,5x,a16,3x,i0)'
    printf,lun,'s','type','float',format=sformat
    printf,lun,'i','ndata',ndata,format=iformat
    printf,lun,'i','height',1,format=iformat
    printf,lun,'i','width',ndata+ntags,format=iformat
    printf,lun,'i','size',4,format=iformat
    printf,lun,'s','machine','SPARC',format=sformat
    printf,lun,'i','nchan',1,format=iformat
    printf,lun,'i','ntags',ntags,format=iformat
    printf,lun,'i','nspec',1,format=iformat
    printf,lun,'s','format','RDF_spectrum_1',format=sformat
    printf,lun,'.'

    ; Write the spectrum and tags (as binary)
    ; -- write all tags as floating point to be compatible with tkplay

    writeu,lun,swap_endian(spec, /swap_if_little_endian)
    for k=0L,ntags-1 do begin
      writeu,lun,swap_endian(float(tags1.(k)), /swap_if_little_endian)
    endfor

    ; Write the pds footer

    for k=0L,nextra-1 do begin
      printf,lun,strtrim(string(extratags[k],format='(a,5x,a16,3x,a,3x,a)'), 2)
    endfor
    printf,lun,'s','target',tname,format=sformat
    printf,lun,'s','polarization',pol,format=sformat
    printf,lun,'.'

    nwrite = nwrite + 1L

  endif
endfor

print,'Wrote ',nwrite,' pds frames, each containing 1 single-channel spectrum, to file ', $
      outfile,format='(a,i0,2a)'

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
