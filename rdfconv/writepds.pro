forward_function qq
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

if n_tags(pds) le 0 then begin
  print, "pds structure not set"
  return
endif

pairpointer = loaded

writeapdsfile,pairpointer,outfile,chan=chan,_extra=_ext

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setpds,show=show,imfile=infile,target=target,pname=pname,ttype=ttype,author=author,editor=editor,level=level,bookmark=bookmark,facet=facet,waves=waves,help=help,reset=reset

common pdsBlock,pds

if keyword_set(help) then begin
print, "setpds [,/show][,name='value']... [/help]"
print, "       /show shows the values"
print, "       /reset changes all to the default
print, "       sets one or more required values for PDS. Names can be:"
print, "       infile: OVERRIDE original rdf file in tags"
print, "       target: OVERRIDE target listed in tags."
print, "       pname: Product Name. REQUIRED, no default."
print, "       ttype: Target Type. Default: 'Asteroid'"
print, "       author: author XOR editor is REQUIRED"
print, "       editor: author XOR editor is REQUIRED"
print, "       level: Product Processing Level, default is 'Calibrated'"
print, "       facet: Science Search Facet, default='Tabulated,Physical Properties'"
print, "       version: Product version. Default is '1.0'"
print, "       waves: Wavelength range, default = 'Microwave'"
print, "       bookmark: override default bookmark string. Depending on xmit_sta
print, "                 tag, defaults to"
print, "                 'AO TX;AO RX;AO RI' or"
print, "                 'DSS14 TX;DSS14 RX'"
print, "       You can clear a value by setting it to ''"
endif
               
if (n_tags(pds) eq 0) or keyword_set(reset) then pds = {infile: '', target:'', pname:'',ttype:'Asteroid', author:'', editor:'', level:'Calibrated', bookmark:'',facet: 'Tabulated,Physical Properties', waves: 'Microwave', version: '1.0'}
if keyword_set(arecibo) then begin
  pds.bookmark = 'AO TX;AO RX;AO RI'
end
if keyword_set(dss14) then begin
  pds.bookmark = 'DSS14 TX;DSS14 RX'
endif
if defined(infile) then pds.infile = infile
if defined(target) then pds.target = target
if defined(pname) then pds.pname = pname
if defined(ttype) then pds.ttype = ttype
if defined(author) then pds.author = author
if defined(editor) then pds.editor = editor
if defined(level) then pds.level = level
if defined(facet) then pds.facet = facet
if defined(version) then pds.version = version
if defined(bookmark) then pds.bookmark = bookmark
if keyword_set(show) then help, pds,/str

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

if keyword_set(help) or isnull(name) or n_params() ne 2 then begin
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
pro writeapdsfile,pairpointer,outfile,tablelun,chan=chan,_extra=_ext
; This will actually do the writing for one file.
; First implement the writing of the data, TODO add a table

;
; PDS csv files. One pair per file. Starts with PDS keywords, then the column descriptions, then the data as three csv columns: frequency,pol1,pol2
; Should we allow one channel? I don't think it's needed, but should 
; almost work automagically.
; append will only mean anything if we make a table. Good idea?
; Need a way to specify bistatic rx station, nonstandard receiver or datataking setup.

common pdsBlock,pds

if n_params() lt 2 or keyword_set(help) then begin
  print, 'writepdsfile,pair,outfile,pds[,tablelun][,/overwrite][,chan=1 or 2][,/help]'
  return
endif

; pairpointter is a pointer to either the loaded pair or a stack entry. I *think* they look the same.
;pds is a structure containg the PDS keywords that aren't natively in the data

; Check which channel(s) to write

if n_elements(chan) eq 0 then begin
  writepol = [1,1]
  threecol = 1
  nchan = 2
endif else if chan eq 1 or chan eq 2 then begin
  writepol = (chan eq 1) ? [1,0] : [0,1]
  threecol = 0
  nchan = 1
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
addcr = string(13B)
if nchan eq 2 then begin
  extratags_write = extratags
  nextra_write = nextra
  addcomma = ',' + addcr
endif else begin
  splitExtra,chan,extratags,nextra,extratags_write,nextra_write
  addcomma = addcr
endelse

; Priority: set value, then rdf value
mytarget = tname
if not isnull(pds.target) then mytarget = pds.target
if isnull(mytarget) then begin
  print, "Target not set in pair or setpds"
  return
end
if isnull(pds.pname) then begin
  print, "Product name is required"
  return
end
if isnull(pds.author) + isnull(pds.editor) ne 1 then begin
  print, "Must have either an author or an editor"
  return
end
;
; Station-specific stuff. Need to mark better
;
xmitsta = getextral(extratags,'xmit_sta')
xmitpol = getextral(extratags,'xmit_poln')
if notnull(xmitpol) then polstring='xmit_poln,'+xmitpol
if xmitsta eq 'Arecibo' then begin
  bookmark='AO TX;AO RX;AO RI'
  if isnull(xmitpol) then polstring='xmit_poln,LCP'
endif else if xmitsta eq 'DSS14' then begin
  bookmark='DSS14 TX;DSS14 RX'
  if isnull(xmitpol) then polstring='xmit_poln,RCP'
endif
if notnull(pds.bookmark) then bookmark=pds.bookmark

; Get the times Assuming jdstart and end are correct for now.
; convert to ISO. Do start last so we have start times in vars
jdend = getextral(extratags,'jdend')
caldat, jdend, mon,day,year,hh,mm,ss
ss = fix(ss + 0.5) ; round to nearest to avoid jd rounding
endstring = string(year,mon,day,hh,mm,ss, format="(I04,'-',I02,'-',I02,'T',I02,':',I02,':',I02)")
jdstart = getextral(extratags,'jdstart')
caldat, jdstart, mon,day,year,hh,mm,ss
ss = fix(ss + 0.5) ; round to nearest to avoid jd rounding
startstring = string(year,mon,day,hh,mm,ss,format="(I04,'-',I02,'-',I02,'T',I02,':',I02,':',I02)") 

err = openOutfile(lun,outfile,'csv',/get_lun,_extra=_ext)
if err ne 0 then return

; Construct and print some PDS header strings
printf, lun, '# Keywords,',addcomma
printf, lun, 'Product Name,',qq(pds.pname),addcomma
printf, lun, 'Product Description,CW spectrum converted from RDF format',addcomma
printf, lun, 'Product Version,',qq(pds.version),addcomma
printf, lun, 'Start Time,', startstring,addcomma
printf, lun, 'Stop Time,', endstring, addcomma
printf, lun, 'Target Name,', qq(mytarget), addcomma
printf, lun, 'Target Type,', qq(pds.ttype), addcomma
if notnull(pds.author) then printf, lun, 'Author List,', qq(pds.author), addcomma
if notnull(pds.editor) then printf, lun, 'Editor List,', qq(pds.editor), addcomma
printf, lun, 'Product Processing Level,', qq(pds.level), addcomma
printf, lun, 'Science Search Facet,', qq(pds.facet), addcomma
printf, lun, 'Product Wavelength Ranges,',qq(pds.waves), addcomma
; next should be replaced with istruments and telescopes
printf, lun, 'Observing System Bookmark,', qq(bookmark), addcomma 
inf = getextral(extratags,'infile')
if notnull(pds.infile) then inf = pds.infile
if notnull(inf) then printf, lun, 'Original CW data file,', qq(inf), addcomma,format='(A,A,A)'
printf, lun, 'Software Version,20210411',addcomma

caldat, systime(/utc,/julian), mon,dd,yy, hh,mm,ss
ss = fix(ss + 0.5)
nowstring = string(yy,mon,dd,hh,mm,ss,format="(I04,'-',I02,'-',I02,'T',I02,':',I02,':',I02)")
printf, lun, 'Creation Date,', nowstring, addcomma

;
; Now tags
;
printf, lun, '# Tags,',addcomma
; Should allow override
skiptags=['iyy','imm','idd','rchour','rcmin','rcsec']
tagnames = strlowcase(tag_names(tags[0]))
for i = 0, ntags-1 do begin
dummy=where(strlowcase(tagnames[i]) eq skiptags, count)
  if count gt 0 then continue
  if threecol then begin
; can do either columns or _n
;    tn = string(tagnames[i], format='(a,"_1,")')
;    printf, lun, tn, tags[0].(i), addcomma
;    tn = string(tagnames[i], format='(a,"_2,")')
;    printf, lun, tn, tags[1].(i), addcomma
     o = string(tagnames[i],arf(tags[0].(i)),arf(tags[1].(i)),format='(A,",",A,",",A)')
  endif else begin
     o = string(tagnames[i],arf(tags[0].(i)),format='(A,",",A)')
  endelse
  printf, lun,o,addcomma
endfor

;
; and extra tags
;
printf, lun, '# Extra Tags,',addcomma
printf, lun, 'File date,', nowstring, addcomma

skipextra = ['xmit_pol','tzcorr','timezone']
for i = 0, nextra-1 do begin
  if (isnull(extratags[i].name) or extratags[i].format eq 't') then continue ; skip tag names
  dummy = where(strlowcase(extratags[i].name) eq skipextra, count)
  if count gt 0 then continue
  printf, lun, extratags[i].name, ',', qq(extratags[i].value), addcomma
endfor
;These were fixed up: keep in skiptags
printf, lun, polstring, addcomma
printf, lun, 'timezone,UTC',addcomma

;Column definitions

printf, lun, "# Column Definitions,", addcomma
if threecol then begin
  printf, lun, 'frequency,pol1,pol2',addcr
  printf, lun, 'real,real,real',addcr
  printf, lun, 'Hz,,',addcr
  printf, lun, ',,',addcr
  printf, lun, 'Offset from Ephemeris frequency,Polarization 1 normalized to unit standard deviation and zero mean in baseline,Polarization 2 normalized to unit standard deviation and zero mean in baseline',addcr
endif else begin
  printf, lun, 'frequency,pol1',addcr
  printf, lun, 'real,real',addcr
  printf, lun, 'Hz,',addcr
  printf, lun, 'Offset from Ephemeris frequency,Polarization 1 normalized to unit standard deviation and zero mean in baseline',addcr
endelse

;
; and now the data
;
printf, lun, '# Data,', addcomma

dfreq = tags[0].dfreq
posfr = tags[0].posfr
xjcen = tags[0].xjcen
freq = posfr*dfreq*(dindgen(ndata) - xjcen)


if threecol then begin
  for i = 0L, ndata-1 do begin
    printf, lun, freq[i],pair[0,i],pair[1,i],addcr, format='(e0,",",e0,",",e0,A1)'
  endfor
endif else begin
  for i = 0L, ndata-1 do begin
    printf, lun, freq[i],pair[0,i],addcr,format='(e0,",",e0,a1)'
  endfor
endelse

close, lun
free_lun, lun
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
; one pds file per pair.  The keywords g, ming and maxg, or ag specify the groups;
; omitting all of these causes all stack pairs to be written.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common pdsBlock,pds

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
  print,'              [,chan=1 or 2][,/help]'
  print,' '
  print,"    'file number and .csv' output file extension is added"
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
  print,'              [,chan=1 or 2][,/help]'
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

if n_tags(pds) le 0 then begin
  print, "pds structure not set"
  return
endif

; Go through the stack pair by pair and write out pairs which
; are included in one of the specified groups

filenum=1

sformat = '(a,5x,a16,3x,a)'
iformat = '(a,5x,a16,3x,i0)'
nwrite = 0L

for n=0L,nstack-1 do begin

  ; See if this pair is included in one of the groups which will be written

  group_el = where(groupvec eq (*stack[n]).group, count)

  if count gt 0 then begin

    ; This pair should be written to disk

    outname = string(outfile,filenum,format="(A,'_',I03)")
    pairpointer = stack[n]
    writeapdsfile,pairpointer,outname,chan=chan,_extra=_ext
    filenum = filenum + 1
    
    nwrite = nwrite + 1L

  endif
endfor

if nchan eq 2 then begin
  print,'Wrote ',nwrite,' pds files, each containing 1 OC and 1 SC spectrum, ', $
        format='(a,i0,a)'
endif else begin
  print,'Wrote ',nwrite,' pds files, each containing 1 ',chanstring[chan-1], $
        ' spectrum', format='(a,i0,3a)'
endelse

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
function qq, str
;
; Protect string with quotation marks if it contains a comma and the
; other kind of quotation mark if it contains a quotation mark,
;

hasd=strmatch(str, '*"*')
hasc=strmatch(str, '*,*')
hass=strmatch(str, "*'*")

if hasd+hass eq 2 then begin
  message, "Your string, "+str+ " has both kinds of quotes. I don't know what to do."
end

if hass+hasd+hasc eq 0 then return, str

if hasd then return, "'" + str + "'" else return, '"' + str + '"'

end

function arf, number
;
;format reals the way I want them.
;integers should be spelled out up to more digits than G allows.
if ((number eq long64(number)) and (number lt 1.e16)) then begin
  text = string(number, format='(i0)')
endif else begin
  text = string(number, format='(g0)')
end
if not strmatch(text,'*.*') then text = text + '.'
text = strcompress(text, /rem)
return, text
end
