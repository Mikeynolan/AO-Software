;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro readcsv,infile,silent=silent, $
            help=help

; Read in the entire contents of an csv-format disk file.
;
; Read in a single pair from a csv-format file in archive format.
;
;
; Aug 2002: If readrdf doesn't find the input file as specified,
;           add an '.csv' extension (if it's not already there) and try again
;
common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack


if n_params() ne 1L or not goodaddpol or keyword_set(help) then begin
  print,' '
  print,"readrdf,infile[,/nopairs][,/loose][,addpol='OC' or 'SC'][,/tkplay][,/silent][,/help]"
  print,' '
  print,"        The default input file extension is '.rdf'"
  print,' '
  print,'        /nopairs forces treatment as single-channel spectra'
  print,'        /loose   accepts spectra with incorrect jcp tags so long as'
  print,"                     nchan=1 and polarization='OC' or 'SC';"
  print,'                 accepts spectra that claim nchan=2 for a frame'
  print,"                     despite the frame having polarization='OC' or 'SC';"
  print,"                 automatically resets huge values of the 'igw' tag to zero"
  print,'        addpol   lets you set the polarization state of an image'
  print,'                     in case that tag is missing from the rdf file'
  print,'        /tkplay  imports CW spectra from tkplay and subtracts 1 from tags'
  print,'                     that are counted from 1 by tkplay but from 0 here'
  print,'        /silent  turns off screen output for individual input frames'
  print,' '
  return
endif else if size(infile, /type) ne 7 then begin
  print,' '
  print,"readrdf,infile[,/nopairs][,/loose][,addpol='OC' or 'SC'][,/tkplay][,/silent][,/help]"
  print,'Make sure that infile is a quoted string!'
  print,' '
  return
endif

; Open the input file; if that doesn't work, try it with an '.rdf' extension

err = openInfile(lun,infile,'csv',/get_lun)
if err ne 0 then return

; Make sure that the input file isn't empty

if eof(lun) then begin
  print,' '
  print,'ERROR in readrdf: ',infile,' is empty'
  print,' '
  close,lun
  free_lun,lun
  return
endif

; Get new group numbers for placing these spectra into the
; single-channel stack and/or images into the image stack.
; (We'll change the CW group number later if it turns out
; that all spectra are components of OC/SC pairs destined
; for the pair stack.)

groupCW = (nstack1 gt 0) ? max(stackgroups1()) + 1L : 1L
groupImage = (nstacki gt 0) ? max(stackgroupsi()) + 1L : 1L

; Get the tag names

tagnames = strlowcase(tag_names(blanktags1()))
tags = blanktags()

; Read an rdf frame each time through the loop

headline = ''   ; initialize as a string or else readf will assume float
tailline = ''   ; same here
userchoice = '' ; ditto (for read)
nframes = 0L
nreadCW = 0L
nreadi = 0L

  ; Initialize some parameters

  dataformat = ''
  ndata = 0L
  height = 0L
  width = 0L
  bytesPerValue = 0L
  polarization = ''
  targetname = '<unspecified>'
  machine = ''
  nchan = 2L
  ntags = 0L
  nspec = 1L
  fileformat = ''
  extratag = {format:'', name:'', value:'', comment:''}
  nextra = 0L
  has_spb_tag = 0
  has_infile_tag = 0
  infile_fullpath = (file_search(infile,/fully_qualify_path))[0]

  initial = 1
  intags=0
  inextra=0
  indata=0

; Initialize extra tags with tag names

nextra = ntags
extratags = replicate(extratag, ntags)
extratags.format = replicate('t', ntags)
extratags.name = tagnames
extratags.value = string(indgen(ntags), format='(i0)')
extratags.comment = ''


  ; Start reading in the header

  if not eof(lun) then begin
    readf,lun,headline,format='(a)'
  endif else begin
    print,' '
    print,'ERROR in readcsv: csv header is truncated'
    print,' '
    close,lun
    free_lun,lun
    return
  endelse

  ; Check that this isn't a FITS file

  fits_keyname = strtrim(strcompress(strmid(headline,0,8)), 2)
  if fits_keyname eq 'SIMPLE' then begin
    print,' '
    print,'ERROR in readcsv: This is a FITS file, use readfits to read it'
    print,' '
    close,lun
    free_lun,lun
    return
  endif

  if stregex(headline, "^ *# *keywords", /foldcase) lt 0 then begin
    print, ' '
    print, "ERROR in readcsv: This is not an archive-format CSV file"
    print, " "
    free_lun, lun
    return
  endif

; Use read_csv to take care of dealing with quoted commas, etc, then parse lines.
; If more than 4 columns, they we're dealing with future stokes data, so it's OK
; if they are numbers after that
; don't need the lun anymore.

free_lun, lun

  csv=read_csv(infile, count=count, types=replicate('string',4))

; First check for markers
    
  markers=['^ *tags','^ *extra *tags','^ *# *column *definition','^ *# *data']
  wmark = lonarr(n_elements(markers))
  for i = 0, n_elements(markers)-1 do begin
   wmark[i] = where(stregex(csv.field1,markers[i], /bool, /fold) gt 0)
  endfor

  if total(wmark) ne 4 then begin
    print, ' '
    print, 'ERROR in csvread: Could find marker strings. Looking for 4, got', total(wmark)
    print, ' '
    return
  endif

; Get start and stop time from top info
  ind = where(stregex(csv.field1, 'start *time', /fold, /bool));
  if ind[0] lt 0 then begin
    print, ' '
    print, 'ERROR in readcsv: Couldn''nt find start time'
    print, ' '
    return
  end
  starttime = csv.field2[ind[0]]
  ind = where(stregex(csv.field1, 'stop *time', /fold, /bool));
  if ind[0] lt 0 then begin
    print, ' '
    print, 'ERROR in readcsv: Couldn''nt find stop time'
    print, ' '
    return
  end
  stoptime = csv.field2[ind[0]]
  timestamptovalues, starttime, day=day,hour=hour,minute=minute,month=month,offset=offset,second=second,year=year
  tags[*].iyy = year
  tags[*].imm = month
  tags[*].idd = day
  tags[*].rchour = hour
  tags[*].rcmin = minute
  tags[*].rcsec = second

  ind = where(stregex(csv.field1, 'target *name', /fold, /bool)
  if ind[0] lt 0 then tname = 'Unknown' else tname = csv.field2[ind[0]]
  
  ; now tags.

  for itag = 0, n_elements(tagnames) -1 do begin
    row = where(strmatch(csv.field1[wmark[0]:wmark[1], tagnames[i], /bool, /fold))
    if row[0] lt 0 then continue
    if n_elements(row) ne 1 then begin
      print, ' '
      print, 'ERROR: found duplicate tags in csv file'
      print, ' '
      return
    endif
    tags[0].(itag) = double(csv.field2[row[0]])
    tags[1].(itag) = double(csv.field3[row[0]])
  endfor

  ; extra tags. Go ahead and put the description in the rdf file, why not
  ; Check for 'target': We'll add it in if it's not already there.
  for itag = wmark[1] + 1, wmark[2] - 1
    extratag.name = csv.field1[itag]
    extratag.value = csv.field2[itag]
    extratag.format = csv.field3[itag]
    extratag.comment = csv.field4[itag]
    extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
    nextra = nextra + 1
  endfor
  
; And the data.
  oc = float(csv.field2[wmark[3]+1:*])
  sc = float(csv.field3[wmark[3]+1:*])
  ndata = n_elements(oc)
  spec = transpose([[oc],[sc]]

  stackStruct = {group:groupCW, spec:spec, tags:tags, extratags:extratags, ndata:ndata, ntags:ntags, nextra:nextra, tname=targetname


  ; Check that the polarization tag(s) are OK;
  ; if so, set the pol_in array elements to 'OC' and/or 'SC'.


  if isCW then begin
    stack1Struc = {group:groupCW, spec:fltarr(ndata), tags:blanktags1(), $
                   extratags:extratags, ndata:ndata, ntags:ntags, $
                   nextra:nextra, tname:targetname, pol:''}
    if nreadCW eq 0 then $
      stackRead1 = ptrarr(height, /allocate_heap) $
    else $
      stackRead1 = [stackRead1, ptrarr(height, /allocate_heap)]
    for n=0L,height-1 do begin
      stack1Struc.spec = spec_in[n,*]
      for k=0L,ntags-1 do stack1Struc.tags.(k) = tags_in[n].(k)
      stack1Struc.pol = pol_in[n]
      *stackRead1[nreadCW+n] = stack1Struc
    endfor
    nreadCW = nreadCW + height
    nframes = nframes + 1L
    if not silent then begin
      print,'Finished reading rdf frame with ',height,' spectra',format='(a,i0,a)'
    endif
  endif else begin
    stackiStruc = {group:groupImage, image:fltarr(width,height), $
                   extratags:extratags, width:width, height:height, $
                   nextra:nextra, tname:targetname, pol:''}
    if nreadi eq 0 then $
      stackReadi = ptrarr(1, /allocate_heap) $
    else $
      stackReadi = [stackReadi, ptrarr(1, /allocate_heap)]
    stackiStruc.image = image_in
    stackiStruc.pol = pol_in[0]
    *stackReadi[nreadi] = stackiStruc
    nreadi = nreadi + 1L
    nframes = nframes + 1L
    if not silent then begin
      print,'Finished reading one ',width,' x ',height,' rdf image', $
            format='(a,i0,a,i0,a)'
    endif
  endelse

; Check whether all CW spectra can be paired

if nreadCW eq 0 or nopairs or (2*(nreadCW/2L) ne nreadCW) then begin

  ; Treat as single spectra, not pairs:
  ; no spectra, or an odd number of spectra, or else the "don't pair" flag was set

  pairSpectra = 0

endif else begin

  ; See if the file is broken into halves: first all the spectra from one
  ; polarization channel, then all the corresponding spectra from the other

  nhalf = nreadCW/2L
  halves = 1
  n = 0L
  while (halves and n lt nhalf) do begin
    halves = ((*stackRead1[n]).pol eq (*stackRead1[0]).pol) and $
             ((*stackRead1[n+nhalf]).pol ne (*stackRead1[0]).pol) and $
             ((*stackRead1[n]).tname eq (*stackRead1[n+nhalf]).tname) and $
             ((*stackRead1[n]).ndata eq (*stackRead1[n+nhalf]).ndata) and $
             ((*stackRead1[n]).ntags eq (*stackRead1[n+nhalf]).ntags) and $
             ((*stackRead1[n]).tags.xjcen eq (*stackRead1[n+nhalf]).tags.xjcen) and $
             ((*stackRead1[n]).tags.dfreq eq (*stackRead1[n+nhalf]).tags.dfreq) and $
             ((*stackRead1[n]).tags.posfr eq (*stackRead1[n+nhalf]).tags.posfr) and $
             ((*stackRead1[n]).tags.kpts eq (*stackRead1[n+nhalf]).tags.kpts) and $
             ((*stackRead1[n]).tags.igw eq (*stackRead1[n+nhalf]).tags.igw) and $
             ((*stackRead1[n]).tags.nfreq eq (*stackRead1[n+nhalf]).tags.nfreq) and $
             ((*stackRead1[n]).tags.frstep eq (*stackRead1[n+nhalf]).tags.frstep) and $
             ((*stackRead1[n]).tags.freq1 eq (*stackRead1[n+nhalf]).tags.freq1)
    n = n + 1L
  endwhile

  ; If that didn't pan out, see if the file consists entirely of alternating
  ; polarization channels: each spectrum immediately preceded or followed by the
  ; corresponding spectrum from the opposite channel

  alternating = not halves
  n = 0L
  while (alternating and n lt nreadCW-1) do begin
    alternating = ((*stackRead1[n]).pol ne (*stackRead1[n+1]).pol) and $
                  ((*stackRead1[n]).tname eq (*stackRead1[n+1]).tname) and $
                  ((*stackRead1[n]).ndata eq (*stackRead1[n+1]).ndata) and $
                  ((*stackRead1[n]).ntags eq (*stackRead1[n+1]).ntags) and $
                  ((*stackRead1[n]).tags.xjcen eq (*stackRead1[n+1]).tags.xjcen) and $
                  ((*stackRead1[n]).tags.dfreq eq (*stackRead1[n+1]).tags.dfreq) and $
                  ((*stackRead1[n]).tags.posfr eq (*stackRead1[n+1]).tags.posfr) and $
                  ((*stackRead1[n]).tags.kpts eq (*stackRead1[n+1]).tags.kpts) and $
                  ((*stackRead1[n]).tags.igw eq (*stackRead1[n+1]).tags.igw) and $
                  ((*stackRead1[n]).tags.nfreq eq (*stackRead1[n+1]).tags.nfreq) and $
                  ((*stackRead1[n]).tags.frstep eq (*stackRead1[n+1]).tags.frstep) and $
                  ((*stackRead1[n]).tags.freq1 eq (*stackRead1[n+1]).tags.freq1)
    n = n + 2L
  endwhile

  pairSpectra = halves or alternating

endelse

; Add the spectra to the appropriate stack
; (stack for pairs, stack1 for single-channel spectra)

if pairSpectra then begin

  ; Pairs: Get a new group number for adding them to the pair stack

  groupCW = (nstack gt 0) ? max(stackgroups()) + 1L : 1L
  
  ; Pair up the spectra just read in

  npairs = nhalf
  stackRead = ptrarr(npairs, /allocate_heap)
  for n=0L,npairs-1 do begin
    if halves then begin
      if (*stackRead1[n]).pol eq 'OC' then begin
        n1 = n
        n2 = n + nhalf
      endif else begin
        n1 = n + nhalf
        n2 = n
      endelse
    endif else begin
      if (*stackRead1[2*n]).pol eq 'OC' then begin
        n1 = 2*n
        n2 = 2*n + 1
      endif else begin
        n1 = 2*n + 1
        n2 = 2*n
      endelse
    endelse
    pair_ntags =  (*stackRead1[n1]).ntags
    pair_tags = blanktags()
    for k=0L,pair_ntags-1 do begin
      pair_tags.(k) = [(*stackRead1[n1]).tags.(k), (*stackRead1[n2]).tags.(k)]
    endfor
    mergeExtra,(*stackRead1[n1]).extratags,(*stackRead1[n2]).extratags, $
               (*stackRead1[n1]).nextra,(*stackRead1[n2]).nextra, $
               pair_extratags,pair_nextra
    stackStruc = {group:groupCW, $
                  spec:fltarr(2,(*stackRead1[n1]).ndata), $
                  tags:pair_tags, $
                  extratags:pair_extratags, $
                  ndata:(*stackRead1[n1]).ndata, $
                  ntags:pair_ntags, $
                  nextra:pair_nextra, $
                  tname:(*stackRead1[n1]).tname}
    stackStruc.spec[0,*] = (*stackRead1[n1]).spec
    stackStruc.spec[1,*] = (*stackRead1[n2]).spec
    *stackRead[n] = stackStruc
  endfor

  ; Add the pairs to the stack

  if nstack eq 0 then $
    stack = ptrarr(npairs, /allocate_heap) $
  else $
    stack = [stack, ptrarr(npairs, /allocate_heap)]
  for n=0L,npairs-1 do *stack[nstack+n] = *stackRead[n]
  nstack = nstack + npairs
  print,'Added ',npairs,' pairs to the pair stack, which now contains ', $
        nstack,' pairs',format='(a,i0,a,i0,a)'

endif else if nreadCW gt 0 then begin

  ; Single-channel spectra: Add them to stack1

  if nstack1 eq 0 then $
    stack1 = ptrarr(nreadCW, /allocate_heap) $
  else $
    stack1 = [stack1, ptrarr(nreadCW, /allocate_heap)]
  for n=0L,nreadCW-1 do *stack1[nstack1+n] = *stackRead1[n]
  nstack1 = nstack1 + nreadCW
  print,'Added ',nreadCW,' spectra to the single-channel stack, ', $
        'which now contains ',nstack1,' spectra',format='(a,i0,2a,i0,a)'

endif

; Add any images to the image stack, while also adding some
; Julian date extra tags to each image (if not already present)

if nreadi gt 0 then begin
  if nstacki eq 0 then $
    stacki = ptrarr(nreadi, /allocate_heap) $
  else $
    stacki = [stacki, ptrarr(nreadi, /allocate_heap)]
  for n=0L,nreadi-1 do begin
    *stacki[nstacki] = *stackReadi[n]
    nstacki = nstacki + 1L
    if isnull(getextrai('jdstart',stack=nstacki)) or $
             isnull(getextrai('jdend',stack=nstacki)) then begin
      year = getextrai('year',stack=nstacki)
      doy = getextrai('iday',stack=nstacki)
      rxstart_secs = getextrai('dtime',stack=nstacki)
      tau = getextrai('int_time',stack=nstacki)
      jdstart = (julday(1L,1L,year) - 0.5D0) + (doy - 1L) + rxstart_secs/86400D0
      jdend = jdstart + tau/86400D0
      setextrai,'d','jdstart',string(jdstart, format='(d13.5)'), $
                comment='RX start Julian date',stack=nstacki
      setextrai,'d','jdend',string(jdend, format='(d13.5)'), $
                comment='RX end Julian date',stack=nstacki
    endif
  endfor
  print,'Added ',nreadi,' images to the image stack, ', $
        'which now contains ',nstacki,' images',format='(a,i0,2a,i0,a)'
endif

; Clean up pointers no longer needed

if nreadCW gt 0 then ptr_free,stackRead1
if nreadi gt 0 then ptr_free,stackReadi
if pairSpectra then ptr_free,stackRead

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cw2dat,outfile,sc=sc,stack=n,help=help,_extra=_ext

; Write one channel of spectral data to a dat (ASCII text) file:
; -- Do this for the loaded pair, or else for stack pair n
;         if the stack keyword is set
; -- Do this for the OC spectrum, unless the /sc flag is set

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,'cw2dat,outfile[,/overwrite][,/append][,stack=n][,/sc][,/help]'
  print,' '
  print,'Write one channel of spectral data to a dat (ASCII text) file'
  print,' '
  print,'   Do this for the loaded pair, or else for stack pair n if the stack keyword is used'
  print,' '
  print,'   Do this for the OC spectrum, unless /sc is set'
  print,' '
  print,"   '.dat' output file extension is added if not already present"
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'cw2dat,outfile[,/overwrite][,/append][,stack=n][,/sc][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in cw2dat: The pair stack is empty, so nothing can be written to disk'
    return
  endif else if n le 0 then begin
    print,'ERROR in cw2dat: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in cw2dat: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in cw2dat: No pair is loaded, so nothing can be written to disk'
  return
endif

; Get the spectrum to be written to disk

if n_elements(n) gt 0 then begin
  pair = (*stack[n-1]).spec
endif else begin
  pair = (*loaded).spec
endelse
if keyword_set(sc) then begin
  spec = reform(pair[1,*])
  chanstring = 'SC'
endif else begin
  spec = reform(pair[0,*])
  chanstring = 'OC'
endelse

; Open the output file

err = openOutfile(lun,outfile,'dat',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the data and close up shop

printf,lun,spec,format='(e13.5)'

print,'Wrote 1 ',chanstring,' spectrum to file ',outfile

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cw2dat1,outfile,stack=n,help=help,_extra=_ext

; Write the spectral data from a single-channel spectrum to a dat (ASCII text) file:
; -- Do this for the loaded single-channel spectrum, or else for
;         stack1 spectrum n if the stack keyword is set

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,'cw2dat1,outfile[,/overwrite][,/append][,stack=n][,/help]'
  print,' '
  print,'Write the spectral data from a single-channel spectrum to a dat (ASCII text) file'
  print,' '
  print,'   Do this for the loaded single-channel spectrum, or else for'
  print,'        stack1 spectrum n if the stack keyword is used'
  print,' '
  print,"   '.dat' output file extension is added if not already present"
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'cw2dat1,outfile[,/overwrite][,/append][,stack=n][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in cw2dat1: The single-channel stack is empty, so nothing can be written to disk'
    return
  endif else if n le 0 then begin
    print,'ERROR in cw2dat1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in cw2dat1: There are only ',nstack1,' single-channel spectra in stack1', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in cw2dat1: No single-channel spectrum is loaded, so nothing can be written to disk'
  return
endif

; Get the spectrum to be written to disk

if n_elements(n) gt 0 then begin
  spec = (*stack1[n-1]).spec
  chanstring = (*stack1[n-1]).pol
endif else begin
  spec = (*loaded1).spec
  chanstring = (*loaded1).pol
endelse

; Open the output file

err = openOutfile(lun,outfile,'dat',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the data and close up shop

printf,lun,spec,format='(e13.5)'

print,'Wrote 1 ',chanstring,' spectrum to file ',outfile

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro im2dat,outfile,stack=n,help=help,_extra=_ext

; Write the data from an image to a dat (ASCII text) file:
; -- Do this for the loaded image, or else for
;         stacki image n if the stack keyword is set

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,'im2dat,outfile[,/overwrite][,/append][,stack=n][,/help]'
  print,' '
  print,'Write the data from an image to a dat (ASCII text) file'
  print,' '
  print,'   Do this for the loaded image, or else for'
  print,'        stacki image n if the stack keyword is used'
  print,' '
  print,"   '.dat' output file extension is added if not already present"
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'im2dat,outfile[,/overwrite][,/append][,stack=n][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in im2dat: The image stack is empty, so nothing can be written to disk'
    return
  endif else if n le 0 then begin
    print,'ERROR in im2dat: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in im2dat: There are only ',nstacki,' images in stacki', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in im2dat: No image is loaded, so nothing can be written to disk'
  return
endif

; Get the image to be written to disk

if n_elements(n) gt 0 then begin
  image = (*stacki[n-1]).image
  chanstring = (*stacki[n-1]).pol
endif else begin
  image = (*loadedi).image
  chanstring = (*loadedi).pol
endelse

; Open the output file

err = openOutfile(lun,outfile,'dat',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the data and close up shop

printf,lun,image,format='(e13.5)'

print,'Wrote 1 ',chanstring,' image to file ',outfile

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function fitsHeaderOK,dataformat,ndims,ndata,width,height,nspec,nchan,isCW

; Check that FITS header tags have reasonable values;
; print error messages if they don't
;
; This function is called by readfits

; Do the various checks

if isnull(dataformat) then begin
  print,'ERROR: BITPIX tag not included in FITS header'
  return, 0
endif else if ndims le 0 or ndims ge 3 then begin
  print,'ERROR: NAXIS = ',ndims,format='(a,i0)'
  print,'       readfits is only set up to handle NAXIS = 1 or 2'
  return, 0
endif

if isCW then begin
  if ndata le 0 then begin
    print,'ERROR: ndata = ',ndata,format='(a,i0)'
    print,'       readfits needs > 0 data points per spectrum'
    return, 0
  endif else if height le 0 then begin
    print,'ERROR: height = ',height,format='(a,i0)'
    print,'       readfits needs > 0 spectra per FITS frame'
    return, 0
  endif else if width ne ndata then begin
    print,'ERROR: width = ',width,format='(a,i0)'
    print,'       readfits expects width = ',ndata,format='(a,i0)'
    return, 0
  endif else if nspec ne height then begin
    print,'ERROR: nspec = ',nspec,format='(a,i0)'
    print,'       readfits expects nspec = height (= ',height,')',format='(a,i0,a)'
    return, 0
  endif else if (nchan ne 1 and nchan ne 2) then begin
    print,'ERROR: nchan = ',nchan,format='(a,i0)'
    print,'       readfits is only set up to handle 1 or 2 channels per FITS frame'
    return, 0
  endif
endif else begin
  if height le 0 then begin
    print,'ERROR: height = ',height,format='(a,i0)'
    print,'       readfits needs > 0 range bins (rows) per FITS image'
    return, 0
  endif else if width le 0 then begin
    print,'ERROR: width = ',width,format='(a,i0)'
    print,'       readfits needs > 0 frequency bins (columns) per FITS image'
    return, 0
  endif
endelse

return, 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro readfits,infile,nopairs=nopairs,addpol=addpol,silent=silent,help=help

; Read in the entire contents of an FITS-format disk file.
;
; If readfits doesn't find the input file as specified,
; add a '.fit' extension (if it's not already there) and try again;
; if that doesn't work, try it with a '.fits' extension
;
; As of June 2005 the section on reading spectra shouldn't be
; taken too seriously, as it's not yet clear how FITS-format
; radar spectra will be written (e.g., whether both spectra of a
; dual-polarization pair will be allowed in a single file)
;
; Jul 2006: Add the addpol keyword, which allows you to set an image's
;           polarization to 'OC' or 'SC' if that tag is missing in the FITS file

common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if not keyword_set(nopairs) then nopairs = 0
if not keyword_set(silent) then silent = 0

if n_elements(addpol) eq 0 then begin
  goodaddpol = 1
  addpol = ''
endif else if n_elements(addpol) gt 1 then begin
  goodaddpol = 0
endif else begin
  goodaddpol = addpol eq 'OC' or addpol eq 'oc' or addpol eq 'SC' or addpol eq 'sc'
endelse

if n_params() ne 1L or not goodaddpol or keyword_set(help) then begin
  print,' '
  print,"readfits,infile[,/nopairs][,addpol='OC or 'SC'][,/silent][,/help]"
  print,' '
  print,"         The default input file extension is '.fit'"
  print,"                   or (if that doesn't work) '.fits'"
  print,' '
  print,'         /nopairs forces treatment as single-channel spectra'
  print,'         addpol   lets you set the polarization state of an image'
  print,'                  in case that tag is missing from the rdf file'
  print,'         /silent  turns off screen output for individual input frames'
  print,' '
  return
endif else if size(infile, /type) ne 7 then begin
  print,' '
  print,"readfits,infile[,/nopairs][,addpol='OC or 'SC'][,/silent][,/help]"
  print,'Make sure that infile is a quoted string!'
  print,' '
  return
endif

; Open the input file; if that doesn't work, try it with a '.fit' extension
; and then with a '.fits' extension

err = openInfile(lun,infile,['fit','fits'],/get_lun)

; Make sure that the input file isn't empty

if eof(lun) then begin
  print,' '
  print,'ERROR in readfits: ',infile,' is empty'
  print,' '
  close,lun
  free_lun,lun
  return
endif

; Get new group numbers for placing these spectra into the
; single-channel stack and/or images into the image stack.
; (We'll change the CW group number later if it turns out
; that all spectra are components of OC/SC pairs destined
; for the pair stack.)

groupCW = (nstack1 gt 0) ? max(stackgroups1()) + 1L : 1L
groupImage = (nstacki gt 0) ? max(stackgroupsi()) + 1L : 1L

; Get the tag names in case there's a problem reading in cw tags

tagnames = strlowcase(tag_names(blanktags1()))

; Read an FITS frame each time through the loop

userchoice = ''   ; initialize as a string or else readf will assume float
year = ''         ;               "                reads        "
mon = ''
day = ''
hour = ''
min = ''
sec = ''
nframes = 0L
nreadCW = 0L
nreadi = 0L
headrecord = bytarr(80)
squote_byte = (byte("'"))[0]

while not eof(lun) do begin

  ; Initialize some parameters

  dataformat = ''
  isCW = 0L
  ndata = 0L
  height = 0L
  width = 0L
  bytesPerValue = 0L
  polarization = ''
  targetname = '<unspecified>'
  machine = ''
  nchan = 0L
  ntags = 0L
  nspec = 0L
  extratag = {format:'', name:'', value:'', comment:''}
  nextra = 0L
  has_infile_tag = 0
  has_eph_row_tag = 0
  has_eph_col_tag = 0
  infile_fullpath = (file_search(infile,/fully_qualify_path))[0]
  machine = 'sparc'
  n_header_bytes = 0L

  ; Start reading in the FITS header

  if not eof(lun) then begin
    readu,lun,headrecord
    n_header_bytes = n_header_bytes + 80L
    headline = string(headrecord)
    keyname = strtrim(strcompress(strmid(headline,0,8)), 2)
    if keyname ne 'SIMPLE' then begin
      print,' '
      print,"ERROR in readfits: FITS header should start with 'SIMPLE'"
      print,' '
      close,lun
      free_lun,lun
      return
    endif
  endif else begin
    print,' '
    print,'ERROR in readfits: FITS file is empty'
    print,' '
    close,lun
    free_lun,lun
    return
  endelse

  while keyname ne 'END' do begin

    ; Look out for blank records used to pad the end of the header block

    is_blank_padding = isnull(keyname)

    ; Start parsing this header record into name, value, comment
    ; by extracting the comment (if any)

    key_has_value = (strmid(headline,8,2) eq '= ')
    if key_has_value then begin

      ; Extract the comment (if it's there) from a line which has
      ; a key name/value pair.  The comment is separated from the
      ; key value by a forward slash ; don't be fooled by a forward
      ; slash which is part of a string key value (i.e, which is
      ; contained between single quotes)

      where_is_squote = where(headrecord eq squote_byte, n_squote)
      n_squote_skip = 2L*(n_squote/2)
      if n_squote_skip eq 0 then begin
        commentstart = strpos(headline,'/')
      endif else begin
        searchstart = where_is_squote[n_squote_skip-1] + 1L
        if searchstart eq strlen(headline) then begin
          commentstart = -1L
        endif else begin
          slashpos = strpos(strmid(headline,searchstart),'/')
          commentstart = (slashpos gt -1) ? (searchstart + slashpos) : -1L
        endelse
      endelse
      if commentstart eq -1 then begin
        keycomment = ''
        keyvaluestring = strtrim(strcompress(strmid(headline,10)), 2)
      endif else begin
        keycomment = strtrim(strmid(headline,commentstart+1), 2)
        keyvaluestring = strtrim(strcompress(strmid(headline,10,commentstart-10)), 2)
      endelse
    endif else begin

      ; If the keyname is 'COMMENT' or 'HISTORY' then there's no key value:
      ; everything which follows the key name is a comment

      if keyname eq 'COMMENT' then begin
        keycomment = strtrim(strcompress(strmid(headline,7)), 2)
        keyvaluestring = ''
      endif else begin
        keycomment = strtrim(strcompress(headline), 2)
        keyvaluestring = ''
      endelse
    endelse
    if notnull(keycomment) then keycomment = '# ' + keycomment
    keyvaluelen = strlen(keyvaluestring)

    ; Now extract the key value (if it's there) and also
    ; convert the key name to its rdf equivalent if the two differ

    if key_has_value then begin

      ; Some keys determine the basic structure of the image
      ; rather than showing up as tags

      if keyname eq 'BITPIX' then begin
        bitsPerValue = long(keyvaluestring)
        case bitsPerValue of
            8: dataformat = 'byte'
           16: dataformat = 'short'
           32: dataformat = 'int'
          -32: dataformat = 'float'
          -64: dataformat = 'double'
        endcase
        bytesPerValue = abs(bitsPerValue)/8L
      endif else if keyname eq 'NAXIS' then begin
        ndims = long(keyvaluestring)
        isCW = (ndims eq 1)
        if ndims lt 1 or ndims gt 2 then begin
          print,' '
          print,'ERROR in readfits: must have NAXIS = 1 or 2'
          print,' '
          close,lun
          free_lun,lun
          return
        endif
      endif else if keyname eq 'NAXIS1' then begin
        if isCW then begin
          ndata = long(keyvaluestring)
        endif else begin
          width = long(keyvaluestring)
        endelse
      endif else if keyname eq 'NAXIS2' then begin
        height = long(keyvaluestring)
      endif else if keyname eq 'POL' then begin
        polarization = strupcase(strtrim(strmid(keyvaluestring,1,keyvaluelen-2), 2))
      endif else if keyname eq 'OBJECT' then begin
        targetname = strtrim(strmid(keyvaluestring,1,keyvaluelen-2), 2)
      endif else if keyname eq 'NCHAN' then begin
        nchan = long(keyvaluestring)
      endif else if keyname eq 'NTAGS' then begin
        ntags = long(keyvaluestring)
      endif else if keyname eq 'NSPEC' then begin
        nspec = long(keyvaluestring)
      endif else if keyname eq 'DATE-OBS' then begin

        ; From here on down we're dealing with tags
        ;
        ; This particular date/time tag is a concatenation of seven tags
        ; (year mon day hour min sec, with UT implied in the tag definition),
        ; so extract the seven values and add those seven tags in addition
        ; to the concatenated tag
        ;
        ; However, don't add the seven if they exist individually and have
        ; already been added

        dateandtime = strtrim(strmid(keyvaluestring,1,keyvaluelen-2), 2)
        reads,dateandtime,year,mon,day,hour,min,sec,format='(a4,5(1x,a2))'
        mon = string(long(mon), format='(i0)')
        day = string(long(day), format='(i0)')
        hour = string(long(hour), format='(i0)')
        min = string(long(min), format='(i0)')
        sec = string(long(sec), format='(i0)')
        extratag = {format:'s', name:'date-obs', value:dateandtime, comment:keycomment}
        extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
        nextra = nextra + 1L
        refers_to_rx_start = (keycomment eq '# UT date of RX start')
        k = where(extratags.name eq 'timezone', count)
        if count eq 0 then begin
          comment = (refers_to_rx_start) ? '# time zone for RX start date/time' : ''
          extratags = [extratags, {format:'s', name:'timezone', value:'UTC', comment:comment}]
          nextra = nextra + 1L
        endif
        k = where(extratags.name eq 'year', count)
        if count eq 0 then begin
          comment = (refers_to_rx_start) ? '# year of RX start date/time (UTC)' : '# (UTC)'
          extratags = [extratags, {format:'i', name:'year',   value:year, comment:comment}]
          nextra = nextra + 1L
        endif
        k = where(extratags.name eq 'month', count)
        if count eq 0 then begin
          comment = (refers_to_rx_start) ? '# month of RX start date/time (UTC)' : '# (UTC)'
          extratags = [extratags, {format:'i', name:'month',  value:mon,  comment:comment}]
          nextra = nextra + 1L
        endif
        k = where(extratags.name eq 'day', count)
        if count eq 0 then begin
          comment = (refers_to_rx_start) ? '# day of RX start date/time (UTC)' : '# (UTC)'
          extratags = [extratags, {format:'i', name:'day',    value:day,  comment:comment}]
          nextra = nextra + 1L
        endif
        k = where(extratags.name eq 'hour', count)
        if count eq 0 then begin
          comment = (refers_to_rx_start) ? '# hour of RX start date/time (UTC)' : '# (UTC)'
          extratags = [extratags, {format:'i', name:'hour',   value:hour, comment:comment}]
          nextra = nextra + 1L
        endif
        k = where(extratags.name eq 'minute', count)
        if count eq 0 then begin
          comment = (refers_to_rx_start) ? '# minute of RX start date/time' : ''
          extratags = [extratags, {format:'i', name:'minute', value:min,  comment:comment}]
          nextra = nextra + 1L
        endif
        k = where(extratags.name eq 'second', count)
        if count eq 0 then begin
          comment = (refers_to_rx_start) ? '# second of RX start date/time' : ''
          extratags = [extratags, {format:'i', name:'second', value:sec,  comment:comment}]
          nextra = nextra + 1L
        endif
      endif else if keyname eq 'TIMEZONE' then begin

        ; This and the next six tags are part of the 'DATE-OBS' tag (see above);
        ; if they already exist individually in the FITS header record, these records
        ; should supersede the ones constructed from the 'DATE-OBS' tag.

        if strmid(keyvaluestring,0,1) eq "'" and $
                      strmid(keyvaluestring,keyvaluelen-1,1) eq "'" then begin
          keyvalue_temp = strmid(keyvaluestring,1,keyvaluelen-2)
          if strlen(keyvalue_temp) eq 0 then begin
            keyvalue = '<null>'
          endif else begin
            keyvalue_temp = strtrim(strcompress(keyvalue_temp), 2)
            keyvalue = (strlen(keyvalue_temp) eq 0) ? '<blank>' : keyvalue_temp
          endelse
        endif
        k = where(extratags.name eq 'timezone', count)
        if count eq 0 then begin
          extratags = [extratags, {format:'s', name:'timezone', value:keyvalue, comment:keycomment}]
          nextra = nextra + 1L
        endif else begin
          extratags[k] = {format:'s', name:'timezone', value:keyvalue, comment:keycomment}
        endelse
      endif else if keyname eq 'YEAR' then begin
        k = where(extratags.name eq 'year', count)
        if count eq 0 then begin
          extratags = [extratags, {format:'i', name:'year', value:keyvaluestring, comment:keycomment}]
          nextra = nextra + 1L
        endif else begin
          extratags[k[0]] = {format:'i', name:'year', value:keyvaluestring, comment:keycomment}
        endelse
      endif else if keyname eq 'MONTH' then begin
        k = where(extratags.name eq 'month', count)
        if count eq 0 then begin
          extratags = [extratags, {format:'i', name:'month', value:keyvaluestring, comment:keycomment}]
          nextra = nextra + 1L
        endif else begin
          extratags[k[0]] = {format:'i', name:'month', value:keyvaluestring, comment:keycomment}
        endelse
      endif else if keyname eq 'DAY' then begin
        k = where(extratags.name eq 'day', count)
        if count eq 0 then begin
          extratags = [extratags, {format:'i', name:'day', value:keyvaluestring, comment:keycomment}]
          nextra = nextra + 1L
        endif else begin
          extratags[k[0]] = {format:'i', name:'day', value:keyvaluestring, comment:keycomment}
        endelse
      endif else if keyname eq 'HOUR' then begin
        k = where(extratags.name eq 'hour', count)
        if count eq 0 then begin
          extratags = [extratags, {format:'i', name:'hour', value:keyvaluestring, comment:keycomment}]
          nextra = nextra + 1L
        endif else begin
          extratags[k[0]] = {format:'i', name:'hour', value:keyvaluestring, comment:keycomment}
        endelse
      endif else if keyname eq 'MINUTE' then begin
        k = where(extratags.name eq 'minute', count)
        if count eq 0 then begin
          extratags = [extratags, {format:'i', name:'minute', value:keyvaluestring, comment:keycomment}]
          nextra = nextra + 1L
        endif else begin
          extratags[k[0]] = {format:'i', name:'minute', value:keyvaluestring, comment:keycomment}
        endelse
      endif else if keyname eq 'SECOND' then begin
        k = where(extratags.name eq 'second', count)
        if count eq 0 then begin
          extratags = [extratags, {format:'i', name:'second', value:keyvaluestring, comment:keycomment}]
          nextra = nextra + 1L
        endif else begin
          extratags[k[0]] = {format:'i', name:'second', value:keyvaluestring, comment:keycomment}
        endelse
      endif else begin
        if keyname eq 'JDSTART' or keyname eq 'JDMEAN' or keyname eq 'JDEND' then begin
          rdfname = strlowcase(keyname)
          keyformat = 'd'
          keyvalue = keyvaluestring
        endif else if keyname eq 'DOY' then begin
          rdfname = 'iday'
          keyformat = 'i'
          keyvalue = keyvaluestring
        endif else if keyname eq 'SFM' then begin
          rdfname = 'dtime'
          keyformat = 'i'
          keyvalue = keyvaluestring
        endif else if keyname eq 'CDELT1' then begin
          rdfname = 'fres'
          keyformat = 'f'
          keyvalue = keyvaluestring
        endif else if keyname eq 'CDELT2' then begin
          rdfname = 'delayunit'
          keyformat = 'f'
          keyvalue = keyvaluestring
        endif else if keyname eq 'EXPOSURE' then begin
          rdfname = 'int_time'
          keyformat = 'f'
          keyvalue = keyvaluestring
        endif else if keyname eq 'CODELEN' then begin
          rdfname = 'code_len'
          keyformat = 'i'
          keyvalue = keyvaluestring
        endif else if keyname eq 'BANDWIDT' then begin
          rdfname = 'bw'
          keyformat = 'f'
          keyvalue = keyvaluestring
        endif else if keyname eq 'FFTLEN' then begin
          rdfname = 'freqs'
          keyformat = 'i'
          keyvalue = keyvaluestring
        endif else if keyname eq 'LOOKS' then begin
          rdfname = 'nlooks'
          keyformat = 'i'
          keyvalue = keyvaluestring
        endif else if keyname eq 'RANGERES' then begin
          rdfname = 'rangeunit'
          keyformat = 'f'
          keyvalue = keyvaluestring
        endif else if keyname eq 'SPB' or keyname eq 'SAMPPBD' then begin
          rdfname = 'samples_per_baud'
          keyformat = 'i'
          keyvalue = keyvaluestring
        endif else if keyname eq 'RPB' or keyname eq 'ROWSPBD' then begin
          rdfname = 'rows_per_baud'
          keyformat = 'i'
          keyvalue = keyvaluestring
        endif else if keyname eq 'CODEPROC' then begin
          rdfname = 'codemethod'
          keyformat = 's'
          keyvalue = strlowcase(strtrim(strmid(keyvaluestring,1,keyvaluelen-2), 2))
        endif else if keyname eq 'EPHEMERI' then begin
          rdfname = 'ephemeris'
          keyformat = 's'
          keyvalue = strtrim(strmid(keyvaluestring,1,keyvaluelen-2), 2)
        endif else if keyname eq 'TXPOWER' then begin
          rdfname = 'tx_power'
          keyformat = 'f'
          keyvalue = keyvaluestring
        endif else if keyname eq 'TSYSSIGS' then begin
          rdfname = 'tsyspersig'
          keyformat = 'f'
          keyvalue = keyvaluestring
        endif else if keyname eq 'SIGCNTS' then begin
          rdfname = 'sigma'
          keyformat = 'f'
          keyvalue = keyvaluestring
        endif else if keyname eq 'INFILE' then begin
          rdfname = 'infile'
          keyformat = 's'
          keyvalue = infile_fullpath
          has_infile_tag = 1
        endif else if keyname eq 'EPH_ROW' then begin
          rdfname = 'eph_row'
          keyformat = 'f'
          keyvalue = keyvaluestring
          has_eph_row_tag = 1
        endif else if keyname eq 'EPH_COL' then begin
          rdfname = 'eph_col'
          keyformat = 'f'
          keyvalue = keyvaluestring
          has_eph_col_tag = 1
        endif else begin
          rdfname = strlowcase(keyname)
          if keyname eq 'LAMBDA' or keyname eq 'TXOFFSET' or keyname eq 'BAUDLEN' $
                                 or keyname eq 'LOOPTIME' or keyname eq 'DC_COL'  $
                                 or keyname eq 'CRPIX1'   or keyname eq 'CRPIX2'  $
                                 then begin
            keyformat = 'f'
            keyvalue = keyvaluestring
          endif else if keyname eq 'FREQSAMP' then begin
            keyformat = 'i'
            keyvalue = keyvaluestring
          endif else if keyvaluestring eq 'T' or keyvaluestring eq 'F' then begin

            ; From here on down we're dealing with tags for which
            ; we can't guess the name in advance

            keyformat = 's'
            keyvalue = keyvaluestring
          endif else if strmid(keyvaluestring,0,1) eq "'" and $
                        strmid(keyvaluestring,keyvaluelen-1,1) eq "'" then begin
            keyformat = 's'
            keyvalue_temp = strmid(keyvaluestring,1,keyvaluelen-2)
            if strlen(keyvalue_temp) eq 0 then begin
              keyvalue = '<null>'
            endif else begin
              keyvalue_temp = strtrim(strcompress(keyvalue_temp), 2)
              keyvalue = (strlen(keyvalue_temp) eq 0) ? '<blank>' : keyvalue_temp
            endelse
          endif else if strpos(keyvaluestring,'.') eq -1 then begin
            keyformat = 'i'
            keyvalue = keyvaluestring
          endif else begin
            keyformat = 'd'
            keyvalue = keyvaluestring
          endelse
        endelse

        extratag.format = keyformat
        extratag.name = strlowcase(rdfname)
        extratag.value = keyvalue
        extratag.comment = keycomment
        extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
        nextra = nextra + 1L
      endelse

    endif else if not is_blank_padding then begin

      ; Just a comment, no key value

      extratag.format = ''
      extratag.name = ''
      extratag.value = ''
      extratag.comment = keycomment
      extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
      nextra = nextra + 1L
    endif

    ; Read the next header line

    if not eof(lun) then begin
      readu,lun,headrecord
      n_header_bytes = n_header_bytes + 80L
      headline = string(headrecord)
      keyname = strtrim(strmid(headline,0,8), 2)
    endif else begin
      print,' '
      print,'ERROR in readfits: FITS header is truncated'
      print,' '
      close,lun
      free_lun,lun
      return
    endelse

  endwhile

  ; Check that the header is OK

  if ndims lt 1 or ndims gt 2 then begin
    print,' '
    print,'ERROR in readfits: must have NAXIS = 1 or 2'
    print,' '
    close,lun
    free_lun,lun
    return
  endif
  if not fitsHeaderOK(dataformat,ndims,ndata, $
                      width,height,nspec,nchan,isCW) then begin
    close,lun
    free_lun,lun
    return
  endif

  ; Get to the end of the header's last 2880-byte block

  npartial = n_header_bytes mod 2880
  nskip = (npartial gt 0) ? (2880 - npartial) : 0L
  if nskip gt 0 then begin
    skiparr = bytarr(nskip)
    readu,lun,skiparr
  endif

  ; Read in the data (and then skip to the end
  ; of the data's last 2880-byte block)

  if isCW then begin
    if dataformat eq 'byte' then begin
      spec1 = bytarr(ndata)
    endif else if dataformat eq 'short' then begin
      spec1 = intarr(ndata)
    endif else if dataformat eq 'int' then begin
      spec1 = lonarr(ndata)
    endif else if dataformat eq 'float' then begin
      spec1 = fltarr(ndata)
    endif else begin
      spec1 = dblarr(ndata)
    endelse
    npartial = ndata*bytesPerValue mod 2880
    nskip = (npartial gt 0) ? (2880 - npartial) : 0L
    if nskip gt 0 then skiparr = bytarr(nskip)
    spec_in = fltarr(height, ndata)
    fltTags = fltarr(ntags)
    tags_in = replicate(blanktags1(), height)
    for n=0L,height-1 do begin
      readu,lun,spec1
      if nskip gt 0 then readu,lun,skiparr
      spec1 = swap_endian(temporary(spec1), /swap_if_little_endian)
      spec_in[n,*] = (dataformat eq 'float') ? spec1 : float(spec1)
      for tagnum=0L,ntags-1 do begin

        ; Make sure there will be no integer overflow when converting cw tags

        if not tagOK(fltTags[tagnum],tagtypes[tagnum]) then begin
          print,' '
          print,'WARNING: cw tag value too large for integer conversion'
          print,'         Spectrum #',n+1,' in FITS frame #',nframes+1L,' has ', $
                     tagnames[tagnum],' = ',fltTags[tagnum],format='(a,i0,a,i0,3a,f0)'
          print,' '
          userprompt = '         Reset to zero and continue? (y/n) '
          read,userchoice,prompt=userprompt
          userchoice = strlowcase(strmid(strtrim(userchoice, 2), 0, 1))
          if userchoice eq 'y' or isnull(userchoice) then begin
            fltTags[tagnum] = 0
            print,'resetting tag value to zero'
          endif else begin
            close,lun
            free_lun,lun
            print,'readfits canceled'
            return
          endelse
        endif
        tags_in[n].(tagnum) = fix(fltTags[tagnum], type=tagtypes[tagnum], /print)
      endfor
    endfor
  endif else begin
    if dataformat eq 'byte' then begin
      image_in = bytarr(width, height)
    endif else if dataformat eq 'short' then begin
      image_in = intarr(width, height)
    endif else if dataformat eq 'int' then begin
      image_in = lonarr(width, height)
    endif else if dataformat eq 'float' then begin
      image_in = fltarr(width, height)
    endif else begin
      image_in = dblarr(width, height)
    endelse
    readu,lun,image_in
    npartial = width*height*bytesPerValue mod 2880
    nskip = (npartial gt 0) ? (2880 - npartial) : 0L
    if nskip gt 0 then begin
      skiparr = bytarr(nskip)
      readu,lun,skiparr
    endif
    image_in = swap_endian(temporary(image_in), /swap_if_little_endian)
    if (dataformat ne 'float') then image_in = float(image_in)
  endelse

  ; Check that the polarization tag(s) are OK;
  ; if so, set the pol_in array elements to 'OC' and/or 'SC'.

  jcp_in = (isCW) ? tags_in[*].jcp : [-1]
  if not polOK(isCW,height,nchan,polarization,jcp_in,pol_in,addpol=addpol) then begin
    close,lun
    free_lun,lun
    return
  endif
  if isCW then tags_in[*].jcp = jcp_in

  ; Add the 'infile' extra tag if necessary

  if not has_infile_tag then begin
    extratag.format = 's'
    extratag.name = 'infile'
    extratag.value = infile_fullpath
    extratag.comment = ''
    extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
    nextra = nextra + 1L
  endif

  ; For images, if the 'eph_row' and/or 'eph_col' extra tags are absent,
  ; get them (if possible) from the 'crpix2' and 'crpix1' tags, since
  ; crpix1 = eph_col + 1 and crpix2 = eph_row + 1

  if not isCW and not has_eph_row_tag then begin
    k = where(extratags.name eq 'crpix2', count)
    if count gt 0 then begin
      extratag.format = 'f'
      extratag.name = 'eph_row'
      extratag.value = extratags[k[0]].value - 1.0
      extratag.comment = 'row in which ephemeris range lies, 0-based'
      extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
      nextra = nextra + 1L
    endif
  endif

  if not isCW and not has_eph_col_tag then begin
    k = where(extratags.name eq 'crpix1', count)
    if count gt 0 then begin
      extratag.format = 'f'
      extratag.name = 'eph_col'
      extratag.value = extratags[k[0]].value - 1.0
      extratag.comment = 'column in which ephemeris Doppler lies, 0-based'
      extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
      nextra = nextra + 1L
    endif
  endif

  ; Make sure there's at least one extra tag (so that the extratags array is defined):
  ; Add the tag names (CW) or a dummy comment (images) if necessary.

  if nextra eq 0 then begin
    if isCW then begin
      nextra = ntags
      extratags = replicate(extratag, ntags)
      extratags.format = replicate('t', ntags)
      extratags.name = strlowcase(tag_names(tags_in))
      extratags.value = string(indgen(ntags), format='(i0)')
      extratags.comment = ''
    endif else begin
      nextra = 1L
      extratags = [extratag]
      extratags[0].format = ''
      extratags[0].name = ''
      extratags[0].value = ''
      extratags[0].comment = '# This is a dummy comment'
    endelse
  endif

  ; Turn each spectrum or image into a structure and add it to a temporary stack

  if isCW then begin
    stack1Struc = {group:groupCW, spec:fltarr(ndata), tags:blanktags1(), $
                   extratags:extratags, ndata:ndata, ntags:ntags, $
                   nextra:nextra, tname:targetname, pol:''}
    if nreadCW eq 0 then $
      stackRead1 = ptrarr(height, /allocate_heap) $
    else $
      stackRead1 = [stackRead1, ptrarr(height, /allocate_heap)]
    for n=0L,height-1 do begin
      stack1Struc.spec = spec_in[n,*]
      for k=0L,ntags-1 do stack1Struc.tags.(k) = tags_in[n].(k)
      stack1Struc.pol = pol_in[n]
      *stackRead1[nreadCW+n] = stack1Struc
    endfor
    nreadCW = nreadCW + height
    nframes = nframes + 1L
    if not silent then begin
      print,'Finished reading FITS frame with ',height,' spectra',format='(a,i0,a)'
    endif
  endif else begin
    stackiStruc = {group:groupImage, image:fltarr(width,height), $
                   extratags:extratags, width:width, height:height, $
                   nextra:nextra, tname:targetname, pol:''}
    if nreadi eq 0 then $
      stackReadi = ptrarr(1, /allocate_heap) $
    else $
      stackReadi = [stackReadi, ptrarr(1, /allocate_heap)]
    stackiStruc.image = image_in
    stackiStruc.pol = pol_in[0]
    *stackReadi[nreadi] = stackiStruc
    nreadi = nreadi + 1L
    nframes = nframes + 1L
    if not silent then begin
      print,'Finished reading one ',width,' x ',height,' FITS image', $
            format='(a,i0,a,i0,a)'
    endif
  endelse

  ; Back to the top of the main loop: Read in another FITS frame

endwhile

; All frames have been read in: close the input file

close,lun
free_lun,lun

; Check whether all CW spectra can be paired

if nreadCW eq 0 or nopairs or (2*(nreadCW/2L) ne nreadCW) then begin

  ; Treat as single spectra, not pairs:
  ; no spectra, or an odd number of spectra, or else the "don't pair" flag was set

  pairSpectra = 0

endif else begin

  ; See if the file is broken into halves: first all the spectra from one
  ; polarization channel, then all the corresponding spectra from the other

  nhalf = nreadCW/2L
  halves = 1
  n = 0L
  while (halves and n lt nhalf) do begin
    halves = ((*stackRead1[n]).pol eq (*stackRead1[0]).pol) and $
             ((*stackRead1[n+nhalf]).pol ne (*stackRead1[0]).pol) and $
             ((*stackRead1[n]).tname eq (*stackRead1[n+nhalf]).tname) and $
             ((*stackRead1[n]).ndata eq (*stackRead1[n+nhalf]).ndata) and $
             ((*stackRead1[n]).ntags eq (*stackRead1[n+nhalf]).ntags) and $
             ((*stackRead1[n]).tags.xjcen eq (*stackRead1[n+nhalf]).tags.xjcen) and $
             ((*stackRead1[n]).tags.dfreq eq (*stackRead1[n+nhalf]).tags.dfreq) and $
             ((*stackRead1[n]).tags.posfr eq (*stackRead1[n+nhalf]).tags.posfr) and $
             ((*stackRead1[n]).tags.kpts eq (*stackRead1[n+nhalf]).tags.kpts) and $
             ((*stackRead1[n]).tags.igw eq (*stackRead1[n+nhalf]).tags.igw) and $
             ((*stackRead1[n]).tags.nfreq eq (*stackRead1[n+nhalf]).tags.nfreq) and $
             ((*stackRead1[n]).tags.frstep eq (*stackRead1[n+nhalf]).tags.frstep) and $
             ((*stackRead1[n]).tags.freq1 eq (*stackRead1[n+nhalf]).tags.freq1)
    n = n + 1L
  endwhile

  ; If that didn't pan out, see if the file consists entirely of alternating
  ; polarization channels: each spectrum immediately preceded or followed by the
  ; corresponding spectrum from the opposite channel

  alternating = not halves
  n = 0L
  while (alternating and n lt nreadCW-1) do begin
    alternating = ((*stackRead1[n]).pol ne (*stackRead1[n+1]).pol) and $
                  ((*stackRead1[n]).tname eq (*stackRead1[n+1]).tname) and $
                  ((*stackRead1[n]).ndata eq (*stackRead1[n+1]).ndata) and $
                  ((*stackRead1[n]).ntags eq (*stackRead1[n+1]).ntags) and $
                  ((*stackRead1[n]).tags.xjcen eq (*stackRead1[n+1]).tags.xjcen) and $
                  ((*stackRead1[n]).tags.dfreq eq (*stackRead1[n+1]).tags.dfreq) and $
                  ((*stackRead1[n]).tags.posfr eq (*stackRead1[n+1]).tags.posfr) and $
                  ((*stackRead1[n]).tags.kpts eq (*stackRead1[n+1]).tags.kpts) and $
                  ((*stackRead1[n]).tags.igw eq (*stackRead1[n+1]).tags.igw) and $
                  ((*stackRead1[n]).tags.nfreq eq (*stackRead1[n+1]).tags.nfreq) and $
                  ((*stackRead1[n]).tags.frstep eq (*stackRead1[n+1]).tags.frstep) and $
                  ((*stackRead1[n]).tags.freq1 eq (*stackRead1[n+1]).tags.freq1)
    n = n + 2L
  endwhile

  pairSpectra = halves or alternating

endelse

; Add the spectra to the appropriate stack
; (stack for pairs, stack1 for single-channel spectra)

if pairSpectra then begin

  ; Pairs: Get a new group number for adding them to the pair stack

  groupCW = (nstack gt 0) ? max(stackgroups()) + 1L : 1L
  
  ; Pair up the spectra just read in

  npairs = nhalf
  stackRead = ptrarr(npairs, /allocate_heap)
  for n=0L,npairs-1 do begin
    if halves then begin
      if (*stackRead1[n]).pol eq 'OC' then begin
        n1 = n
        n2 = n + nhalf
      endif else begin
        n1 = n + nhalf
        n2 = n
      endelse
    endif else begin
      if (*stackRead1[2*n]).pol eq 'OC' then begin
        n1 = 2*n
        n2 = 2*n + 1
      endif else begin
        n1 = 2*n + 1
        n2 = 2*n
      endelse
    endelse
    pair_ntags =  (*stackRead1[n1]).ntags
    pair_tags = blanktags()
    for k=0L,pair_ntags-1 do begin
      pair_tags.(k) = [(*stackRead1[n1]).tags.(k), (*stackRead1[n2]).tags.(k)]
    endfor
    mergeExtra,(*stackRead1[n1]).extratags,(*stackRead1[n2]).extratags, $
               (*stackRead1[n1]).nextra,(*stackRead1[n2]).nextra, $
               pair_extratags,pair_nextra
    stackStruc = {group:groupCW, $
                  spec:fltarr(2,(*stackRead1[n1]).ndata), $
                  tags:pair_tags, $
                  extratags:pair_extratags, $
                  ndata:(*stackRead1[n1]).ndata, $
                  ntags:pair_ntags, $
                  nextra:pair_nextra, $
                  tname:(*stackRead1[n1]).tname}
    stackStruc.spec[0,*] = (*stackRead1[n1]).spec
    stackStruc.spec[1,*] = (*stackRead1[n2]).spec
    *stackRead[n] = stackStruc
  endfor

  ; Add the pairs to the stack

  if nstack eq 0 then $
    stack = ptrarr(npairs, /allocate_heap) $
  else $
    stack = [stack, ptrarr(npairs, /allocate_heap)]
  for n=0L,npairs-1 do *stack[nstack+n] = *stackRead[n]
  nstack = nstack + npairs
  print,'Added ',npairs,' pairs to the pair stack, which now contains ', $
        nstack,' pairs',format='(a,i0,a,i0,a)'

endif else if nreadCW gt 0 then begin

  ; Single-channel spectra: Add them to stack1

  if nstack1 eq 0 then $
    stack1 = ptrarr(nreadCW, /allocate_heap) $
  else $
    stack1 = [stack1, ptrarr(nreadCW, /allocate_heap)]
  for n=0L,nreadCW-1 do *stack1[nstack1+n] = *stackRead1[n]
  nstack1 = nstack1 + nreadCW
  print,'Added ',nreadCW,' spectra to the single-channel stack, ', $
        'which now contains ',nstack1,' spectra',format='(a,i0,2a,i0,a)'

endif

; Add any images to the image stack, while also adding some
; Julian date extra tags to each image (if not already present)

if nreadi gt 0 then begin
  if nstacki eq 0 then $
    stacki = ptrarr(nreadi, /allocate_heap) $
  else $
    stacki = [stacki, ptrarr(nreadi, /allocate_heap)]
  for n=0L,nreadi-1 do begin
    *stacki[nstacki] = *stackReadi[n]
    nstacki = nstacki + 1L
    if isnull(getextrai('jdstart',stack=nstacki)) or $
             isnull(getextrai('jdend',stack=nstacki)) then begin
      year = getextrai('year',stack=nstacki)
      doy = getextrai('iday',stack=nstacki)
      rxstart_secs = getextrai('dtime',stack=nstacki)
      tau = getextrai('int_time',stack=nstacki)
      if notnull(year) and notnull(doy) and notnull(rxstart_secs) $
                                        and notnull(tau) then begin
        jdstart = (julday(1L,1L,year) - 0.5D0) + (doy - 1L) + rxstart_secs/86400D0
        jdend = jdstart + tau/86400D0
        setextrai,'d','jdstart',string(jdstart, format='(d13.5)'), $
                  comment='RX start Julian date',stack=nstacki
        setextrai,'d','jdend',string(jdend, format='(d13.5)'), $
                  comment='RX end Julian date',stack=nstacki
      endif
    endif
  endfor
  print,'Added ',nreadi,' images to the image stack, ', $
        'which now contains ',nstacki,' images',format='(a,i0,2a,i0,a)'
endif

; Clean up pointers no longer needed

if nreadCW gt 0 then ptr_free,stackRead1
if nreadi gt 0 then ptr_free,stackReadi
if pairSpectra then ptr_free,stackRead

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro processtag,tagname,header,header_byte,usedtag,nrec,comment=comment,skip=skip, $
               _extra=_ext

; Look for a tag in the string-format FITS header; it if it's there, convert it to
; one or more properly formatted FITS records and add it/them to the byte-format header
;
; If /skip is set then just mark the tag as used without adding it to the header

if n_elements(comment) eq 0 then begin
  k = where(header[0,*] eq tagname)
endif else begin
  k = where(header[0,*] eq tagname and header[2,*] eq comment)
endelse
if k ne -1 then begin
  newrecord = fitsheadrecord(header[0,k], header[1,k], header[2,k], _extra=_ext)
  if notnull(newrecord) and not keyword_set(skip) then begin
    n_new = (size(newrecord))[2]
    header_byte = transpose( [ transpose(header_byte), transpose(newrecord) ] )
    nrec = nrec + n_new
  endif
  usedtag[k] = 1
endif
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writefitsi,outfile,help=help,_extra=_ext

; Write the loaded image to disk as an FITS file
;
; This routine is vastly longer and clunkier than it needs to be,
; all so that the header records can be output in a sensible order
; (i.e., the same order used by Mike Nolan's "writefits" Perl script)

common loadedBlock,loadedi,loaded1,loaded

; The following four arrays involve tags output by the "writefits" and "writerdf" Perl scripts
; -- more specifically, tags for which the names and/or comments differ depending on which
; output format (rdf vs. FITS) is used.

rdftagnames =  ['iday',             'dtime',         'fres',         'delayunit', 'int_time',  $
                'code_len',         'bw',            'freqs',        'nlooks',    'rangeunit', $
                'samples_per_baud', 'rows_per_baud', 'codemethod',   'ephemeris', 'tx_power',  $
                'tsyspersig',       'sigma',         'coherent_avg', 'codebw',    'baudlen',   $
                'lambda',           'txoffset',      'freqsamp',     'looptime',  'eph_row',   $
                'eph_col',          'dc_col',        'jdstart',      'jdmean',    'jdend',     $
                'data_type',        'elevation',     'gain_rxtx']
fitstagnames = ['DOY',              'SFM',           'CDELT1',       'CDELT2',    'EXPOSURE',  $
                'CODELEN',          'BANDWIDT',      'FFTLEN',       'LOOKS',     'RANGERES',  $
                'SAMPPBD',          'ROWSPBD',       'CODEPROC',     'EPHEMERI',  'TXPOWER',   $
                'TSYSSIGS',         'SIGCNTS',       'COHAVG',       'CODEBW',    'BAUDLEN',   $
                'LAMBDA',           'TXOFFSET',      'FREQSAMP',     'LOOPTIME',  'EPH_ROW',   $
                'EPH_COL',          'DC_COL',        'JDSTART',      'JDMEAN',    'JDEND',     $
                'DATATYPE',         'ELEV',          'GAINRXTX']
rdfcomments =  ['Day Of Year (UTC)', $
                'Receive start time (UTC seconds)', $
                '(Hz)', $
                'in microseconds / pixel', $
                'Exposure time included in image (s)', $
                '', $
                'Of original unvignetted image (Hz)', $
                'Transform length', $
                '', $
                'in meters / pixel', $
                '', $
                '', $
                '', $
                '', $
                'TX power in kW', $
                'Fraction of Tsys in 1 sigma from background subtraction', $
                'Raw counts in 1 sigma 1 sigma from background subtraction', $
                'number of codes coherently averaged before Doppler FFT', $
                'Bandwidth before coherent averaging (Hz)', $
                'in microseconds', $
                '(m)', $
                '(Hz)', $
                'If FREQSAMP < FREQS, the transform was zero-filled.', $
                'Closed-loop time included in EPH_ROW', $
                'row in which ephemeris range lies, 0-based', $
                'column in which ephemeris Doppler lies, 0-based', $
                'column in which zero Doppler (DC) lies, 0-based', $
                '', $
                '', $
                '', $
                '', $
                'mean elevation in deg', $
                'mean rx_gain * tx_gain']
fitscomments = ['RX start Day of Year (UT January 1 = 1)', $
                '[s] RX start seconds from UT midnight',   $
                '[Hz]', $
                '[us]', $
                '[s] Receive time processed in this image', $
                '', $
                '[Hz] Bandwidth of unvignetted image', $
                'FFT length used in the frequency processing', $
                'Incoherently summed spectra', $
                '[m] Range extent per pixel', $
                'Samples Per Baud', $
                'Rows per baud', $
                'Code processing method', $
                '', $
                'Fraction of Tsys in 1 sigma from background subtraction', $
                '[kW] Transmitted power', $
                'Raw counts in 1 sigma from background subtraction', $
                'Number of codes coherently avgd before FFT', $
                '[Hz] Bandwidth before coherent averaging', $
                '[us]', $
                '[m] wavelength', $
                '[Hz] TX offset', $
                'Input Frequency Samples', $
                '[us] Closed-loop time used in axis calculation', $
                'Row in which ephemeris range lies, zero-based', $
                'Column... ephemeris Doppler lies, zero-based', $
                'Column... zero Doppler (DC) lies, zero-based', $
                'RX start Julian date', $
                'mid-RX Julian date', $
                'RX end Julian date', $
                '', $
                '[deg] mean elevation', $
                'mean rx_gain * tx_gain']

; Give a help message if necessary

if n_params() ne 1 or keyword_set(help) then begin
  print,'writefitsi,outfile[,/overwrite][,/append][,/help]'
  print,"           '.fits' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writefitsi,outfile[,/overwrite][,/append][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in writefitsi: No image is loaded,', $
        ' so nothing can be written to disk',format='(2a)'
  return
endif

; Get some elements of the loaded image

image = (*loadedi).image
height = (*loadedi).height
width = (*loadedi).width
extratags = (*loadedi).extratags
nextra = (*loadedi).nextra
tname = (*loadedi).tname
pol = (*loadedi).pol

; Initialize the string-format FITS header with basic information
;
; Note that the ORDER of these records is pretty random (beyond the first eight records);
; we'll get them properly ordered later on when we switch to byte format

nhead = 21L
header = strarr(3,nhead)
header[*, 0] = ['SIMPLE', 'T', 'file does conform to FITS standard']
header[*, 1] = ['BITPIX', '-32', 'number of bits per data pixel']
header[*, 2] = ['NAXIS', '2', 'number of data axes']
header[*, 3] = ['NAXIS1', string(width, format='(i0)'), 'length of data axis 1']
header[*, 4] = ['NAXIS2', string(height, format='(i0)'), 'length of data axis 2']
header[*, 5] = ['EXTEND', 'T', 'FITS dataset may contain extensions']
header[*, 6] = ['COMMENT', '', "  FITS (Flexible Image Transport System) format is defined in 'Astronomy"]
header[*, 7] = ['COMMENT', '', "  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"]
header[*, 8] = ['OBJECT', tname, 'Target Name']
header[*, 9] = ['POL', pol, 'Polarization']
header[*,10] = ['BUNIT', 'sigmas', 'Standard deviations of the noise power']
header[*,11] = ['WCSAXES', '2', '']
header[*,12] = ['WCSNAME', 'delay-Doppler', '']
header[*,13] = ['CRVAL1', '0.', '[Hz]']
header[*,14] = ['CUNIT1', 'Hz', '']
header[*,15] = ['CTYPE1', 'DOPPLER', 'Doppler Shift']
header[*,16] = ['CRVAL2', '0.', '[us]']
header[*,17] = ['CUNIT2', 'us', '']
header[*,18] = ['CTYPE2', 'DELAY', 'Round-Trip Time Delay']
caldat_roundsec,systime(/julian,/utc),month,day,year,hour,minute,second
date = string(year, format='(i4)') + '-' + string(month, format='(i2.2)') + '-'   $
       + string(day, format='(i2.2)') + 'T' + string(hour, format='(i2.2)') + ':' $
       + string(minute, format='(i2.2)') + ':' + string(second, format='(i2.2)')
header[*,19] = ['DATE', date, 'file creation date (YYYY-MM-DDThh:mm:ss UT)']
header[*,20] = ['HISTORY', '', 'Radar image rewritten by IDL procedure writefitsi']

; Add the 'DATE-OBS' header record (and the associated comment records)
; if the necessary information is present in this image's extra tags

date_obs = getextrai('date-obs')
if isnull(date_obs) then begin
  year = getextrai('year')
  month = getextrai('month')
  day = getextrai('day')
  hour = getextrai('hour')
  minute = getextrai('minute')
  second = getextrai('second')
  if notnull(year) and notnull(month) and notnull(day) $
                   and notnull(hour) and notnull(minute) and notnull(second) then begin
    date_obs = string(year, format='(i4)') + '-' + string(month, format='(i2.2)') + '-'   $
               + string(day, format='(i2.2)') + 'T' + string(hour, format='(i2.2)') + ':' $
               + string(minute, format='(i2.2)') + ':' + string(second, format='(i2.2)')
    year_comment = getextrai('year',/comment)
    sdev = getextrai('sdev')
    date_obs_refers_to_rx_start = ((strpos(year_comment,'mean RX') eq -1) and isnull(sdev))
  endif
endif else begin
  date_obs_comment = getextrai('date-obs',/comment)
  date_obs_refers_to_rx_start = (strpos(date_obs_comment,'RX start') ne -1)
endelse
if notnull(date_obs) then begin
  if date_obs_refers_to_rx_start then begin
    newrecords = strarr(3,4)
    newrecords[*,0] = ['DATE-OBS', date_obs, 'UT date of RX start']
    newrecords[*,1] = ['COMMENT', '', 'The RX start time is a reasonable approximation to the experiment']
    newrecords[*,2] = ['COMMENT', '', 'mid-time. A more accurate mid-time is DATE-OBS + (EXPTIME-RTT)/2']
    newrecords[*,3] = ['COMMENT', '', 'which is typically about 3s earlier']
    header = transpose( [ transpose(header), transpose(newrecords) ] )
    nhead = nhead + 4
  endif else begin
    newrecord = ['DATE-OBS', date_obs, 'UT mean RX date']
    header = transpose( [ transpose(header), transpose(newrecord) ] )
    nhead = nhead + 1
  endelse
endif

; Add or overwrite the 'CRPIX1' and 'CRPIX2' header records if the necessary information
; is present in this image's extra tags
;
; CRPIX1 is just eph_col + 1, and CRPIX2 = eph_row + 1; but we need to do this math carefully
; (by treating only the integer portion of the input string as a number) so that we don't
; inadvertently change the number of decimal places.

k = where(extratags.name eq 'eph_col')
if k ne -1 then begin
  eph_col_str = extratags[k].value
  decimalpoint = strpos(eph_col_str,'.')
  int_crpix1 = long(strmid(eph_col_str,0,decimalpoint)) + 1
  crpix1 = string(int_crpix1, format='(i0)') + strmid(eph_col_str,decimalpoint)
  newrecord = ['CRPIX1', crpix1, 'Center of first px is 1.0']
  header = transpose( [ transpose(header), transpose(newrecord) ] )
  nhead = nhead + 1
  skip_crpix1 = 1
endif else begin
  skip_crpix1 = 0
endelse
k = where(extratags.name eq 'eph_row')
if k ne -1 then begin
  eph_row_str = extratags[k].value
  decimalpoint = strpos(eph_row_str,'.')
  int_crpix2 = long(strmid(eph_row_str,0,decimalpoint)) + 1
  crpix2 = string(int_crpix2, format='(i0)') + strmid(eph_row_str,decimalpoint)
  newrecord = ['CRPIX2', crpix2, 'Center of first px is 1.0']
  header = transpose( [ transpose(header), transpose(newrecord) ] )
  nhead = nhead + 1
  skip_crpix2 = 1
endif else begin
  skip_crpix2 = 0
endelse

; Go through all the extra tags and add to the string-format header
; any tags that we haven't already added.

for n=0L,nextra-1 do begin

  ; Check for rdf tags whose names and/or comments differ from the names and comments
  ; of the corresponding FITS tags (as written by the "writerdf" and "writefits" Perl scripts)
  ;
  ; Change the tag names to the FITS tag names; if the comments still have
  ; the standard rdf values, change them to the standard FITS comments
  ;
  ; Leave the 'DOY' tag comment unchanged: the standard rdf comment used to be silent with
  ; respect to mean RX vs. RX start, so we can't always know which one belongs in the comment.
  ;
  ; At this point we also convert all tag names to upper case, and we convert
  ; records with null-string names to 'HISTORY' and 'COMMENT' records

  k = where(rdftagnames eq extratags[n].name)
  if k eq -1 then begin
    if notnull(extratags[n].name) then begin
      tagname = strupcase(extratags[n].name)
      comment = strmid(extratags[n].comment, 2)
    endif else if strmid(extratags[n].comment,0,10) eq '# HISTORY ' then begin
      tagname = 'HISTORY'
      comment = strmid(extratags[n].comment, 10)
    endif else begin
      tagname = 'COMMENT'
      comment = strmid(extratags[n].comment, 2)
    endelse
    value = extratags[n].value
  endif else begin
    tagname = fitstagnames[k]
    value = extratags[n].value
    if tagname eq 'DOY' then begin
      comment = extratags[n].comment
    endif else begin
      comment = (extratags[n].comment eq rdfcomments[k]) ? fitscomments[k] : extratags[n].comment
    endelse
  endelse

  ; We won't use the existing 'DATE' tag because we'll create a new one later;
  ; we won't output the 'INFILE' tag because it can be longer than 80 bytes
  ; (given that it includes the full path);
  ; we won't use the existing 'CRPIX1' and 'CRPIX2' tags if we already got them
  ; from the 'eph_col' and 'eph_row' rdf tags

  usetag = (tagname ne 'DATE' and tagname ne 'INFILE')    $
           and not (tagname eq 'CRPIX1' and skip_crpix1)  $
           and not (tagname eq 'CRPIX2' and skip_crpix2)

  ; If we haven't already added this tag to the string-format header, do so now

  k = where(strtrim(header[0,*],2) eq tagname, count)
  if count gt 0 then begin
    for m=0L,count-1 do begin
      if (strtrim(header[1,k[m]],2) eq strtrim(value,2)) and $
         (strtrim(header[2,k[m]],2) eq strtrim(comment,2)) then usetag = 0
    endfor
  endif
  if usetag then begin
    header = transpose( [ transpose(header), transpose( [tagname, value, comment] ) ] )
    nhead = nhead + 1
  endif
endfor

; Add some comment records that are too long to fit in the same record as the
; associated key and value

k1 = where(header[0,*] eq 'FREQSAMP')
k2 = where(header[0,*] eq 'COMMENT' and $
           header[2,*] eq 'If FREQSAMP < FFTLEN, then the transform was zero-filled')
if k1 ne -1 and k2 eq -1 then begin
  newrecord = ['COMMENT', '', 'If FREQSAMP < FFTLEN, then the transform was zero-filled']
  header = transpose( [ transpose(header), transpose(newrecord) ] )
  nhead = nhead + 1
endif
k1 = where(header[0,*] eq 'CODEPROC')
k2 = where(header[0,*] eq 'COMMENT' and $
           header[2,*] eq 'CODEPROC is the code processing method: values are short, long_orig, or')
if k1 ne -1 and k2 eq -1 then begin
  newrecords = strarr(3,6)
  newrecords[*,0] = ['COMMENT', '', 'CODEPROC is the code processing method: values are short, long_orig, or']
  newrecords[*,1] = ['COMMENT', '', 'long_mod. short means short-code data processed so that each Doppler fft']
  newrecords[*,2] = ['COMMENT', '', ' includes data from all samples per baud; long_mod means long-code data']
  newrecords[*,3] = ['COMMENT', '', 'processed in the same way; long_orig means long-code data where each sam']
  newrecords[*,4] = ['COMMENT', '', 'ple per baud is processed separately and the separate images are then in']
  newrecords[*,5] = ['COMMENT', '', 'terleaved in delay.']
  header = transpose( [ transpose(header), transpose(newrecords) ] )
  nhead = nhead + 6
endif else if k2 ne -1 then begin
  header[2,k2+2] = ' ' + header[2,k2+2]
endif
k1 = where(header[0,*] eq 'SPB')
k2 = where(header[0,*] eq 'COMMENT' and $
           header[2,*] eq 'SPB is the number of input samples taken per baud.')
k3 = where(header[0,*] eq 'RPB')
k4 = where(header[0,*] eq 'COMMENT' and $
           header[2,*] eq 'RPB is the number of output image rows saved per baud.')
k5 = where(header[0,*] eq 'COMMENT' and $
           header[2,*] eq 'RPB and SPB may differ either because data were thrown away, or because')
if k1 ne -1 and k2 eq -1 then begin
  newrecord = ['COMMENT', '', 'SPB is the number of input samples taken per baud.']
  header = transpose( [ transpose(header), transpose(newrecord) ] )
  nhead = nhead + 1
endif
if k3 ne -1 and k4 eq -1 then begin
  newrecord = ['COMMENT', '', 'RPB is the number of output image rows saved per baud.']
  header = transpose( [ transpose(header), transpose(newrecord) ] )
  nhead = nhead + 1
endif
if k1 ne -1 and k3 ne -1 and k5 eq -1 then begin
  newrecords = strarr(3,2)
  newrecords[*,0] = ['COMMENT', '', 'RPB and SPB may differ either because data were thrown away, or because']
  newrecords[*,1] = ['COMMENT', '', 'the extra samples were used in the frequency processing.']
  header = transpose( [ transpose(header), transpose(newrecords) ] )
  nhead = nhead + 2
endif

; Add the 'END' record marking the end of the FITS header

newrecord = ['END', '', '']
header = transpose( [ transpose(header), transpose(newrecord) ] )
nhead = nhead + 1

; Initialize a byte array which will contain properly formatted 80-byte header records, ordered as the
; "writefits" Perl script would order them; start with just eight records, since we know for sure that
; the first eight string-format records will fit in one 80-byte record each
;
; Also initialize a vector which will tell us which of the string-format records we have already
; converted to byte format

nrec = 8L
header_byte = bytarr(80,nrec)
for n=0L,nrec-1 do header_byte[*,n] = fitsheadrecord(header[0,n], header[1,n], header[2,n])
usedtag = intarr(nhead)
usedtag[0:nrec-1] = 1

; Add all other standard header records (those output by the "writefits" Perl script) to the byte array

processtag, 'DATE-OBS', header, header_byte, usedtag, nrec, /quotes
skip = (date_obs_refers_to_rx_start) ? 0 : 1
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='The RX start time is a reasonable approximation to the experiment',skip=skip
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='mid-time. A more accurate mid-time is DATE-OBS + (EXPTIME-RTT)/2',skip=skip
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='which is typically about 3s earlier',skip=skip
processtag, 'DOY', header, header_byte, usedtag, nrec
processtag, 'SFM', header, header_byte, usedtag, nrec
processtag, 'JDSTART', header, header_byte, usedtag, nrec
processtag, 'JDMEAN', header, header_byte, usedtag, nrec
processtag, 'JDEND', header, header_byte, usedtag, nrec
processtag, 'OBJECT', header, header_byte, usedtag, nrec, /quotes
processtag, 'BUNIT', header, header_byte, usedtag, nrec, /quotes
processtag, 'WCSAXES', header, header_byte, usedtag, nrec
processtag, 'WCSNAME', header, header_byte, usedtag, nrec, /quotes
processtag, 'CRVAL1', header, header_byte, usedtag, nrec
processtag, 'CRPIX1', header, header_byte, usedtag, nrec
processtag, 'CDELT1', header, header_byte, usedtag, nrec
processtag, 'CUNIT1', header, header_byte, usedtag, nrec, /quotes
processtag, 'CTYPE1', header, header_byte, usedtag, nrec, /quotes
processtag, 'CRVAL2', header, header_byte, usedtag, nrec
processtag, 'CRPIX2', header, header_byte, usedtag, nrec
processtag, 'CDELT2', header, header_byte, usedtag, nrec
processtag, 'CUNIT2', header, header_byte, usedtag, nrec, /quotes
processtag, 'CTYPE2', header, header_byte, usedtag, nrec, /quotes
processtag, 'PC1_1', header, header_byte, usedtag, nrec
processtag, 'PC1_2', header, header_byte, usedtag, nrec
processtag, 'PC2_2', header, header_byte, usedtag, nrec
processtag, 'EXPOSURE', header, header_byte, usedtag, nrec
processtag, 'LAMBDA', header, header_byte, usedtag, nrec
processtag, 'TXOFFSET', header, header_byte, usedtag, nrec
processtag, 'CODELEN', header, header_byte, usedtag, nrec
processtag, 'BAUDLEN', header, header_byte, usedtag, nrec
processtag, 'LOOPTIME', header, header_byte, usedtag, nrec
processtag, 'COHAVG', header, header_byte, usedtag, nrec
processtag, 'CODEBW', header, header_byte, usedtag, nrec
processtag, 'BANDWIDT', header, header_byte, usedtag, nrec
processtag, 'EPH_ROW', header, header_byte, usedtag, nrec
processtag, 'EPH_COL', header, header_byte, usedtag, nrec
processtag, 'DC_COL', header, header_byte, usedtag, nrec
processtag, 'FFTLEN', header, header_byte, usedtag, nrec
processtag, 'LOOKS', header, header_byte, usedtag, nrec
processtag, 'FREQSAMP', header, header_byte, usedtag, nrec
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='If FREQSAMP < FFTLEN, then the transform was zero-filled'
processtag, 'RANGERES', header, header_byte, usedtag, nrec
processtag, 'SPB', header, header_byte, usedtag, nrec
processtag, 'RPB', header, header_byte, usedtag, nrec
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='SPB is the number of input samples taken per baud.'
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='RPB is the number of output image rows saved per baud.'
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='RPB and SPB may differ either because data were thrown away, or because'
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='the extra samples were used in the frequency processing.'
processtag, 'CODEPROC', header, header_byte, usedtag, nrec, /quotes
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='CODEPROC is the code processing method: values are short, long_orig, or'
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='long_mod. short means short-code data processed so that each Doppler fft'
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment=' includes data from all samples per baud; long_mod means long-code data'
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='processed in the same way; long_orig means long-code data where each sam'
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='ple per baud is processed separately and the separate images are then in'
processtag, 'COMMENT', header, header_byte, usedtag, nrec, $
           comment='terleaved in delay.'
processtag, 'EPHEMERI', header, header_byte, usedtag, nrec, /quotes
processtag, 'POL', header, header_byte, usedtag, nrec, /quotes
processtag, 'DATE', header, header_byte, usedtag, nrec, /quotes
processtag, 'HISTORY', header, header_byte, usedtag, nrec, $
           comment='Radar image rewritten by IDL procedure writefitsi'

; Add any remaining string-format header records to the byte array, finishing with the END record

for k=0L,nhead-1 do begin
  if not usedtag[k] then begin
    if header[0,k] eq 'TIMEZONE' or header[0,k] eq 'CHECKSUM' or header[0,k] eq 'DATASUM' then begin
      newrecord = fitsheadrecord(header[0,k], header[1,k], header[2,k], /quotes)
    endif else begin
      newrecord = fitsheadrecord(header[0,k], header[1,k], header[2,k])
    endelse
    if notnull(newrecord) then begin
      n_new = (size(newrecord))[2]
      header_byte = transpose( [ transpose(header_byte), transpose(newrecord) ] )
      nrec = nrec + n_new
    endif
  endif
endfor

; Open the output file

err = openOutfile(lun,outfile,'fits',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the header

writeu,lun,header_byte

; Pad the header with blanks to make it an integer multiple of 2880 bytes (36 records)

n_padbytes = (36*ceil(nrec/36.0) - nrec)*80L
if n_padbytes gt 0 then begin
  blank = 32B
  writeu,lun,replicate(blank,n_padbytes)
endif

; Write the image

writeu,lun,swap_endian(image, /swap_if_little_endian)

; Pad the image with ASCII nulls to make it an integer multiple of 2880 bytes

n_imagebytes = width*height*4
n_padbytes = 2880*ceil(n_imagebytes/2880.0) - n_imagebytes
if n_padbytes gt 0 then begin
  null = 0B
  writeu,lun,replicate(null,n_padbytes)
endif

; Close up shop

close,lun
free_lun,lun
print,'Wrote 1 ',pol,' image to FITS file ',outfile,format='(4a)'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
