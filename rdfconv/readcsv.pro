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

; Make sure it's a string to avoid confusion
headline = ''

; Get new group numbers for placing these spectra into the
; stack and/or images into the image stack. We increment if
; the number of channels changes
groupCW = (nstack gt 0) ? max(stackgroups()) + 0L : 0L
oldndata = -1L

extratag = {format:'', name:'', value:'', comment:''}


if n_params() ne 1L or keyword_set(help) then begin
  print,' '
  print,"readcsv,infile[,/silent][,/help]"
  print,' '
  print,"        infile can be a wildcarded string. Each file will be added to the stack"
  print,' '
  print,'        /silent  turns off screen output for individual input frames'
  print,' '
  return
endif else if size(infile, /type) ne 7 then begin
  print,' '
  print,"readrdf,infile[,/silent][,/help]"
  print,'Make sure that infile is a quoted string!'
  print,' '
  return
endif

filelist = file_search(infile, /FULLY_QUALIFY_PATH, count=filecount)

if filecount lt 1L then begin
  print," "
  print,"ERROR: No files matched ",infile
  print," "
  return
endif

for ifile = 0L, filecount-1 do begin
  infile = filelist[ifile]
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

; Get the tag names

  tagnames = strlowcase(tag_names(blanktags1()))
  tags = blanktags()
  ntags = n_elements(tagnames)

  ; Initialize some parameters

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
    print,infile
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
    print, infile
    print,' '
    close,lun
    free_lun,lun
    return
  endif

  if stregex(headline, "^ *# *keywords", /fold) lt 0 then begin
    print, ' '
    print, "ERROR in readcsv: This is not an archive-format CSV file:"
    print, infile
    print, " "
    free_lun, lun
    return
  endif

; Use read_csv to take care of dealing with quoted commas, etc, then parse lines.
; If more than 4 columns, they we're dealing with future stokes data, so it's OK
; if they are numbers after that
; don't need the lun anymore.

  free_lun, lun

  if not keyword_set(silent) then print, "Reading csv file ", infile
  csv=read_csv(infile, count=count, types=replicate('string',4))

; First check for markers
    
  markers=['^ *tags','^ *extra *tags','^ *# *column *definition','^ *# *data']
  wmark = lonarr(n_elements(markers))
  for i = 0, n_elements(markers)-1 do begin
    wmark[i] = where(stregex(csv.field1,markers[i], /bool, /fold) gt 0)
    if wmark[i] le 0 then begin
      print," "
      print,"Could not find marker ", markers[i]
      print, " "
      return
    endif
  endfor

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

  ind = where(stregex(csv.field1, 'target *name', /fold, /bool))
  if ind[0] lt 0 then tname = 'Unknown' else tname = csv.field2[ind[0]]
  
  ; now tags.
  for itag = 0, n_elements(tagnames) -1 do begin
    row = where(strmatch(csv.field1[wmark[0]+1:wmark[1]-1], tagnames[itag], /fold))
    if row[0] lt 0 then continue
    if n_elements(row) ne 1 then begin
      print, ' '
      print, 'ERROR: found duplicate tags in csv file'
      print, ' '
      return
    endif
    tags[0].(itag) = double(csv.field2[row[0]+wmark[0]+1])
    tags[1].(itag) = double(csv.field3[row[0]+wmark[0]+1])
  endfor

  ; extra tags. Go ahead and put the description in the rdf file, why not
  ; Check for 'target': We'll add it in if it's not already there.
  for itag = wmark[1] + 1, wmark[2] - 1 do begin
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
  if ndata ne oldndata then begin
    groupCW = groupCW + 1L
    oldndata = ndata
  endif
  spec = transpose([[oc],[sc]])

  stackStruc = {group:groupCW, spec:spec, tags:tags, $
    extratags:extratags, ndata:ndata, ntags:ntags, nextra:nextra, $
    tname:tname}


 if nstack eq 0 then begin
    stack = ptrarr(1, /allocate_heap)
    *stack[0] = stackStruc
  endif else begin
    stack = [stack, ptr_new(stackStruc)]
  endelse
  nstack = nstack + 1L
  if not keyword_set(silent) then $
    print,'Added pair #',nstack,' to the pair stack',format='(a,i0,a)'

endfor ; loop over files

return
end
