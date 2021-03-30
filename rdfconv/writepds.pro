forward_function getextrai

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function openInfile,lun,infile,insuffix,get_lun=get_lun,silent=silent

; Open an input file
;    infile:   on input,  infile is the file name (or stem);
;              on output, infile is the exact filename opened,
;                         or the null string if there was a problem
;    insuffix: the default suffix to be added to the stem (infile)
;    lun:      the logical unit opened, or -999 if there was a problem
;
; Return the error number from the "openr" command
;
; 2006 Jun 28: Allow insuffix to be a vector of suffices to try

; Assign "failure" values to output arguments;
; these will later be replaced if the operation succeeds

infile_input = infile
infile = ''

; Validate input arguments

if n_params() eq 2 then begin
  insuffix = ''
endif else if n_params() ne 3 then begin
  print,' '
  print,'err = openInfile(lun,infile[,suffix][,/get_lun])'
  print,'ERROR in openInfile: invalid number of arguments'
  print,' '
  return, 1
endif

if not keyword_set(get_lun) then get_lun = 0

if n_elements(lun) eq 0 and not get_lun then begin
  print,' '
  print,'err = openInfile(lun,infile[,suffix][,/get_lun][,/silent])'
  print,'ERROR in openInfile: lun is undefined and /get_lun is not set'
  print,' '
  return, 1
endif

; Create a sorted list of unique filenames to try -- starting with infile
; (i.e., filename given in full)

infile_try = [infile_input]
for k=0L,n_elements(insuffix)-1 do begin

  ; Define a variable which is true if the default suffix is the null string
  ; or if the file name already ends in the default suffix

  sufstring = '\.' + insuffix[k] + '$'
  uselessSuffix = isnull(insuffix[k]) or stregex(infile_input, sufstring) ne -1
  if not uselessSuffix then begin
    if strmid(insuffix[k],0,1) eq '.' xor strmid(infile_input,0,/reverse) eq '.' then begin
      fullname = infile_input + insuffix[k]
    endif else if strmid(insuffix[k],0,1) ne '.' then begin
      fullname = infile_input + '.' + insuffix[k]
    endif else begin
      fullname = strmid(infile_input,0,strlen(infile_input)-1) + insuffix[k]
    endelse
    dummy = where(infile_try eq fullname, count)
    if count eq 0 then infile_try = [infile_try, fullname]
  endif
endfor
n_try = n_elements(infile_try)

; Go through the list of filenames to find one that works

nfiles_save = lonarr(n_try)
found_unique_file = 0
k = -1L
repeat begin
  k = k + 1

  ; Try to find a unique file name to open

  filenames = findfile(infile_try[k], count=nfiles)
  if nfiles eq 1 then begin
    found_unique_file = 1
    unique_file_to_try = filenames[0]
  endif else begin
    nfiles_save[k] = nfiles
  endelse
endrep until (found_unique_file or k eq n_try-1)

; Let the user know if the file wasn't found

if not found_unique_file then begin
  print,' '
  formatstring = (nfiles_save[0] lt 100) ? '(a,i2,3a)' : '(a,i0,3a)'
  print,'ERROR in openInfile: ',nfiles_save[0],' files ',infile_try[0],' exist', $
        format=formatstring
  for k=1L,n_try-1 do begin
    formatstring = (nfiles_save[k] lt 100) ? '(a,i2,3a)' : '(a,i0,3a)'
    print,'                     ',nfiles_save[k],' files ',infile_try[k],' exist', $
          format=formatstring
  endfor
  print,' '
  return, 1
endif

; Try to open the file

openr,lun,unique_file_to_try,error=err,get_lun=get_lun
if err eq 0 then begin
  infile = unique_file_to_try
  if not keyword_set(silent) then print,'Successfully opened file ',infile
endif else begin
  print,' '
  print,'ERROR in openInfile: unsuccessful attempt to open file ',unique_file_to_try
  print,' '
endelse

return, err
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function openOutfile,lun,outfile,outsuffix,get_lun=get_lun,overwrite=overwrite,append=append

; Open an output file
;    outfile:   on input,  outfile is the file name (or stem);
;               on output, outfile is the exact filename opened,
;                          or the null string if there was a problem
;    outsuffix: the default suffix to be added to the stem (outfile)
;    lun:       the logical unit number opened, or -999 if there was a problem
;
; Return the error number from the "openw" command

; Assign "failure" values to output arguments;
; these will later be replaced if the operation succeeds

outfile_input = outfile
outfile = ''

; Validate input arguments

if n_params() eq 2 then begin
  outsuffix = ''
endif else if n_params() ne 3 then begin
  print,' '
  print,'err = openOutfile(lun,outfile[,suffix][,/get_lun][,/overwrite or /append])'
  print,'ERROR in openOutfile: invalid number of arguments'
  print,' '
  return, 1
endif

if not keyword_set(get_lun) then get_lun = 0

if keyword_set(overwrite) and keyword_set(append) then begin
  print,"ERROR in openOutfile: Can't set both overwrite and append keywords"
  return, 1
endif else if n_elements(lun) eq 0 and not get_lun then begin
  print,' '
  print,'err = openOutfile(lun,outfile[,suffix][,/get_lun][,/overwrite or /append])'
  print,'ERROR in openInfile: lun is undefined and /get_lun is not set'
  print,' '
  return, 1
endif

; Define the file name to open

if isnull(outsuffix) then begin
  outfile_try = outfile_input
endif else begin
  sufstring = '\.' + outsuffix + '$'
  if stregex(outfile_input, sufstring) eq -1 then begin
    outfile_try = outfile_input + '.' + outsuffix
  endif else begin
    outfile_try = outfile_input
  endelse
endelse

; Try to open the file

if keyword_set(overwrite) then begin
  openw,lun,outfile_try,error=err,get_lun=get_lun
endif else if keyword_set(append) then begin
  openw,lun,outfile_try,error=err,get_lun=get_lun,/append
endif else begin
  filenames = findfile(outfile_try, count=nfiles)
  if nfiles eq 0 then begin
    openw,lun,outfile_try,error=err,get_lun=get_lun
  endif else begin
    userchoice = '' ; define as a string for the read procedure
    userprompt = 'File ' + outfile_try + ' already exists: overwrite (o), append (a), or cancel (c)? '
    read,userchoice,prompt=userprompt
    userchoice = strlowcase(strmid(strtrim(userchoice), 0, 1))
    if userchoice eq 'o' then begin
      openw,lun,outfile_try,error=err,get_lun=get_lun
    endif else if userchoice eq 'a' then begin
      openw,lun,outfile_try,error=err,get_lun=get_lun,/append
    endif else begin
      print,'Write operation canceled'
      return, 1
    endelse
  endelse
endelse

if err eq 0 then begin
  outfile = outfile_try
endif else begin
  print,' '
  print,'ERROR in openOutfile: unsuccessful attempt to open file ',outfile_try
  print,' '
endelse

return, err
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function addsuffix,instring,oldsuffix,newsuffix,nodot=nodot

; Take an input string (filename) and append a new suffix.  If the input
; string already ends in the new suffix, leave it unchanged; if instead it
; ends in a specified old suffix, remove that old suffix before appending
; the new suffix.
;
; Setting /nodot forces addsuffix NOT to prepend a period to all suffices
; 
; Written 2005 June 29 by CM

if keyword_set(nodot) then begin
  oldsuffix_use = oldsuffix
  newsuffix_use = newsuffix
endif else begin
  oldsuffix_use = '.' + oldsuffix
  newsuffix_use = '.' + newsuffix
endelse

inLen = strlen(instring)
oldsufLen = strlen(oldsuffix_use)
newsufLen = strlen(newsuffix_use)
outLen = inLen + newsufLen    ; for starters

; Check whether or not the input string ends in the old suffix

found_oldsuffix = 0
if (oldsufLen gt 0 and oldsufLen le inLen) then begin
  found_oldsuffix = (strmid(instring, (inLen - oldsufLen)) eq oldsuffix_use)
  if found_oldsuffix then outLen = inLen - oldsufLen + newsufLen
endif

; If the input string doesn't end in the old suffix,
; check whether or not it already ends in the new suffix

if (not found_oldsuffix and newsufLen gt 0 and newsufLen le inLen) then begin
  found_newsuffix = (strmid(instring, (inLen - newsufLen - 1)) eq newsuffix_use)
  if found_newsuffix then outLen = inLen
endif

; Create the output string: copy n_copy characters from the
; beginning of instring to the beginning of outstring, then
; concatenate the new suffix to the end of outstring.

n_copy = outLen - newsufLen
outstring = strmid(instring, 0, n_copy) + newsuffix_use

return, outstring
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function fitsheadrecord,key,value,comment,quotes=quotes

; Create an 80-byte FITS header record and return it as a byte array;
; return the null string if there's any problem
;
; /quotes places the value in single quotes

; Deal with the key

keyuse = strupcase(key)
keylen = strlen(keyuse)
if keylen eq 0 then begin
  keyuse = 'COMMENT '
  keylen = 8
endif else if keylen lt 8 then begin
  keyuse = keyuse + strjoin(replicate(' ',8-keylen))
endif else if keylen gt 8 then begin
  print,"WARNING in fitsheadrecord: key name '",keyuse,"' will be truncated"
  keyuse = strmid(keyuse,0,8)
endif
headrecord = keyuse

; Deal with the value; check that it won't extend this record beyond 80 bytes

if keyuse ne 'COMMENT ' and keyuse ne 'HISTORY ' and keyuse ne 'END     ' then begin
  valtype = size(value, /type)
  if (valtype ge 1 and valtype le 3) or (valtype ge 12 and valtype le 15) then begin
    valstring = string(value, format='(i0)')
  endif else if valtype eq 4 then begin
    valstring = string(value, format='(f0)')
  endif else if valtype eq 5 then begin
    valstring = string(value, format='(d0)')
  endif else if valtype eq 7 then begin
    valstring = value
  endif else begin
    valstring = strtrim(string(value), 2)
  endelse
  vallen = strlen(valstring)
  if keyword_set(quotes) then begin
    if vallen lt 8 then begin
      valstring = "'" + valstring + strjoin(replicate(' ',8-vallen)) + "'"
      vallen = 10
    endif else begin
      valstring = "'" + valstring + "'"
      vallen = vallen + 2
    endelse
  endif
  if vallen gt 72 then begin
    print,"WARNING in fitsheadrecord: key '",strtrim(keyuse,2), $
          "' will be omitted because the value is too long",format='(3a)'
    print,'        (value = ',valstring,')'
    return,''
  endif
  if vallen lt 20 then begin
    padding = strjoin(replicate(' ',20-vallen))
    valstring = (keyword_set(quotes)) ? valstring + padding : padding + valstring
  endif
  headrecord = headrecord + '= ' + valstring
  if notnull(comment) then headrecord = headrecord + ' / '
endif

; Deal with the comment, continuing onto an extra record(s)
; if the comment would otherwise extend this record beyond 80 bytes

nrecords = 1L
if isnull(comment) then begin
  nrecords = 1L
  headrecord = [headrecord]
endif else begin
  if strmid(comment,0,2) eq '# ' then begin
    commentuse = strmid(comment,2)
  endif else begin
    commentuse = comment
  endelse
  nskipbytes = (strlen(headrecord) < 50)
  maxcommentlen = 80 - nskipbytes
  startstring = 'COMMENT' + strjoin(replicate(' ',nskipbytes-7))
  commentlen = strlen(commentuse)
  if commentlen le maxcommentlen then begin
    headrecord = [headrecord + commentuse]
    commentuse = ''
  endif else begin
    commentpiece = strmid(commentuse,0,maxcommentlen)
    k = strpos(commentpiece,' ',/reverse_search)
    if k eq -1 then k = maxcommentlen - 1
    headrecord = [headrecord + strmid(commentuse,0,k)]
    commentuse = strmid(commentuse,k+1)
  endelse
  commentlen = strlen(commentuse)
  while commentlen gt 0 do begin
    if commentlen le maxcommentlen then begin
      continuation = startstring + commentuse
      commentuse = ''
    endif else begin
      commentpiece = strmid(commentuse,0,maxcommentlen)
      k = strpos(commentpiece,' ',/reverse_search)
      if k eq -1 then k = maxcommentlen - 1
      continuation = startstring + strmid(commentuse,0,k)
      commentuse = strmid(commentuse,k+1)
    endelse
    nrecords = nrecords + 1
    headrecord = [headrecord, continuation]
    commentlen = strlen(commentuse)
  endwhile
endelse

; Convert from string to byte

blank = 32B
headrecord_byte = replicate(blank,80,nrecords)
for n=0L,nrecords-1 do begin
  headrecordlen = strlen(headrecord[n])
  headrecord_byte[0:headrecordlen-1,n] = byte(headrecord[n])
endfor
return,headrecord_byte
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro plotsum,smooth=smooth,xrange=xrange,maxf=maxf,full=full,mu=mu,maxmu=maxmu,  $
            chan=chan,separate=separate,rcs=rcs,help=help,_extra=_ext

; Plot an OC/SC pair
;
; Modified 2005 Jan 05 by CM: For some reason the "linestyle" keyword to oplot
;     is ignored unless the window is at least partly hidden while being drawn,
;     so I hide it, draw it, and then show it.
;
; Modified 2006 Jun 25 by CM: Frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded

; Check that frequency range was correctly input (if input at all)

if n_elements(xrange) gt 0 then begin
  if n_elements(maxf) gt 0 or keyword_set(full) or n_elements(xrange) ne 2 then begin
    xrangeOK = 0
  endif else begin
    xrangeOK = (xrange[0] lt xrange[1])
  endelse
endif else if n_elements(maxf) gt 0 then begin
  if keyword_set(full) or n_elements(maxf) ne 1 then begin
    xrangeOK = 0
  endif else begin
    xrangeOK = (maxf gt 0)
  endelse
endif else begin
  xrangeOK = 1
endelse

if n_params() ne 0 or not xrangeOK or keyword_set(help) then begin
  print,' '
  print,'plotsum[,smooth=smooth][,xrange=[minf,maxf]][,maxf=maxf][,/full][,mu=mu]  $'
  print,'       [,maxmu=maxmu][,chan=1 or 2][,/separate][,/rcs][,plot keywords]    $'
  print,'       [,/help]'
  print,' '
  print,'       Default frequency range is [-2000,2000] Hz or else the full spectrum,'
  print,'           whichever is narrower'
  print,' '
  print,'       The next three keywords are mutually exclusive:'
  print,' '
  print,'           For xrange, must have maxf > minf even if', $
                        ' frequency increases leftward',format='(2a)'
  print,'           maxf=100 is the same as xrange=[-100,100]'
  print,'           /full plots the full spectrum'
  print,' '
  print,'       mu = 1 plots polarization ratio SC/OC'
  print,'            (values <= 0 in *either* spectrum are assigned SC/OC = 0)'
  print,'       mu > 1 plots SC/OC after low-pass filtering each spectrum'
  print,'            (larger mu --> stronger filtering)'
  print,'       maxmu is the maximum polarization ratio to plot'
  print,'            (default ~ the smaller of 1.2 or the maximum ratio in the data)'
  print,' '
  print,'       /separate produces separate OC and SC plots; otherwise the'
  print,'            OC and SC spectra are superimposed within the same plot'
  print,' '
  print,'       /rcs expresses signal strength as radar cross section density in'
  print,'            km^2 per Hz; otherwise the unit is the rms noise'
  print,' '
  return
endif else if (*loaded).ndata le 2 then begin
  print,"ERROR in plotsum: There's no loaded pair to plot"
  return
endif

; Check which channel(s) to write

if not keyword_set(chan) then begin
  plotpol = [1,1]
endif else if chan eq 1 or chan eq 2 then begin
  plotpol = (chan eq 1) ? [1,0] : [0,1]
  if keyword_set(separate) then begin
    print,"ERROR in plotsum: Can't set /separate if 'chan' keyword is used"
    return
  endif
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse
nplots = (keyword_set(separate)) ? total(plotpol) : 1L

; Get some elements of the loaded pair

f = (*loaded).freq
oc = reform((*loaded).spec[0,*])
sc = reform((*loaded).spec[1,*])
tags = (*loaded).tags
oc_sdev = tags[0].sdev
sc_sdev = tags[1].sdev
df = tags[0].dfreq
posfr = tags[0].posfr
bin0 = tags[0].xjcen
ndata = (*loaded).ndata

; Express in absolute units if desired

if keyword_set(rcs) then begin
  oc = oc * oc_sdev / df
  sc = sc * sc_sdev / df
endif

; Smooth the spectra if desired

oc2 = oc
sc2 = sc
if n_elements(smooth) gt 0 then begin
  if smooth ge 2 then begin
    oc2 = smooth(oc2,smooth)
    sc2 = smooth(sc2,smooth)
    if not keyword_set(rcs) then begin
      oc2 = oc2*sqrt(smooth)
      sc2 = sc2*sqrt(smooth)
    endif
  endif
endif

; Prepare the polarization ratio ("mu") plot if desired

pmultisave = !p.multi

domu = 0
if n_elements(mu) gt 0 then begin
  domu = 1
  nplots = nplots + 1
  if mu gt 1 then begin
    sc3 = convol(sc2,digital_filter(0,1.0/mu,50,mu*2))
    oc3 = convol(oc2,digital_filter(0,1.0/mu,50,mu*2))
  endif else begin
    sc3 = sc2
    oc3 = oc2
  endelse
  positivemask = sc3 gt 0 and oc3 gt 0
  numer = sc3
  denom = oc3 > 1.0e-10
  if not keyword_set(rcs) then begin
    numer = numer*sc_sdev
    denom = denom*oc_sdev
  endif
  ratio = positivemask*numer/denom
endif

!p.multi = [0,1,nplots]
erase
if !d.window ge 0 then wshow,!d.window,0  ; temporary fix: hide the window while drawing

; Get the y-limits for the spectral plot(s)

if n_elements(xrange) gt 0 then begin
  xr = (posfr eq 1) ? xrange : reverse(xrange)
endif else if n_elements(maxf) gt 0 then begin
  xr = (posfr eq 1) ? [-maxf, maxf] : [maxf, -maxf]
endif else if keyword_set(full) then begin
  xr = [(f[0] - 0.5*posfr*df), (f[ndata-1] + 0.5*posfr*df)]
endif else if posfr eq 1 then begin
  minf = -2000.0 > (f[0] - 0.5*df)
  maxf =  2000.0 < (f[ndata-1] + 0.5*df)
  xr = [minf, maxf]
endif else begin
  minf = -2000.0 > (f[ndata-1] - 0.5*df)
  maxf =  2000.0 < (f[0] + 0.5*df)
  xr = [maxf, minf]
endelse
bin1 = 0L > (bin0 + posfr*round(xr[0]/df))
bin2 = (ndata - 1) < (bin0 + posfr*round(xr[1]/df))

if total(plotpol) eq 2 then begin
  plottedvalues = [oc2[bin1:bin2], sc2[bin1:bin2]]
endif else if chan eq 1 then begin
  plottedvalues = oc2[bin1:bin2]
endif else begin
  plottedvalues = sc2[bin1:bin2]
endelse
maxy = max(plottedvalues)
miny = min(plottedvalues)
yspan = maxy - miny
maxy = maxy + 0.12*yspan
miny = miny - 0.08*yspan

; Get the y-limits for the +/- one-standard-deviation line

if keyword_set(rcs) then begin
  oc_sd_miny = (-oc_sdev/df) > miny
  oc_sd_maxy = ( oc_sdev/df) < maxy
  sc_sd_miny = (-sc_sdev/df) > miny
  sc_sd_maxy = ( sc_sdev/df) < maxy
endif else begin
  oc_sd_miny = -1.0 > miny
  oc_sd_maxy =  1.0 < maxy
  sc_sd_miny = -1.0 > miny
  sc_sd_maxy =  1.0 < maxy
endelse

; Get the maximum and minimum frequencies (for plotting the baseline).

if posfr eq 1 then begin
  fmin = f[0] - 0.5*df
  fmax = f[ndata-1] + 0.5*df
endif else begin
  fmin = f[ndata-1] - 0.5*df
  fmax = f[0] + 0.5*df
endelse

; Create the spectral plot(s)

xlabel = xr[0] + 0.9*(xr[1] - xr[0])
ylabel = miny + 0.85*(maxy - miny)
if nplots eq 1 then begin
  axis_charsize = 1.2
  label_charsize = 2
endif else if nplots eq 2 then begin
  axis_charsize = 1.1
  label_charsize = 1.5
endif else begin
  axis_charsize = 2
  label_charsize = 1.2
endelse
xtitlestring = strarr(3)
xtitlestring[nplots-1] = 'Doppler frequency  (Hz)'
if keyword_set(rcs) then begin
  ytitlestring = 'cross section density  (km^2 / Hz)'
endif else begin
  ytitlestring = 'noise standard deviations'
endelse
if plotpol[0] then begin
  plot,f,oc2,xrange=xr,yrange=[miny,maxy], $
       xtitle=xtitlestring[0],ytitle=ytitlestring, $
       charsize=axis_charsize,_extra=_ext
  if miny lt 0 and maxy gt 0 then oplot,[fmin,fmax],[0,0]
  if (miny lt -1 and maxy ge 0) or (miny le 0 and maxy gt 1) then begin
    oplot,[0,0],[oc_sd_miny,oc_sd_maxy],thick=2
  endif
  if (keyword_set(separate) or not plotpol[1]) then begin
    xyouts,xlabel,ylabel,'OC',alignment=0.5,charsize=label_charsize
  endif else begin
    oplot,f,sc2,linestyle=1
    xyouts,xlabel,ylabel,'OC, SC',alignment=0.5,charsize=label_charsize
  endelse
endif
if plotpol[1] and (keyword_set(separate) or not plotpol[0]) then begin
  plot,f,sc2,xrange=xr,yrange=[miny,maxy], $
       xtitle=xtitlestring[1],ytitle=ytitlestring, $
       charsize=axis_charsize,_extra=_ext
  if miny lt 0 and maxy gt 0 then oplot,[fmin,fmax],[0,0]
  if (miny lt -1 and maxy ge 0) or (miny le 0 and maxy gt 1) then begin
    oplot,[0,0],[sc_sd_miny,sc_sd_maxy],thick=2
  endif
  xyouts,xlabel,ylabel,'SC',alignment=0.5,charsize=label_charsize
endif

; Create the polarization ratio plot if requested

if domu then begin

  ;  lim1 = tags[0].jsnr1
  ;  lim2 = tags[0].jsnr2
  ;  ratio = smooth(ratio[lim1:lim2],10)
  ;  fmu_plot = f[lim1:lim2]

  if n_elements(maxmu) gt 0 then begin
    maxy = maxmu
    miny = -0.08*maxmu
  endif else begin
    plottedvalues = ratio[bin1:bin2]
    maxy = max(plottedvalues) < 1.2
    miny = 0
    yspan = maxy - miny
    maxy = maxy + 0.12*yspan
    miny = miny - 0.08*yspan
  endelse
  plot,f,ratio,xrange=xr,yrange=[miny,maxy], $
       xtitle=xtitlestring[2],charsize=axis_charsize,_extra=_ext
  oplot,[fmin,fmax],[0,0],linestyle=1
  ylabel = miny + 0.85*(maxy - miny)
  xyouts,xlabel,ylabel,'SC/OC',alignment=0.5,charsize=label_charsize
endif

if !d.window ge 0 then wshow  ; temporary fix: show the window now that the plot is completed

!p.multi = pmultisave

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro plotsum1,smooth=smooth,xrange=xrange,maxf=maxf,full=full,rcs=rcs,help=help,_extra=_ext

; Plot a single-channel spectrum
;
; Modified 2006 Jun 25 by CM: Frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded

; Check that frequency range was correctly input (if input at all)

if n_elements(xrange) gt 0 then begin
  if n_elements(maxf) gt 0 or keyword_set(full) or n_elements(xrange) ne 2 then begin
    xrangeOK = 0
  endif else begin
    xrangeOK = (xrange[0] lt xrange[1])
  endelse
endif else if n_elements(maxf) gt 0 then begin
  if keyword_set(full) or n_elements(maxf) ne 1 then begin
    xrangeOK = 0
  endif else begin
    xrangeOK = (maxf gt 0)
  endelse
endif else begin
  xrangeOK = 1
endelse

if n_params() ne 0 or not xrangeOK or keyword_set(help) then begin
  print,' '
  print,'plotsum1[,smooth=smooth][,xrange=[minf,maxf]][,maxf=maxf][,/full][,/rcs]'
  print,'        [,plot keywords][,/help]'
  print,' '
  print,'        Default frequency range is [-2000,2000] Hz or else the full spectrum,'
  print,'            whichever is narrower'
  print,' '
  print,'        The xrange, maxf, and full keywords are mutually exclusive:'
  print,' '
  print,'            For xrange, must have maxf > minf even if', $
                         ' frequency increases leftward',format='(2a)'
  print,'            maxf=100 is the same as xrange=[-100,100]'
  print,'            /full plots the full spectrum'
  print,' '
  print,'        /rcs expresses signal strength as radar cross section density in'
  print,'             km^2 per Hz; otherwise the unit is the rms noise'
  print,' '
  return
endif else if (*loaded1).ndata le 2 then begin
  print,"ERROR in plotsum1: There's no loaded single-channel spectrum to plot"
  return
endif

; Set up the plotting space

nplots = 1
pmultisave = !p.multi
!p.multi = [0,1,nplots]
erase

; Get some elements of the loaded single-channel spectrum

f = (*loaded1).freq
spec = (*loaded1).spec
tags1 = (*loaded1).tags
sdev = tags1.sdev
df = tags1.dfreq
posfr = tags1.posfr
bin0 = tags1.xjcen
ndata = (*loaded1).ndata
chanstring = (*loaded1).pol

; Express in absolute units if desired

if keyword_set(rcs) then spec = spec * sdev / df

; Smooth the spectra if desired

spec2 = spec
if n_elements(smooth) gt 0 then begin
  if smooth ge 2 then begin
    spec2 = smooth(spec2,smooth)
    if not keyword_set(rcs) then spec2 = spec2*sqrt(smooth)
  endif
endif

; Get the y-limits for the spectral plot(s)

if n_elements(xrange) gt 0 then begin
  xr = (posfr eq 1) ? xrange : reverse(xrange)
endif else if n_elements(maxf) gt 0 then begin
  xr = (posfr eq 1) ? [-maxf, maxf] : [maxf, -maxf]
endif else if keyword_set(full) then begin
  xr = [(f[0] - 0.5*posfr*df), (f[ndata-1] + 0.5*posfr*df)]
endif else if posfr eq 1 then begin
  minf = -2000.0 > (f[0] - 0.5*df)
  maxf =  2000.0 < (f[ndata-1] + 0.5*df)
  xr = [minf, maxf]
endif else begin
  minf = -2000.0 > (f[ndata-1] - 0.5*df)
  maxf =  2000.0 < (f[0] + 0.5*df)
  xr = [maxf, minf]
endelse
bin1 = 0L > (bin0 + posfr*round(xr[0]/df))
bin2 = (ndata - 1) < (bin0 + posfr*round(xr[1]/df))

plottedvalues = spec2[bin1:bin2]
maxy = max(plottedvalues)
miny = min(plottedvalues)
yspan = maxy - miny
maxy = maxy + 0.12*yspan
miny = miny - 0.08*yspan

; Get the y-limits for the +/- one-standard-deviation line

if keyword_set(rcs) then begin
  sd_miny = (-sdev/df) > miny
  sd_maxy = ( sdev/df) < maxy
endif else begin
  sd_miny = -1.0 > miny
  sd_maxy =  1.0 < maxy
endelse

; Get the maximum and minimum frequencies (for plotting the baseline).

if posfr eq 1 then begin
  fmin = f[0] - 0.5*df
  fmax = f[ndata-1] + 0.5*df
endif else begin
  fmin = f[ndata-1] - 0.5*df
  fmax = f[0] + 0.5*df
endelse

; Create the spectral plot

xlabel = xr[0] + 0.9*(xr[1] - xr[0])
ylabel = miny + 0.85*(maxy - miny)
axis_charsize = 1.2
label_charsize = 2
if keyword_set(rcs) then begin
  ytitlestring = 'cross section density  (km^2 / Hz)'
endif else begin
  ytitlestring = 'noise standard deviations'
endelse
plot,f,spec2,xrange=xr,yrange=[miny,maxy], $
     xtitle='Doppler frequency  (Hz)',ytitle=ytitlestring, $
     charsize=axis_charsize,_extra=_ext
if miny lt 0 and maxy gt 0 then oplot,[fmin,fmax],[0,0]
if (miny lt -1 and maxy ge 0) or (miny le 0 and maxy gt 1) then begin
  oplot,[0,0],[sd_miny,sd_maxy],thick=2
endif
xyouts,xlabel,ylabel,chanstring,alignment=0.5,charsize=label_charsize

!p.multi = pmultisave

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showim,bins=bins,window_scale=window_scale,aspect=aspect,interp=interp, $
           doprange=doprange,delrange=delrange,levels=levels,sqroot=sqroot, $
           log=log,gamma=gamma,help=help,_extra=_ext

; Display the loaded image.
;
; This routine is a modification of Mike Nolan's "image4" procedure

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'showim[,/bins][,/window_scale][,/aspect][,/interp]              $'
  print,'      [,doprange=doprange][,delrange=delrange][,levels=levels]  $
  print,'      [,/sqroot][,/log][,gamma=gamma][,contour keywords][,/help]'
  print,' '
  print,'Display part or all of a delay-Doppler image:'
  print,' '
  print,'      /bins sets the units on the x and y axes to image bins (counting from 0);'
  print,'            otherwise the units are Doppler in Hz (x) and delay in usec (y),'
  print,'            measured relative to the ephemeris predictions contained in the'
  print,"            image's 'eph_col' and 'eph_row' tags."
  print,' '
  print,'      /window_scale scales the window size to the image size;'
  print,'            otherwise the image size is scaled to the window size.'
  print,'            This keyword is ignored when outputting to devices with'
  print,'            scalable pixels (e.g., PostScript).'
  print,' '
  print,"      /aspect retains the image's aspect ratio.  Square pixels are assumed."
  print,'            Note that setting /window_scale automatically retains the aspect ratio.'
  print,' '
  print,'      /interp causes bilinear interpolation to be used if the image is resized.'
  print,' '
  print,'      doprange=[150,250] displays only the Doppler range 150 Hz to 250 Hz'
  print,'            (or only Doppler bins 150-250 if /bins is set).'
  print,'      doprange=100 is the same as doprange=[-100,100].'
  print,' '
  print,'      delrange=[-20.5,-15.5] displays only the delay range -20.5 usec to -15.5 usec'
  print,"            (or produces an error if /bins is set, since bin numbers can't be < 0)."
  print,' '
  print,'      levels=[1.0,4.5] scales the display so that pixel values <= 1.0 map to black'
  print,'            and pixel values >= 4.5 map to bright white.'
  print,'            (default: minimum, maximum pixel values --> black, bright white)'
  print,'      levels=3.8 is the same as levels=[0.0,3.8].'
  print,' '
  print,'      The next three scaling keywords are mutually exclusive:
  print,'            /sqroot displays the square root of pixel values.'
  print,"                 (If levels aren't specified, pixels < 0.0 are first reset to 0.0.)"
  print,'            /log displays the base-10 logarithm of pixel values.'
  print,"                 (If levels aren't specified, pixels < 1.0 are first reset to 1.0.)"
  print,'            gamma is the gamma correction parameter'
  print,' '
  return
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in showim: No image is loaded, so nothing can be displayed'
  return
endif else if (keyword_set(sqroot) and (keyword_set(log) or n_elements(gamma) gt 0))  $
           or (keyword_set(log) and n_elements(gamma) gt 0) then begin
  print,'ERROR in showim: /sqroot, /log, and gamma=gamma are mutually exclusive keywords'
  return
endif

; Get some elements of the loaded image

image = (*loadedi).image
width = (*loadedi).width
height = (*loadedi).height

; Get the units on the two axes

if keyword_set(bins) then begin
  x = findgen(width)
  y = findgen(height)
  xtitlestring = 'Doppler bin number'
  ytitlestring = 'Delay bin number'
endif else begin
  eph_col = getextrai('eph_col')
  eph_row = getextrai('eph_row')
  dfreq = getextrai('fres')
  delayunit = getextrai('delayunit')
  baud = getextrai('baudlen')
  if isnull(delayunit) and notnull(baud) then begin
    spb = getextrai('samples_per_baud')
    if isnull(spb) then spb = 1L
    rows_per_baud = getextrai('rows_per_baud')
    if isnull(rows_per_baud) then rows_per_baud = spb
    delayunit = baud/rows_per_baud
  endif
  if notnull(eph_col) and notnull(eph_row) and notnull(dfreq) and notnull(delayunit) then begin
    x = dfreq*(findgen(width) - eph_col)       ; Doppler (Hz)
    y = delayunit*(findgen(height) - eph_row)  ; delay (usec)
    xtitlestring = 'Doppler (Hz)'
    ytitlestring = 'Delay (' + string("265B) + 's)'  ; includes lower-case "mu"
  endif else begin
    print,' '
    print,'WARNING in showim: Axis units will be bin number rather than Hz and usec'
    print,"                   -- need the 'eph_col' and 'eph_row' and 'fres' tags,"
    print,"                      plus either the 'delayunit' or 'baudlen' tag"
    print,' '
    x = findgen(width)
    y = findgen(height)
    xtitlestring = 'Doppler bin number'
    ytitlestring = 'Delay bin number'
  endelse
endelse
dx = x[1] - x[0]
dy = y[1] - y[0]
axis_charsize = 1.2

; Get the x-axis limits

if n_elements(doprange) gt 0 then begin
  if n_elements(doprange) ne 1 and n_elements(doprange) ne 2 then begin
    print,'ERROR in showim: doprange must be a scalar or a 2-element vector'
    return
  endif else if n_elements(doprange) eq 1 and doprange[0] le 0 then begin
    print,'ERROR in showim: scalar doprange must be > 0'
    return
  endif
  xrange = (n_elements(doprange) eq 1) ? [-doprange[0], doprange[0]] : doprange
  if xrange[1] le xrange[0] then begin
    print,'ERROR in showim: must have doprange[1] > doprange[0]'
    return
  endif else if xrange[0] gt max(x) or xrange[1] lt min(x) then begin
    print,"ERROR in showim: doprange is entirely outside the image's Doppler limits"
    return
  endif
  lx = max( [ max( where((x-0.5*dx) le xrange[0]) ) , 0L ] )
  ux = min( where((x+0.5*dx) ge xrange[1]) )
  if ux eq -1 then ux = width - 1
  if ux eq lx then begin
    print,'ERROR in showim: doprange includes only one Doppler column'
    return
  endif
  x = x[lx:ux]
  width_use = ux - lx + 1
endif else begin
  lx = 0L
  ux = width - 1
  width_use = width
endelse
x = [x, x[n_elements(x)-1]+dx] - 0.5*dx

; Get the y-axis limits

if n_elements(delrange) gt 0 then begin
  if n_elements(delrange) ne 2 then begin
    print,"ERROR in showim: delrange must have 2 elements"
    return
  endif else if delrange[1] le delrange[0] then begin
    print,"ERROR in showim: must have delrange[1] > delrange[0]"
    return
  endif else if delrange[0] gt max(y) or delrange[1] lt min(y) then begin
    print,"ERROR in showim: delrange is entirely outside the image's delay limits"
    return
  endif
  yrange = delrange
  ly = max( [ max( where((y-0.5*dy) le yrange[0]) ) , 0L ] )
  uy = min( where((y+0.5*dy) ge yrange[1]) )
  if uy eq -1 then uy = height - 1
  if uy eq ly then begin
    print,"ERROR in showim: delrange includes only one delay row"
    return
  endif
  y = y[ly:uy]
  height_use = uy - ly + 1
endif else begin
  ly = 0L
  uy = height - 1
  height_use = height
endelse
y = [y, y[n_elements(y)-1]+dy] - 0.5*dy

; Truncate the image in x and y

image = image[lx:ux,ly:uy]

; Rescale the image if specified

maxpix = max(image, min=minpix)
if maxpix eq minpix then begin
  print,'CANCEL DISPLAY in showim: All pixels to be displayed have the same value (',  $
        maxpix,')',format='(a,f0,a)'
  return
endif

if n_elements(levels) gt 0 then begin
  if n_elements(levels) eq 1 then begin
    blacklevel = 0.0
    whitelevel = 1.0*levels[0]
    if whitelevel le blacklevel then begin
      print,'ERROR in showim: levels must be > 0 if scalar'
      return
    endif
  endif else if n_elements(levels) eq 2 then begin
    blacklevel = 1.0*levels[0]
    whitelevel = 1.0*levels[1]
    if whitelevel le blacklevel then begin
      print,'ERROR in showim: must have levels[1] > levels[0]'
      return
    endif
  endif else begin
    print,'ERROR in showim: levels must be a scalar or a 2-element vector'
    return
  endelse
  if blacklevel ge maxpix or whitelevel le minpix then begin
    print,"ERROR in showim: levels is entirely outside the image's pixel-value range"
    return
  endif else if keyword_set(sqroot) and blacklevel lt 0.0 then begin
    print,'ERROR in showim: levels must be non-negative if /sqroot is set'
    return
  endif else if keyword_set(log) and blacklevel le 0.0 then begin
    print,'ERROR in showim: levels must be positive if /log is set'
    return
  endif
endif else if keyword_set(sqroot) then begin
  blacklevel = minpix > 0.0
  whitelevel = maxpix
  if whitelevel le blacklevel then begin
    print,'ERROR in showim: max pixel value <= 0.0, must use levels keyword or unset /sqroot'
    return
  endif
endif else if keyword_set(log) then begin
  blacklevel = minpix > 1.0
  whitelevel = maxpix
  if whitelevel le blacklevel then begin
    print,'ERROR in showim: max pixel value <= 1.0, must use levels keyword or unset /log'
    return
  endif
endif else begin
  blacklevel = minpix
  whitelevel = maxpix
endelse

image = (image > blacklevel) < whitelevel
if keyword_set(sqroot) then begin
  image = sqrt(image)
  blacklevel = sqrt(blacklevel)
  whitelevel = sqrt(whitelevel)
endif else if keyword_set(log) then begin
  image = alog10(image)
  blacklevel = alog10(blacklevel)
  whitelevel = alog10(whitelevel)
endif else if n_elements(gamma) gt 0 then begin
  if gamma le 0.0 then begin
    print,'ERROR in showim: gamma must be positive'
    return
  endif
  image = (whitelevel - blacklevel) *                                     $
          ( (image - blacklevel)/(whitelevel - blacklevel) )^(1.0/gamma)  $
          + blacklevel
endif
byte_image = bytscl(image, min=blacklevel, max=whitelevel, top=!d.table_size-1)

; Print the minimum and maximum (pre-scaled) pixel values to be displayed

if abs(minpix) le 999999L and abs(maxpix) le 999999L then begin
  formatstring = '(f9.2)'
endif else begin
  formatstring = '(e12.5)'
endelse
print,'Minimum and maximum displayed pixel values = ',         $
      strtrim(string(minpix,format=formatstring), 2),' and ',  $
      strtrim(string(maxpix,format=formatstring), 2),          $
      format='(4a)'

; Set the window used by contour

contour,[[0,0],[1,1]],/nodata,xstyle=4,ystyle=4,xmargin=[12,3],  $
        charsize=axis_charsize,_extra=_ext

; Get the image size and the initial values for the window size

px = !x.window * !d.x_vsize     ; Get size of window in device units
py = !y.window * !d.y_vsize
swx = px[1] - px[0]             ; Window size in x in device units
swy = py[1] - py[0]             ; Window size in y in device units
six = float(width_use)          ; Image sizes
siy = float(height_use)
aspi = six/siy                  ; Image aspect ratio
aspw = swx/swy                  ; Window aspect ratio
f = aspi/aspw                   ; Ratio of aspect ratios

; Output the image; the method depends on whether or not pixels are scalable

if (!d.flags and 1) ne 0 then begin

  ; Scalable pixels (e.g., PostScript)

  ; If /aspect is set, adjust window size to retain image aspect ratio

  if keyword_set(aspect) then begin
    if f ge 1.0 then swy = swy/f else swx = swx*f
  endif
    
  ; Output image

  tv,byte_image,px[0],py[0],xsize=swx,ysize=swy,/device

endif else begin

  ; Not scalable pixels

  if keyword_set(window_scale) then begin

    ; Output image, then set window size equal to image size

    tv,byte_image,px[0],py[0]
    swx = six
    swy = siy

  endif else begin

    ; If /aspect is set, adjust window size to retain image aspect ratio

    if keyword_set(aspect) then begin
      if f ge 1.0 then swy = swy/f else swx = swx*f
    endif

    ; Output resampled image

    tv,poly_2d(byte_image, [[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],  $
                           keyword_set(interp),swx,swy),              $
       px[0],py[0]

  endelse

endelse

; Do the contour

mx = !d.n_colors - 1            ; Brightest color
colors = [mx,mx,mx,0,0,0]       ; Color vectors
contour,intarr(width_use+1,height_use+1,/noz),x,y,/noerase,/nodata,/xst,/yst, $
        pos=[px[0],py[0],px[0]+swx,py[0]+swy],/dev,c_color=colors,            $
        xtitle=xtitlestring,ytitle=ytitlestring,charsize=axis_charsize,       $
        xmargin=[12,3],_extra=_ext
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro impix,bins=bins,help=help

; Continuously prints pixel values for the displayed loaded image.

common loadedBlock,loadedi,loaded1,loaded

if keyword_set(help) or n_params() gt 0 then begin
  print,' '
  print,'impix[,/bins][,/help]'
  print,' '
  print,'Left-click to begin, then continuously print pixel values'
  print,'      for the displayed loaded image as you move the cursor;'
  print,'      click the middle or right mouse button to quit.'
  print,' '
  print,"Set /bins IFF you used 'showim,/bins' to display the image'
  print,' '
  print,"NOTE: If you don't (un)set /bins correctly, or if the current'
  print,"      loaded image isn't the same as the displayed image,"
  print,'      impix will return incorrect pixel values!'
  print,' '
  return
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in showdimsi: No image is loaded'
  return
endif

; Get the image and the Doppler and delay vectors

image = (*loadedi).image
width = (*loadedi).width
height = (*loadedi).height
dfreq = getextrai('fres')
eph_col = getextrai('eph_col')
eph_row = getextrai('eph_row')
delayunit = getextrai('delayunit')
spb = getextrai('samples_per_baud')
rows_per_baud = getextrai('rows_per_baud')
baud = getextrai('baudlen')
if isnull(delayunit) and notnull(baud) then begin
  if isnull(spb) then spb = 1L
  if isnull(rows_per_baud) then rows_per_baud = spb
  delayunit = baud/rows_per_baud
endif
if notnull(eph_col) and notnull(eph_row) and notnull(dfreq) and notnull(delayunit) then begin
  f = dfreq*(findgen(width) - eph_col)       ; Doppler (Hz)
  d = delayunit*(findgen(height) - eph_row)  ; delay (usec)
endif else begin
  print,"ERROR in impix: Need the 'eph_col' and 'eph_row' and 'fres' tags,"
  print,"                plus either the 'delayunit' or 'baudlen' tag"
  return
endelse

; Get some format strings for nice output

maxpix = max(image, min=minpix)
len = strlen(strtrim(string(maxpix,format='(f9.2)'), 2))    $
      > strlen(strtrim(string(minpix,format='(f9.2)'), 2))
formatstring1 = '(f' + string(len,format='(i0)') + '.2)'
len = strlen(string(width-1,format='(i0)'))
formatstring2 = '(i' + string(len,format='(i0)') + ')'
len = strlen(string(height-1,format='(i0)'))
formatstring3 = '(i' + string(len,format='(i0)') + ')'
len = strlen(strtrim(string(f[0],format='(f9.2)'), 2))          $
      > strlen(strtrim(string(f[width-1],format='(f9.2)'), 2))
formatstring4 = '(f' + string(len,format='(i0)') + '.2)'
len = strlen(strtrim(string(d[0],format='(f9.2)'), 2))          $
      > strlen(strtrim(string(d[height-1],format='(f9.2)'), 2))
formatstring5 = '(f' + string(len,format='(i0)') + '.2)'

; Display pixel values in a loop until the middle or right mouse button is clicked

print,' '
print,'Left-click on a pixel to begin,',                           $
      ' then move the cursor to print pixel values.',format='(2a)'
print,'Right-click or middle-click to end.'
print,' '
dopbin_last = -1
delbin_last = -1
cursor,x,y,/down
while (!mouse.button le 1) do begin
  if keyword_set(bins) then begin
    dopbin = round(x)
    delbin = round(y)
  endif else begin
    dopbin = round(x/dfreq + eph_col)
    delbin = round(y/delayunit + eph_row)
  endelse
  if dopbin ge 0 and dopbin lt width and delbin ge 0 and delbin lt height  $
                 and (dopbin ne dopbin_last or delbin ne delbin_last)      $
                 then begin
    print,'Pixel value = ',                                             $
          string(image[dopbin,delbin],format=formatstring1),'  at  (',  $
          string(dopbin,format=formatstring2),', ',                     $
          string(delbin,format=formatstring3),')  (',                   $
          string(f[dopbin],format=formatstring4),' Hz, ',               $
          string(d[delbin],format=formatstring5),' us)',                $
          format='(11a)'
    dopbin_last = dopbin
    delbin_last = delbin
  endif
  cursor,x,y,/change
endwhile

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writepds,outfile,chan=chan,tkplay=tkplay,help=help,_extra=_ext

; Write the loaded pair to disk as an pds file

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 1 or keyword_set(help) then begin
  print,'writepds,outfile[,/overwrite][,/append][,chan=1 or 2][,/tkplay][,/help]'
  print,"         '.csv' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writepds,outfile[,/overwrite][,/append][,chan=1 or 2][,/tkplay][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loaded).ndata le 2 then begin
  print,'ERROR in writepds: No pair is loaded, so nothing can be written to disk'
  return
endif

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

; Get some elements of the loaded spectrum

pair = (*loaded).spec
extratags = (*loaded).extratags
tags = (*loaded).tags
ndata = (*loaded).ndata
ntags = (*loaded).ntags
tname = (*loaded).tname
nextra = (*loaded).nextra
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
    print,'ERROR in writepds: jcp tags are already [0,1] rather than [1,2]'
    return
  endif
  extratags_write = [extratags_write, extratags_write[0]]
  extratags_write[nextra_write].format = ''
  extratags_write[nextra_write].name = ''
  extratags_write[nextra_write].value = ''
  extratags_write[nextra_write].comment = $
               '# Conv to   tkplay format ' + systime(/utc) + ' UT'
  nextra_write = nextra_write + 1L
endif

; Open the output file

err = openOutfile(lun,outfile,'csv',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the pds header

sformat = '(a,5x,a16,3x,a)'
iformat = '(a,5x,a16,3x,i0)'
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
    for n=0L,ntags-1 do begin
      writeu,lun,swap_endian(float(tags[ch-1].(n)), /swap_if_little_endian)
    endfor
  endif
endfor

; Write the pds footer

for n=0L,nextra_write-1 do begin
  printf,lun,strtrim(string(extratags_write[n],format='(a,5x,a16,3x,a,3x,a)'), 2)
endfor
printf,lun,'s','target',tname,format=sformat
chanstring = ['OC','SC']
if nchan eq 1 then printf,lun,'s','polarization',chanstring[chan-1],format=sformat
printf,lun,'.'

if nchan eq 2 then begin
  print,'Wrote 1 pds frame, containing 1 OC and 1 SC spectrum, to file ',outfile
endif else begin
  print,'Wrote 1 pds frame, containing 1 ',chanstring[chan-1], $
        ' spectrum, to file ',outfile,format='(4a)'
endelse

close,lun
free_lun,lun

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
