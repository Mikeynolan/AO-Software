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
common channelBlock, chanstrings, maxchan

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
  print,'       /separate produces separate plots for each pol; otherwise the'
  print,'            spectra are superimposed within the same plot'
  print,' '
  print,'       /rcs expresses signal strength as radar cross section density in'
  print,'            km^2 per Hz; otherwise the unit is the rms noise'
  print,'       chan can be a single number or a vector, for multiple superimposed plots.'
  print,' '
  return
endif else if (*loaded).ndata le 2 then begin
  print,"ERROR in plotsum: There's no loaded pair to plot"
  return
endif

; Check which channel(s) to write

; Get some elements of the loaded pair

f = (*loaded).freq
specs = (*loaded).spec
oc = reform((*loaded).spec[0,*])
sc = reform((*loaded).spec[1,*])
tags = (*loaded).tags
npol = n_elements(tags)
sdev = tags.sdev
oc_sdev = tags[0].sdev
sc_sdev = tags[1].sdev
xx_sdev = sqrt(sdev[0] * sdev[1])
df = tags[0].dfreq
posfr = tags[0].posfr
bin0 = tags[0].xjcen
ndata = (*loaded).ndata
ocsconly = 1 ; Flag if anything but OC or SC is requested

plotpol = intarr(npol)
if n_elements(chan) eq 0 then begin ; default is OC,SC
  plotpol[0:1] = 1
endif else begin
  if n_elements(chan) eq 1 then begin
    if chan ge 1 && chan le npol then begin
      plotpol[chan-1] = 1
      if chan gt 2 then ocsconly = 0
    endif else begin
      print, 'ERROR in plotsum: Chan must be 1 .. npol'
      return
    endelse
  endif else begin
    for i = 0, n_elements(chan)-1 do begin
      if chan[i] ge 1 && chan[i] le npol then begin
        plotpol[chan[i]-1] = 1
        if chan[i] gt 2 then ocsconly = 0
      endif else begin
        print, 'ERROR in plotsum: All chans must be 1 .. npol'
        return
      endelse
    endfor
  endelse
endelse

nplots = (keyword_set(separate)) ? total(plotpol) : 1L

; Express in absolute units if desired

; For stokes, need rcs
if ~ ocsconly then rcs = 1

if keyword_set(rcs) then begin
  specrcs = specs * (sdev # replicate(1,ndata)) * df
  oc = oc * oc_sdev / df
  sc = sc * sc_sdev / df
  specs = specrcs
endif

; Smooth the spectra if desired

oc2 = oc
sc2 = sc
if n_elements(smooth) gt 0 then begin
  specsm = specs
  if smooth ge 2 then begin
    for i=0, npol-1 do begin
      specsm[i,*] = smooth(specs[i,*], smooth)
    endfor
    oc2 = smooth(oc2,smooth)
    sc2 = smooth(sc2,smooth)
    if not keyword_set(rcs) then begin
      oc2 = oc2*sqrt(smooth)
      sc2 = sc2*sqrt(smooth)
    endif
  endif
  specs = specsm
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

if total(plotpol) lt 1 then begin
  printf, 'No channels to plot'
  return
endif
plottedvalues = []
for i = 0, npol-1 do begin
  if plotpol[i] then plottedvalues = [plottedvalues, reform(specs[i,bin1:bin2])]
endfor
maxy = max(plottedvalues)
miny = min(plottedvalues)
yspan = maxy - miny
maxy = maxy + 0.12*yspan
miny = miny - 0.08*yspan

; Get the y-limits for the +/- one-standard-deviation line

sdmins = intarr(npol)
sdmaxs = intarr(npol)
for i = 0, npol-1 do begin
  if keyword_set(rcs) then begin
    sdmins[i] = (-sdev[i]/df) > miny
    sdmaxs[i] = ( sdev[i]/df) < maxy
  endif else begin
    sdmins[i] =  -1 > miny
    sdmaxs[i] =   1 > maxy
  endelse
endfor
    
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
xtitlestring = strarr(nplots)
xtitlestring[nplots-1] = 'Doppler frequency  (Hz)'
if keyword_set(rcs) then begin
  ytitlestring = 'cross section density  (km^2 / Hz)'
endif else begin
  ytitlestring = 'noise standard deviations'
endelse

firstplot = 1
polstrings = ""
for i = 0, npol-1 do begin
  if plotpol[i] then begin
    polstrings = polstrings + (firstplot ? '' : ', ') + chanstrings[i]
    firstplot = 0
  endif
endfor
  
firstplot = 1
ls = 0
for i= 0, npol-1 do begin
  if plotpol[i] then begin
    if firstplot || keyword_set(separate) then begin
      firstplot = 0
      plot,f,specs[i,*],xrange=xr,yrange=[miny,maxy], $
           xtitle=xtitlestring[0],ytitle=ytitlestring, $
           charsize=axis_charsize,_extra=_ext
      if miny lt 0 and maxy gt 0 then oplot,[fmin,fmax],[0,0]
      if (miny lt -1 and maxy ge 0) or (miny le 0 and maxy gt 1) then begin
        oplot,[0,0],[sdmins[i], sdmaxs[i]],thick=2
      endif 
      if keyword_set(separate) then begin
        xyouts,xlabel,ylabel,chanstrings[i],alignment=0.5,charsize=label_charsize
      endif else begin
        xyouts,xlabel,ylabel,polstrings,alignment=0.5,charsize=label_charsize
      endelse
    endif else begin
      oplot, f, specs[i,*], linestyle=ls
    endelse
    ls = ls + 1
  endif
endfor

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
       xtitle=xtitlestring[nplots-1],charsize=axis_charsize,_extra=_ext
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

; For stokes, need rcs
if keyword_set(chan) and chan lt 0 then rcs = 1

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
  print,'      The next three scaling keywords are mutually exclusive:'
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
    ytitlestring = 'Delay (' + string("265B) + 's)'  ; includes lower-case "mu" "
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
pro writerdf,outfile,chan=chan,tkplay=tkplay,help=help,_extra=_ext

; Write the loaded pair to disk as an rdf file

common loadedBlock,loadedi,loaded1,loaded
common channelBlock, chanstrings, maxchan

if n_params() ne 1 or keyword_set(help) then begin
  print,'writerdf,outfile[,/overwrite][,/append][,chan=1 or 2][,/tkplay][,/help]'
  print,"         '.rdf' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writerdf,outfile[,/overwrite][,/append][,chan=1 or 2][,/tkplay][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loaded).ndata le 2 then begin
  print,'ERROR in writerdf: No pair is loaded, so nothing can be written to disk'
  return
endif


; Get some elements of the loaded spectrum

pair = (*loaded).spec
extratags = (*loaded).extratags
tags = (*loaded).tags
ndata = (*loaded).ndata
ntags = (*loaded).ntags
tname = (*loaded).tname
nextra = (*loaded).nextra
npol = n_elements((*loaded).tags)
; Check which channel(s) to write

writepol = intarr(npol)
if n_elements(chan) eq 0 then begin
  writepol = writepol + 1
endif else if chan ge 1 && chan le npol then begin
  writepol[chan-1] = 1
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
  return
endelse
nchan = total(writepol)
if nchan gt 1 then begin
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
    print,'ERROR in writerdf: jcp tags are already [0,1] rather than [1,2]'
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

err = openOutfile(lun,outfile,'rdf',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the rdf header

sformat = '(a,5x,a16,3x,a)'
iformat = '(a,5x,a16,3x,i0)'
printf,lun,'s','type','float',format=sformat
printf,lun,'i','ndata',ndata,format=iformat
printf,lun,'i','height',nchan,format=iformat
printf,lun,'i','width',ndata+ntags,format=iformat
printf,lun,'i','size',4,format=iformat
printf,lun,'s','machine','SPARC',format=sformat ; === bigendian IEEE
printf,lun,'i','nchan',nchan,format=iformat
printf,lun,'i','ntags',ntags,format=iformat
printf,lun,'i','nspec',nchan,format=iformat
printf,lun,'s','format','RDF_spectrum_1',format=sformat
printf,lun,'.'

; Write the spectrum and tags (as binary)
; -- write all tags as floating point to be compatible with tkplay

for ch=1,npol do begin
  if writepol[ch-1] then begin
    writeu,lun,swap_endian(pair[ch-1,*], /swap_if_little_endian)
    for n=0L,ntags-1 do begin
      writeu,lun,swap_endian(float(tags[ch-1].(n)), /swap_if_little_endian)
    endfor
  endif
endfor

; Write the rdf footer

for n=0L,nextra_write-1 do begin
  printf,lun,strtrim(string(extratags_write[n],format='(a,5x,a16,3x,a,3x,a)'), 2)
endfor
printf,lun,'s','target',tname,format=sformat
polstring = ""
first=1
for n=0, npol-1 do begin
  if (writepol[n]) then polstring = polstring + (first ? "" : ", ") + chanstrings[n]
  first = 0
endfor
 
if nchan eq 1 then printf,lun,'s','polarization',chanstrings[chan-1],format=sformat
printf,lun,'.'

if nchan ge 2 then begin
  print,'Wrote 1 rdf frame, containing '+polstring+' spectra, to file ',outfile
endif else begin
  print,'Wrote 1 rdf frame, containing 1 ',chanstring[chan-1], $
        ' spectrum, to file ',outfile,format='(4a)'
endelse

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writerdf1,outfile,tkplay=tkplay,help=help,_extra=_ext

; Write the loaded single-channel spectrum to disk as an rdf file

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 1 or keyword_set(help) then begin
  print,'writerdf1,outfile[,/overwrite][,/append][,/tkplay][,/help]'
  print,"          '.rdf' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writerdf1,outfile[,/overwrite][,/append][,/tkplay][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loaded1).ndata le 2 then begin
  print,'ERROR in writerdf1: No single-channel spectrum is loaded,', $
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
    print,'ERROR in writerdf1: jcp tag for OC is already 0 rather than 1'
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

err = openOutfile(lun,outfile,'rdf',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the rdf header

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

; Write the rdf footer

for n=0L,nextra-1 do begin
  printf,lun,strtrim(string(extratags[n],format='(a,5x,a16,3x,a,3x,a)'), 2)
endfor
printf,lun,'s','target',tname,format=sformat
printf,lun,'s','polarization',pol,format=sformat
printf,lun,'.'

print,'Wrote 1 rdf frame, containing 1 ',pol, $
      ' spectrum, to file ',outfile,format='(4a)'

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writerdfi,outfile,help=help,_extra=_ext

; Write the loaded image to disk as an rdf file

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 1 or keyword_set(help) then begin
  print,'writerdfi,outfile[,/overwrite][,/append][,/help]'
  print,"          '.rdf' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writerdfi,outfile[,/overwrite][,/append][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in writerdfi: No image is loaded,', $
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

; Open the output file

err = openOutfile(lun,outfile,'rdf',/get_lun,_extra=_ext)
if err ne 0 then return

; Write the rdf header

sformat = '(a,5x,a16,3x,a)'
iformat = '(a,5x,a16,3x,i0)'
printf,lun,'s','type','float',format=sformat
printf,lun,'i','height',height,format=iformat
printf,lun,'i','width',width,format=iformat
printf,lun,'i','size',4,format=iformat
printf,lun,'s','machine','SPARC',format=sformat
printf,lun,'s','format','RDF_range_1',format=sformat
printf,lun,'.'

; Write the image (as binary)

linefeed = 10B
writeu,lun,swap_endian(image, /swap_if_little_endian),   $
           swap_endian(linefeed, /swap_if_little_endian)

; Write the rdf footer

for n=0L,nextra-1 do begin
  printf,lun,strtrim(string(extratags[n],format='(a,5x,a16,3x,a,3x,a)'), 2)
endfor
printf,lun,'s','target',tname,format=sformat
printf,lun,'s','polarization',pol,format=sformat
printf,lun,'.'

print,'Wrote 1 rdf frame, containing 1 ',pol, $
      ' image, to file ',outfile,format='(4a)'

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writestackrdf,outfile,group=g,ming=ming,maxg=maxg,arrayg=ag, $
                  chan=chan,tkplay=tkplay,help=help,_extra=_ext

; Take stack pairs in one or more groups and write them to disk as an rdf file,
; one rdf frame per pair.  The keywords g, ming and maxg, or ag specify the groups;
; omitting all of these causes all stack pairs to be written.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

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
  print,'writestackrdf,outfile[,/overwrite][,/append]'
  print,'              [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'              [,chan=1 or 2][,/tkplay][,/help]'
  print,' '
  print,"    '.rdf' output file extension is added if not already present"
  print,' '
  print,'    Setting keywords g, ming and maxg, or ag writes all pairs in ', $
        'the specified group(s) to an rdf file',format='(2a)'
  print,'    -- these three keyword choices are mutually exclusive'
  print,'    Calling writestackrdf with none of those keywords set writes ', $
        'all stack pairs to an rdf file',format='(2a)'
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writestackrdf,outfile[,/overwrite][,/append]'
  print,'              [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'              [,chan=1 or 2][,/tkplay][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if nstack eq 0 then begin
  print,'ERROR in writestackrdf: The pair stack is empty'
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
    print,'ERROR in writestackrdf: Must specify range with ming <= maxg'
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
    print,'writestackrdf: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
endif else if n_elements(ag) gt 0 then begin
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroups then begin
    print,"ERROR in writestackrdf: Don't include duplicate group numbers ", $
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
    print,'ERROR in writestackrdf: Group #',groupvec[j], $
          ' is not present in the pair stack',format='(a,i0,a)'
    return
  endif
endfor

; Open the output file

err = openOutfile(lun,outfile,'rdf',/get_lun,_extra=_ext)
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
    npol = n_elements(tags)
    writepol = intarr(npol)
    if n_elements(chan) eq 0 then begin
      writepol = writepol + 1
    endif else if chan ge 1 && chan le npol then begin
      writepol[chan-1] = 1
    endif else begin
      print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
      return
    endelse
    nchan = total(writepol)
    if nchan ge 2 then begin
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
        print,'ERROR in writestackrdf: jcp tags are already [0,1] rather than [1,2]'
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

    ; Write the rdf header

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

    for ch=1,npol do begin
      if writepol[ch-1] then begin
        writeu,lun,swap_endian(pair[ch-1,*], /swap_if_little_endian)
        for k=0L,ntags-1 do begin
          writeu,lun,swap_endian(float(tags[ch-1].(k)), /swap_if_little_endian)
        endfor
      endif
    endfor

    ; Write the rdf footer

    for k=0L,nextra_write-1 do begin
      printf,lun,strtrim(string(extratags_write[k],format='(a,5x,a16,3x,a,3x,a)'), 2)
    endfor
    printf,lun,'s','target',tname,format=sformat
    polstring = ""
    first=1
    for ch=0, npol-1 do begin
      if (writepol[ch]) then polstring = polstring + (first?"":", ") + chanstrings[ch]
      first = 0
    endfor

    if nchan eq 1 then printf,lun,'s','polarization',chanstrings[chan-1],format=sformat
    printf,lun,'.'

    nwrite = nwrite + 1L

  endif
endfor

if nchan ge 2 then begin
  print,'Wrote ',nwrite,' rdf frames, the last containing '+polstring+' spectra, to file ', $
        outfile,format='(a,i0,2a)'
endif else begin
  print,'Wrote ',nwrite,' rdf frames, each containing 1 ',chanstrings[chan-1], $
        ' spectrum, to file ',outfile,format='(a,i0,4a)'
endelse

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writestackrdf1,outfile,group=g,ming=ming,maxg=maxg,arrayg=ag, $
                   chan=chan,tkplay=tkplay,help=help,_extra=_ext

; Take single-channel stack spectra in one or more groups and write them to disk
; as an rdf file, one rdf frame per spectrum.  The keywords g, ming and maxg,
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
  print,'writestackrdf1,outfile[,/overwrite][,/append]'
  print,'               [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'               [,chan=1 or 2][,/tkplay][,/help]'
  print,' '
  print,"    '.rdf' output file extension is added if not already present"
  print,' '
  print,'    Setting keywords g, ming and maxg, or ag writes all single-channel ', $
        'spectra in the specified group(s) to an rdf file',format='(2a)'
  print,'    -- these three keyword choices are mutually exclusive'
  print,'    Calling writestackrdf1 with none of those keywords set writes ', $
        'all stack1 spectra to an rdf file',format='(2a)'
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writestackrdf1,outfile[,/overwrite][,/append]'
  print,'               [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'               [,chan=1 or 2][,/tkplay][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if nstack1 eq 0 then begin
  print,'ERROR in writestackrdf1: The single-channel stack is empty'
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
    print,'ERROR in writestackrdf1: Must specify range with ming <= maxg'
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
    print,'writestackrdf1: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
endif else if n_elements(ag) gt 0 then begin
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroups then begin
    print,"ERROR in writestackrdf1: Don't include duplicate group numbers ", $
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
    print,'ERROR in writestackrdf1: Group #',groupvec[j], $
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
  print,'writestackrdf1: No spectra meet these group/polarization criteria'
  return
endif

; Open the output file

err = openOutfile(lun,outfile,'rdf',/get_lun,_extra=_ext)
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
        print,'ERROR in writestackrdf1: jcp tag for OC is already 0 rather than 1'
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

    ; Write the rdf header

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

    ; Write the rdf footer

    for k=0L,nextra-1 do begin
      printf,lun,strtrim(string(extratags[k],format='(a,5x,a16,3x,a,3x,a)'), 2)
    endfor
    printf,lun,'s','target',tname,format=sformat
    printf,lun,'s','polarization',pol,format=sformat
    printf,lun,'.'

    nwrite = nwrite + 1L

  endif
endfor

print,'Wrote ',nwrite,' rdf frames, each containing 1 single-channel spectrum, to file ', $
      outfile,format='(a,i0,2a)'

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writestackrdfi,outfile,group=g,ming=ming,maxg=maxg,arrayg=ag, $
                   chan=chan,help=help,_extra=_ext

; Take stacki images in one or more groups and write them to disk
; as an rdf file, one rdf frame per spectrum.  The keywords g, ming and maxg,
; or ag specify the groups; omitting all of these causes all stacki images to 
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
  print,'writestackrdfi,outfile[,/overwrite][,/append]'
  print,'               [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'               [,chan=1 or 2][,/help]'
  print,' '
  print,"    '.rdf' output file extension is added if not already present"
  print,' '
  print,'    Setting keywords g, ming and maxg, or ag writes all ', $
        'images in the specified group(s) to an rdf file',format='(2a)'
  print,'    -- these three keyword choices are mutually exclusive'
  print,'    Calling writestackrdfi with none of those keywords set writes ', $
        'all stacki images to an rdf file',format='(2a)'
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writestackrdfi,outfile[,/overwrite][,/append]'
  print,'               [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'               [,chan=1 or 2][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if nstacki eq 0 then begin
  print,'ERROR in writestackrdfi: The image stack is empty'
  return
endif

; Create a sorted array of the group numbers which will be written

allgroups = stackgroupsi()
n_all = n_elements(allgroups)

if n_elements(g) gt 0 then begin
  ngroups = 1L
  groupvec = [long(g)]
endif else if n_elements(ming) gt 0 then begin
  if ming gt maxg then begin
    print,'ERROR in writestackrdfi: Must specify range with ming <= maxg'
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
    print,'writestackrdfi: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
endif else if n_elements(ag) gt 0 then begin
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroups then begin
    print,"ERROR in writestackrdfi: Don't include duplicate group numbers ", $
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
    print,'ERROR in writestackrdfi: Group #',groupvec[j], $
          ' is not present in the image stack',format='(a,i0,a)'
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
; one image will be written (not a foregone conclusion
; if the chan keyword is set)

n = 0L
count = 0L

while n lt nstacki and count eq 0 do begin

  ; See if stacki image #(n+1) will be written

  group_el = where( (groupvec eq (*stacki[n]).group) and                 $
                    ( (writepol[0] and (*stacki[n]).pol eq 'OC') or      $
                      (writepol[1] and (*stacki[n]).pol eq 'SC')     ),  $
                    count)
  n = n + 1
endwhile

if count eq 0 then begin
  print,'writestackrdfi: No images meet these group/polarization criteria'
  return
endif

; Open the output file

err = openOutfile(lun,outfile,'rdf',/get_lun,_extra=_ext)
if err ne 0 then return

; Go through the image stack, image by image and write out 
; images which are included in one of the specified groups

nwrite = 0L

for n=0L,nstacki-1 do begin

  ; See if this image is included in one of the groups which will be written

  group_el = where(groupvec eq (*stacki[n]).group, count)

  if count gt 0 then begin

    ; This image should be written to disk

    image = (*stacki[n]).image
    height = (*stacki[n]).height
    width = (*stacki[n]).width
    extratags = (*stacki[n]).extratags
    nextra = (*stacki[n]).nextra
    tname = (*stacki[n]).tname
    pol = (*stacki[n]).pol

    ; Write the rdf header

    sformat = '(a,5x,a16,3x,a)'
    iformat = '(a,5x,a16,3x,i0)'
    printf,lun,'s','type','float',format=sformat
    printf,lun,'i','height',height,format=iformat
    printf,lun,'i','width',width,format=iformat
    printf,lun,'i','size',4,format=iformat
    printf,lun,'s','machine','SPARC',format=sformat
    printf,lun,'s','format','RDF_range_1',format=sformat
    printf,lun,'.'

    ; Write the image (as binary)

    linefeed = 10B
    writeu,lun,swap_endian(image, /swap_if_little_endian),   $
               swap_endian(linefeed, /swap_if_little_endian)

    ; Write the rdf footer

    for k=0L,nextra-1 do begin
      printf,lun,strtrim(string(extratags[k],format='(a,5x,a16,3x,a,3x,a)'), 2)
    endfor
    printf,lun,'s','target',tname,format=sformat
    printf,lun,'s','polarization',pol,format=sformat
    printf,lun,'.'

    nwrite = nwrite + 1L

  endif
endfor

print,'Wrote ',nwrite,' rdf frames, each containing 1 image, to file ', $
      outfile,format='(a,i0,2a)'

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writetkplay,outfile,chan=chan,help=help,_extra=_ext

; Write the loaded pair to disk as an rdf file readable by tkplay:
; one rdf frame containing the OC spectra and a second rdf frame
; containing the SC spectrum.
;
; Use writerdf to do the work

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 1 or keyword_set(help) then begin
  print,'writetkplay,outfile[,/overwrite][,/append][,chan=1 or 2][,/help]'
  print,"            '.rdf' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writetkplay,outfile[,/overwrite][,/append][,chan=1 or 2][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loaded).ndata le 2 then begin
  print,'ERROR in writetkplay: No pair is loaded, so nothing can be written to disk'
  return
endif

; Call writerdf as needed to write one frame per polarization channel

if n_elements(chan) gt 0 then begin
  writerdf,outfile,chan=chan,/tkplay,_extra=_ext
endif else begin
  writerdf,outfile,chan=1,/tkplay,_extra=_ext
  writerdf,outfile,chan=2,/tkplay,/append
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writetkplay1,outfile,help=help,_extra=_ext

; Write the loaded single-channel spectrum to disk as an rdf file
; readable by tkplay: one rdf frame containing the spectrum
; (either OC or SC).
;
; This is trivial -- that is, identical to what writerdf1 already does.
; I've included it solely for name compatibility with writetkplay,
; writestacktkplay, and writestacktkplay1.

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 1 or keyword_set(help) then begin
  print,'writetkplay1,outfile[,/overwrite][,/append][,/help]'
  print,"             '.rdf' output file extension is added if not already present"
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writetkplay1,outfile[,/overwrite][,/append][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if (*loaded1).ndata le 2 then begin
  print,'ERROR in writetkplay1: No single-channel spectrum is loaded,', $
        ' so nothing can be written to disk',format='(2a)'
  return
endif

; Call writerdf1, which does exactly what we're looking for

writerdf1,outfile,/tkplay,_extra=_ext

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writestacktkplay,outfile,help=help,_extra=_ext

; Take stack pairs in one or more groups and write them to disk as an rdf file
; readable by tkplay: one rdf frame containing all OC spectra and a second rdf
; frame containing all SC spectra.  The keywords g, ming and maxg, or ag specify
; the groups; omitting all of these causes all stack pairs to be written.
;
; Let writestacktkplay1 do most of the work (including input validation)
; rather than copying and pasting mostly identical code

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'writestacktkplay,outfile[,/overwrite][,/append]'
  print,'                 [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'                 [,chan=1 or 2][,/help]'
  print,' '
  print,"    '.rdf' output file extension is added if not already present"
  print,' '
  print,'    Setting keywords g, ming and maxg, or ag writes all pairs in ', $
        'the specified group(s) to an rdf file',format='(2a)'
  print,'    -- these three keyword choices are mutually exclusive'
  print,'    Calling writestacktkplay with none of those keywords set writes ', $
        'all stack pairs to an rdf file',format='(2a)'
  print,' '
  return
endif

if nstack eq 0 then begin
  print,'ERROR in writestacktkplay: The pair stack is empty'
  return
endif

; Store the single-channel stack so that we can use that space

nstack1_store = nstack1
if nstack1_store gt 0 then begin
  storeStack1 = ptrarr(nstack1_store, /allocate_heap)
  for n=0L,nstack1_store-1 do *storeStack1[n] = *stack1[n]
endif
resetstack1

; Send all channels of all stack pairs to the single-channel stack

split,/stack,/silent

; Now call writestacktkplay1 to deal with it

writestacktkplay1,outfile,_extra=_ext

; Replace the original contents of the single-channel stack

resetstack1
nstack1 = nstack1_store
if nstack1 gt 0 then begin
  stack1 = ptrarr(nstack1, /allocate_heap)
  for n=0L,nstack1-1 do *stack1[n] = *storeStack1[n]
  ptr_free,storeStack1
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro writestacktkplay1,outfile,group=g,ming=ming,maxg=maxg,arrayg=ag, $
                      chan=chan,help=help,_extra=_ext

; Take single-channel stack spectra in one or more groups and write them to disk
; as an rdf file readable by tkplay: one rdf frame containing all OC spectra and
; a second rdf frame containing all SC spectra.  The keywords g, ming and maxg,
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
  print,'writestacktkplay1,outfile[,/overwrite][,/append]'
  print,'                  [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'                  [,chan=1 or 2][,/help]'
  print,' '
  print,"    '.rdf' output file extension is added if not already present"
  print,' '
  print,'    Setting keywords g, ming and maxg, or ag writes all single-channel ', $
        'spectra in the specified group(s) to an rdf file',format='(2a)'
  print,'    -- these three keyword choices are mutually exclusive'
  print,'    Calling writestacktkplay1 with none of those keywords set writes ', $
        'all stack1 spectra to an rdf file',format='(2a)'
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'writestacktkplay1,outfile[,/overwrite][,/append]'
  print,'                  [,group=g][,ming=ming,maxg=maxg][,arrayg=ag]'
  print,'                  [,chan=1 or 2][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif

if nstack1 eq 0 then begin
  print,'ERROR in writestacktkplay1: The single-channel stack is empty'
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
    print,'ERROR in writestacktkplay1: Must specify range with ming <= maxg'
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
    print,'writestacktkplay1: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
endif else if n_elements(ag) gt 0 then begin
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroups then begin
    print,"ERROR in writestacktkplay1: Don't include duplicate group numbers ", $
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
    print,'ERROR in writestacktkplay1: Group #',groupvec[j], $
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
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both). Stokes not implemented'
  return
endelse

; A tedious check:
; Go through the single-channel stack one polarization channel at a time,
; checking that all spectra to be written for that channel (rdf frame)
; have the same essential parameters (target name, ndata, ntags)

chanstring = ['OC','SC']
nwrite = lonarr(2)

for ch=1,2 do begin
  if writepol[ch-1] then begin
    for n=0L,nstack1-1 do begin

      ; See if this spectrum is included in one of the groups which will be written

      group_el = where( (groupvec eq (*stack1[n]).group) and $
                        (chanstring[ch-1] eq (*stack1[n]).pol), count)

      if count gt 0 then begin

        ; This spectrum will be written to disk later if all goes well

        tname = (*stack1[n]).tname
        ndata = (*stack1[n]).ndata
        ntags = (*stack1[n]).ntags

        if nwrite[ch-1] eq 0 then begin
          tname_first = strlowcase(tname)
          ndata_first = ndata
          ntags_first = ntags
        endif else if strlowcase(tname) ne tname_first then begin
          print,"ERROR in writestacktkplay1: Can't write ",chanstring[ch-1], $
                       " spectra with different",format='(3a)'
          print,'                            target names to the same rdf frame'
          return
        endif else if ndata ne ndata_first then begin
          print,"ERROR in writestacktkplay1: Can't write ",chanstring[ch-1], $
                       " spectra with different",format='(3a)'
          print,'                            # of spectral points to the same rdf frame'
          return
        endif else if ntags ne ntags_first then begin
          print,"ERROR in writestacktkplay1: Can't write ",chanstring[ch-1], $
                       " spectra with different",format='(3a)'
          print,'                            # of cw tags to the same rdf frame'
          return
        endif

        nwrite[ch-1] = nwrite[ch-1] + 1L

      endif
    endfor
  endif
endfor

; Quit if there are no spectra to be written

if total(nwrite) eq 0 then begin
  print,'writestacktkplay1: No spectra meet these group/polarization criteria'
  return
endif

; Open the output file

err = openOutfile(lun,outfile,'rdf',/get_lun,_extra=_ext)
if err ne 0 then return

; Go through the single-channel stack one polarization channel at a time,
; writing out all spectra in that channel which are included in one of
; the specified groups

sformat = '(a,5x,a16,3x,a)'
iformat = '(a,5x,a16,3x,i0)'

for ch=1,2 do begin
  if writepol[ch-1] then begin
    writtenToFrame = 0L
    for n=0L,nstack1-1 do begin

      ; See if this spectrum is included in one of the groups which will be written

      group_el = where( (groupvec eq (*stack1[n]).group) and $
                        (chanstring[ch-1] eq (*stack1[n]).pol), count)

      if count gt 0 then begin

        ; This spectrum should be written to disk

        writtenToFrame = writtenToFrame + 1L
        spec = (*stack1[n]).spec
        tags1 = (*stack1[n]).tags
        ndata = (*stack1[n]).ndata
        ntags = (*stack1[n]).ntags
        tname = (*stack1[n]).tname

        ; Fix three tag values (xjcen, jsnr1, jsnr2) so that the output is consistent
        ; with tkplay, which counts these tags from 1 rather than from 0
        ;
        ; Also fix the jcp tag, which tkplay counts (sometimes) from 0 rather than from 1

        tags1.xjcen = tags1.xjcen + 1L
        tags1.jsnr1 = tags1.jsnr1 + 1L
        tags1.jsnr2 = tags1.jsnr2 + 1L
        tags1.jcp = tags1.jcp - 1L
        if tags1.jcp lt 0 then begin
          print,'ERROR in writestacktkplay1: jcp tag for OC is already 0 rather than 1'
          return
        endif

        ; If this is the first spectrum written to this frame,
        ; write the frame header

        if writtenToFrame eq 1 then begin
          printf,lun,'s','type','float',format=sformat
          printf,lun,'i','ndata',ndata,format=iformat
          printf,lun,'i','height',nwrite[ch-1],format=iformat
          printf,lun,'i','width',ndata+ntags,format=iformat
          printf,lun,'i','size',4,format=iformat
          printf,lun,'s','machine','SPARC',format=sformat
          printf,lun,'i','nchan',1,format=iformat
          printf,lun,'i','ntags',ntags,format=iformat
          printf,lun,'i','nspec',nwrite[ch-1],format=iformat
          printf,lun,'s','format','RDF_spectrum_1',format=sformat
          printf,lun,'.'
        endif

        ; Write the spectrum and tags (as binary)
        ; -- write all tags as floating point to be compatible with tkplay

        writeu,lun,swap_endian(spec, /swap_if_little_endian)
        for k=0L,ntags-1 do begin
          writeu,lun,swap_endian(float(tags1.(k)), /swap_if_little_endian)
        endfor

        ; If this is the last spectrum written to this frame,
        ; write the rdf footer, including a set of extra tags which is
        ; common to all spectra which were just written
        ;
        ; Include a comment that the format was converted

        if writtenToFrame eq nwrite[ch-1] then begin
          mergeExtraStack1,ch,extratags1,nextra1,groupvec=groupvec
          extratags1 = [extratags1, extratags1[0]]
          extratags1[nextra1].format = ''
          extratags1[nextra1].name = ''
          extratags1[nextra1].value = ''
          extratags1[nextra1].comment = $
                   '# Conv to   tkplay format ' + systime(/utc) + ' UT'
          nextra1 = nextra1 + 1L
          for k=0L,nextra1-1 do begin
            printf,lun,strtrim(string(extratags1[k],format='(a,5x,a16,3x,a,3x,a)'), 2)
          endfor
          printf,lun,'s','target',tname,format=sformat
          printf,lun,'s','polarization',chanstring[ch-1],format=sformat
          printf,lun,'.'
        endif

      endif
    endfor

    if nwrite[ch-1] gt 0 then begin
      print,'Wrote 1 rdf frame, containing ',nwrite[ch-1],' ',chanstring[ch-1], $
            ' spectra, to file ',outfile,format='(a,i0,4a)'
    endif

  endif
endfor

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function rdfHeaderOK,fileformat,machine,dataformat,bytesPerValue, $
                     ndata,ntags,width,height,nspec,nchan,isCW

; Check that rdf header tags have reasonable values;
; print error messages if they don't
;
; This function is called by readrdf

common loadedBlock,loadedi,loaded1,loaded
common channelBlock, chanstrings, maxchan

; Find out how many cw tags are assumed by these reduction procedures

ntagsAssumed = (*loaded1).ntags

; Do the various checks

if fileformat eq 'rdf_spectrum_1' then begin
  isCW = 1
endif else if fileformat eq 'rdf_range_1' or fileformat eq 'rdf_image' then begin
  isCW = 0
  fileformat = 'rdf_range_1'
endif else begin
  print,"ERROR: format = '",fileformat,"'"
  print,'       readrdf is only set up to handle files of types'
  print,'               rdf_spectrum_1 or rdf_range_1 / rdf_image'
  return, 0
endelse

if machine ne 'sparc' then begin
  print,"ERROR: machine = '",machine,"'"
  print,"       readrdf is only set up to read files written on a 'sparc' machine"
  return, 0
endif else if dataformat ne 'float' then begin
  print,"ERROR: type = '",dataformat,"'"
  print,"       readrdf is only set up to handle type = 'float'"
  return, 0
endif else if bytesPerValue ne 4 then begin
  print,'ERROR: size = ',bytesPerValue,format='(a,i0)'
  print,'       readrdf is only set up to handle 4 bytes ', $
        'per data point or tag value',format='(2a)'
  return, 0
endif

if isCW then begin
  if ndata le 0 then begin
    print,'ERROR: ndata = ',ndata,format='(a,i0)'
    print,'       readrdf needs > 0 data points per spectrum'
    return, 0
  endif else if ntags ne ntagsAssumed then begin
    print,'ERROR: ntags = ',ntags,format='(a,i0)'
    print,'       readrdf is only set up for ',ntagsAssumed,' tags per spectrum', $
          format='(a,i0,a)'
    return, 0
  endif else if height le 0 then begin
    print,'ERROR: height = ',height,format='(a,i0)'
    print,'       readrdf needs > 0 spectra per rdf frame'
    return, 0
  endif else if width ne (ndata + ntags) then begin
    print,'ERROR: width = ',width,format='(a,i0)'
    print,'       readrdf expects width = ',(ndata + ntags),format='(a,i0)'
    print,'       (',ndata,' data points + ',ntags,' tags per spectrum)', $
          format='(a,i0,a,i0,a)'
    return, 0
  endif else if nspec ne height then begin
    print,'ERROR: nspec = ',nspec,format='(a,i0)'
    print,'       readrdf expects nspec = height (= ',height,')',format='(a,i0,a)'
    return, 0
  endif else if (nchan lt 1 || nchan gt maxchan) then begin
    print,'ERROR: nchan = ',nchan,format='(a,i0)'
    print,'       readrdf is only set up to handle 1 to ', strtrim(string(maxchan),2), ' channels per rdf frame'
    return, 0
  endif
endif else begin
  if height le 0 then begin
    print,'ERROR: height = ',height,format='(a,i0)'
    print,'       readrdf needs > 0 range bins (rows) per rdf image'
    return, 0
  endif else if width le 0 then begin
    print,'ERROR: width = ',width,format='(a,i0)'
    print,'       readrdf needs > 0 frequency bins (columns) per rdf image'
    return, 0
  endif
endelse

return, 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function polOK,isCW,height,nchan,polarization,jcp_in,pol_in,loose=loose,addpol=addpol

; For images:
;   Check that polarization = 'OC' or 'SC' and set
;   pol_in to a one-element array equal to this value.
;
; For CW spectra:
;   If the polarization and jcp tags aren't in accord, return false;
;   otherwise, set the pol_in elements to the appropriate values
;   ('OC' and/or 'SC').  Procedure used depends on whether there's
;   one channel (nchan = 1) or two being read in.  Simple but tedious....
;
;   For rdf frames with nchan = 2, there must be an even number of spectra
;   which alternate jcp values regularly throughout the frame. For
;   example, if nchan = 2 and nspec = 4 then jcp = (1,2,1,2) or (2,1,2,1)
;   would be OK, but not (1,2,2,1) or (2,2,2,1).
;
; Feb 2002: Add the loose keyword, which (for nchan=1) accepts CW spectra
;           so long as the polarization tag is properly set; this can be
;           necessary for reading tkplay-generated files, which at least
;           sometimes have the jcp tag incorrectly set (to 0 and 1 rather
;           than 1 and 2 for OC and SC, respectively).
;
; Feb 2004: Add the addpol keyword, which allows you to set an image's
;           polarization to 'OC' or 'SC' if that tag is missing in the rdf file
;
; Aug 2010: When /loose is set, accept a CW frame that claims nchan = 2 even
;           though polarization is 'OC' or 'SC': treat it as if nchan = 1
;
; This function is called by readrdf and readfits

common channelBlock, chanstrings, maxchan

if isCW and nchan eq 2 and (polarization eq 'OC' or polarization eq 'SC') $
        and keyword_set(loose) then begin
  nchan = 1
endif

if (not isCW) and (polarization eq 'OC' or polarization eq 'SC') then begin
  jcp_in = [-1]
  pol_in = [polarization]
  return, 1
endif else if (not isCW) and isnull(polarization) and notnull(addpol) then begin
  jcp_in = [-1]
  pol_in = [strupcase(addpol)]
  return, 1
endif else if isCW and polarization eq 'OC' then begin
  if nchan eq 1 then begin
    n_bad = where(jcp_in ne 1, count)
    if keyword_set(loose) or count eq 0L then begin
      jcp_in = replicate(1, height)
      pol_in = replicate('OC', height)
      return, 1
    endif else begin
      print,"ERROR in polOK: polarization = 'OC' but spectrum #",n_bad[0]+1, $
            " has jcp = ",jcp_in[n_bad[0]],format='(a,i0,a,i0)'
      print,'(To ignore the jcp tag, use readrdf,infile,/loose)'
      print,'(To import from tkplay, use readrdf,infile,/tkplay)'
      return, 0
    endelse
  endif else begin
    print,"ERROR in polOK: polarization = 'OC' when there are 2 channels per frame"
    print,'(To treat the frame as single-channel data, use readrdf,infile,/loose)'
    return, 0
  endelse
endif else if isCW and polarization eq 'SC' then begin
  if nchan eq 1 then begin
    n_bad = where(jcp_in ne 2, count)
    if keyword_set(loose) or count eq 0L then begin
      jcp_in = replicate(2, height)
      pol_in = replicate('SC', height)
      return, 1
    endif else begin
      print,"ERROR in polOK: polarization = 'SC' but spectrum #",n_bad[0]+1, $
            " has jcp = ",jcp_in[n_bad[0]],format='(a,i0,a,i0)'
      print,'(To ignore the jcp tag, use readrdf,infile,/loose)'
      print,'(To import from tkplay, use readrdf,infile,/tkplay)'
      return, 0
    endelse
  endif else begin
    print,"ERROR in polOK: polarization = 'SC' when there are 2 channels per frame"
    print,'(To treat the frame as single-channel data, use readrdf,infile,/loose)'
    return, 0
  endelse
endif else if isCW and isnull(polarization) then begin
  firstjcp = jcp_in[0]
  if nchan eq 1 then begin

    ; Just check that all jcp tags are the same (1 or 2)

    n_bad = where(jcp_in ne firstjcp, count)
    if count eq 0 then begin
      if firstjcp eq 1 then begin
        pol_in = replicate('OC', height)
        return, 1
      endif else if firstjcp eq 2 then begin
        pol_in = replicate('SC', height)
        return, 1
      endif else begin
        print,'ERROR in polOK: Illegal jcp tag value of ',firstjcp,format='(a,i0)'
        return, 0
      endelse
    endif else begin
      print,'ERROR in polOK: nchan = 1 but not all jcp tags are equal'
      return, 0
    endelse
  endif else begin

    ; nchan = 2: Insist that there's an even number of spectra with
    ;            jcp tags alternating regularly between 1 and 2

    n = 0L
    alternating = 1
    if nchan eq 2 then begin
      while alternating and (n lt height-1) do begin
        alternating = (jcp_in[n] eq 1 and jcp_in[n+1] eq 2) or $
                      (jcp_in[n] eq 2 and jcp_in[n+1] eq 1)
        n = n + 1L
      endwhile
    endif else begin
      print, "There are ", strtrim(string(nchan),2), "channels: assuming correct alternating pols"
    endelse
    if alternating then begin
      pol_in = strarr(height)
      for n=0L,height-1 do begin
        pol_in[n] = chanstrings[jcp_in[n]-1]
      endfor
      return, 1
    endif else begin
      print,"ERROR in polOK: There are 2 channels per frame but jcp tags don't ", $
            "alternate regularly between 1 and 2",format='(2a)'
      print,'      jcp values are ',jcp_in
      return, 0
    endelse
  endelse
endif else begin
  print,"ERROR in polOK: polarization = '",polarization,"'"
  print,"      readrdf and readfits are only set up to handle values 'OC' and 'SC'"
  print,"      -- try using the 'addpol' keyword"
  return, 0
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function tagOK,tag_float,tagtype

; Check that a tag value (just read from disk in float format) can be
; correctly converted to an integer-type format without overflow
;
; Here are the tag types:
;   0=undef, 1=byte, 2=int, 3=long, 4=float, 5=double, 6=complex, 7=string,
;   8=structure, 9=complex double, 10=pointer, 11=object reference,
;   12=unsigned int, 13=unsigned long, 14=64-bit int, 15=unsigned 64-bit int

case tagtype of
  0:    return, 0
  1:    return, (tag_float gt -1 and tag_float lt 256)
  2:    return, (tag_float gt -32769L and tag_float lt 32768L)
  3:    return, (tag_float gt -2147483649LL and tag_float lt 2147483648LL)
  4:    return, 1
  5:    return, 1
  6:    return, 0
  7:    return, 1
  8:    return, 0
  9:    return, 0
  10:   return, 0
  11:   return, 0
  12:   return, (tag_float gt -1L and tag_float lt 65536L)
  13:   return, (tag_float gt -1LL and tag_float lt 4294967296LL)
  14:   return, (tag_float gt -9223372036854775809D0 and tag_float lt 9223372036854775808ULL)
  15:   return, (tag_float gt -1LL and tag_float lt 18446744073709551616D0)
  else: return, 0
endcase

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro readrdf,infile,nopairs=nopairs,loose=loose,addpol=addpol,tkplay=tkplay,silent=silent, $
            help=help

; Read in the entire contents of an rdf-format disk file.
;
; Once the file has been read, the procedure checks whether the file consists
; entirely of two-channel pairs.  It ignores the rdf frame groupings (nspec, nchan)
; in making this check.  Instead it allows two possible polarization sequences:
; a file broken into two halves (all of the spectra from one channel followed by all
; of the spectra, ordered in the same way, from the other); or alternating channels.
; The two halves of each putative pair must have different polarizations but the
; same values of certain other tags (such as target name, number of spectral
; points, and zero-frequency bin).  Note that these other tags must match within a
; pair but not necessarily between different pairs.
;
; Note also that in the alternating case the channels do NOT have to alternate
; regularly throughout the file: The sequence (OC,SC,SC,OC) could count as two pairs
; if all the right tags matched.  (However, when reading an rdf frame with nchan = 2,
; function polOK does require polarizations to alternate regularly within that frame.)
;
; If the file contains nothing but two-channel pairs, the spectra are grouped
; into pairs and are added to the pair stack (stack).  Otherwise no pairs are
; created: All individual spectra are added to the single-channel stack (stack1).
;
; Feb 2002: Add the loose keyword, which allows tkplay-generated nchan=1 rdf files
;           to be read even though the jcp tag has been incorrectly set (i.e., just
;           use the polarization tag from the rdf header or footer to decide which
;           polarization channel we're dealing with)
;
; Aug 2002: If readrdf doesn't find the input file as specified,
;           add an '.rdf' extension (if it's not already there) and try again
;
; Feb 2004: Add the addpol keyword, which allows you to set an image's
;           polarization to 'OC' or 'SC' if that tag is missing in the rdf file
;
; Jun 2005: Add the tkplay keyword, which subtracts 1 from tags xjcen, jsnr1, and
;           jsnr2 for CW spectra imported from tkplay: that package counts these
;           tags from 1, whereas this package counts them from 0
;
; Jun 2006: For the tkplay keyword, also add 1 to the jcp tag (so it's 1 for OC
;           and 2 for SC rather than 0 and 1, as tkplay *sometimes* sets them)
;
; Aug 2010: When /loose is set, (a) treat an allegedly 2-channel CW frame as
;           single-channel data if polarization is set to 'OC' or 'SC', and
;           (b) reset impossibly huge values of the CW 'igw' tag to zero
;
; Feb 2022  Modified to read in multi-channel data. To avoid pain, will only
;           be implemented for 'alternating' rdf files
;

common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if not keyword_set(nopairs) then nopairs = 0
if not keyword_set(loose) then loose = 0
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

err = openInfile(lun,infile,'rdf',/get_lun)
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

; Get the tag names in case there's a problem reading in cw tags

tagnames = strlowcase(tag_names(blanktags1()))

; Read an rdf frame each time through the loop

headline = ''   ; initialize as a string or else readf will assume float
tailline = ''   ; same here
userchoice = '' ; ditto (for read)
nframes = 0L
nreadCW = 0L
nreadi = 0L

while not eof(lun) do begin

  ; Initialize some parameters

  dataformat = ''
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
  fileformat = ''
  extratag = {format:'', name:'', value:'', comment:''}
  nextra = 0L
  has_spb_tag = 0
  has_infile_tag = 0
  infile_fullpath = (file_search(infile,/fully_qualify_path))[0]

  ; Start reading in the rdf header

  if not eof(lun) then begin
    readf,lun,headline,format='(a)'
  endif else begin
    print,' '
    print,'ERROR in readrdf: rdf header is truncated'
    print,' '
    close,lun
    free_lun,lun
    return
  endelse

  ; Check that this isn't a FITS file

  fits_keyname = strtrim(strcompress(strmid(headline,0,8)), 2)
  if fits_keyname eq 'SIMPLE' then begin
    print,' '
    print,'ERROR in readrdf: This is a FITS file, use readfits to read it'
    print,' '
    close,lun
    free_lun,lun
    return
  endif

  while headline ne '.' do begin

    ; Parse this header line into format, name, value, comment;
    ; allow value to include spaces (e.g., '1 Ceres')

    commentstart = strpos(headline,'#')
    if commentstart eq -1 then begin
      comment = ''
      headline = strtrim(strcompress(headline), 2)
    endif else begin
      comment = '# ' + strtrim(strmid(headline,commentstart+1), 2)
      headline = strtrim(strcompress(strmid(headline,0,commentstart)), 2)
    endelse
    headpieces = stregex(headline,'^([^ ]+) ([^ ]+) (.+)$', /subexpr, /extract)
    headformat = headpieces[1]
    name = strlowcase(headpieces[2])
    namepieces = stregex(name,'^(.+)(_[o|s]c)$', /subexpr, /extract)
    if notnull(namepieces[2]) then name = namepieces[1] + strupcase(namepieces[2])
    if strlen(name) gt 16 then $
      print,'WARNING in readrdf: Extra tag name ',name,' will be truncated ', $
            'when displayed or written to disk',format='(3a)'
    value = strtrim(headpieces[3], 2)

    case name of
      'type'        : dataformat = strlowcase(value)
      'ndata'       : ndata = long(value)
      'height'      : height = long(value)
      'width'       : width = long(value)
      'size'        : bytesPerValue = long(value)
      'polarization': polarization = strupcase(value)
      'target'      : targetname = value
      'object'      : targetname = value
      'machine'     : machine = strlowcase(value)
      'nchan'       : nchan = long(value)
      'ntags'       : ntags = long(value)
      'nspec'       : nspec = long(value)
      'format'      : fileformat = strlowcase(value)
      else          : begin
                        if notnull(name) or notnull(comment) then begin
                          if name eq 'infile' then begin
                            value = infile_fullpath
                            has_infile_tag = 1
                          endif else if name eq 'smpb' and not has_spb_tag then begin
                            name = 'samples_per_baud'
                            has_spb_tag = 1
                          endif
                          extratag.format = headformat
                          extratag.name = name
                          extratag.value = value
                          extratag.comment = comment
                          extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
                          nextra = nextra + 1L
                        endif
                      endcase
    endcase

    ; Read the next header line

    if not eof(lun) then begin
      readf,lun,headline,format='(a)'
    endif else begin
      print,' '
      print,'ERROR in readrdf: rdf header is truncated'
      print,' '
      close,lun
      free_lun,lun
      return
    endelse

  endwhile

  ; Check that the header is OK

  if not rdfHeaderOK(fileformat,machine,dataformat,bytesPerValue, $
                     ndata,ntags,width,height,nspec,nchan,isCW) then begin
    close,lun
    free_lun,lun
    return
  endif

  ; Read in the data and (for CW spectra) tags

  if isCW then begin
    spec1 = fltarr(ndata)
    spec_in = fltarr(height, ndata)
    fltTags = fltarr(ntags)
    tags_in = replicate(blanktags1(), height)
    for n=0L,height-1 do begin
      readu,lun,spec1,fltTags
      spec1 = swap_endian(temporary(spec1), /swap_if_little_endian)
      fltTags = swap_endian(temporary(fltTags), /swap_if_little_endian)
      spec_in[n,*] = spec1
      for tagnum=0L,ntags-1 do begin

        ; Make sure there will be no integer overflow when converting cw tags

        if not tagOK(fltTags[tagnum],tagtypes[tagnum]) then begin
          if tagnames[tagnum] eq 'igw' and keyword_set(loose) then begin
            fltTags[tagnum] = 0
          endif else begin
            print,' '
            print,'WARNING: cw tag value too large for integer conversion'
            print,'         Spectrum #',n+1,' in rdf frame #',nframes+1L,' has ', $
                       tagnames[tagnum],' = ',fltTags[tagnum],format='(a,i0,a,i0,3a,f0)'
            print,' '
            userprompt = '         Reset to zero and continue? (y/n) '
            read,userchoice,prompt=userprompt
            userchoice = strlowcase(strmid(strtrim(userchoice), 0, 1))
            if userchoice eq 'y' or isnull(userchoice) then begin
              fltTags[tagnum] = 0
              print,'resetting tag value to zero'
            endif else begin
              close,lun
              free_lun,lun
              print,'readrdf canceled'
              return
            endelse
          endelse
        endif
        tags_in[n].(tagnum) = fix(fltTags[tagnum], type=tagtypes[tagnum], /print)
      endfor
    endfor
  endif else begin
    image_in = fltarr(width, height)
    linefeed = 10B
    readu,lun,image_in,linefeed
    image_in = swap_endian(temporary(image_in), /swap_if_little_endian)
  endelse

  ; Read in the rdf tail (footer)

  if not eof(lun) then begin
    readf,lun,tailline,format='(a)'
  endif else begin
    print,' '
    print,'ERROR in readrdf: rdf footer is truncated'
    print,' '
    close,lun
    free_lun,lun
    return
  endelse
  while tailline ne '.' do begin
    commentstart = strpos(tailline,'#')
    if commentstart eq -1 then begin
      comment = ''
      tailline = strtrim(strcompress(tailline), 2)
    endif else begin
      comment = '# ' + strtrim(strmid(tailline,commentstart+1), 2)
      tailline = strtrim(strcompress(strmid(tailline,0,commentstart)), 2)
    endelse
    tailpieces = stregex(tailline,'^([^ ]+) ([^ ]+) (.+)$', /subexpr, /extract)
    tailformat = tailpieces[1]
    name = strlowcase(tailpieces[2])
    namepieces = stregex(name,'^(.+)(_[o|s]c)$', /subexpr, /extract)
    if notnull(namepieces[2]) then name = namepieces[1] + strupcase(namepieces[2])
    if strlen(name) gt 16 then $
      print,'WARNING in readrdf: Extra tag name ',name,' will be truncated ', $
            'when displayed or written to disk',format='(3a)'
    value = tailpieces[3]
    case name of
      'target'      : targetname = value
      'object'      : targetname = value
      'polarization': polarization = strupcase(value)
      else          : begin
                        if notnull(name) or notnull(comment) then begin
                          if name eq 'infile' then begin
                            value = infile_fullpath
                            has_infile_tag = 1
                          endif else if name eq 'smpb' and not has_spb_tag then begin
                            name = 'samples_per_baud'
                            has_spb_tag = 1
                          endif
                          extratag.format = tailformat
                          extratag.name = name
                          extratag.value = value
                          extratag.comment = comment
                          extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
                          nextra = nextra + 1L
                        endif
                      endcase
    endcase
    if not eof(lun) then begin
      readf,lun,tailline,format='(a)'
    endif else begin
      print,' '
      print,'ERROR in readrdf: rdf footer is truncated'
      print,' '
      close,lun
      free_lun,lun
      return
    endelse
  endwhile

  ; If requested, fix the jcp tag for spectra imported from tkplay, which
  ; *sometimes* counts this tag from 0 rather than from 1

  if isCW and keyword_set(tkplay) then begin
    tags_in[*].jcp = tags_in[*].jcp + 1L
    if total(tags_in[*].jcp gt 2) gt 0 then begin
      print,'ERROR in readrdf: jcp tags are already [1,2] rather than [0,1]'
      return
    endif
  endif

  ; Check that the polarization tag(s) are OK;
  ; if so, set the pol_in array elements to 'OC' and/or 'SC'.

  jcp_in = (isCW) ? tags_in[*].jcp : [-1]
  if not polOK(isCW,height,nchan,polarization,jcp_in,pol_in,  $
               loose=loose,addpol=addpol) then begin
    close,lun
    free_lun,lun
    return
  endif
  if isCW and keyword_set(loose) then tags_in[*].jcp = jcp_in

  ; Add the 'infile' extra tag if necessary

  if not has_infile_tag then begin
    extratag.format = 's'
    extratag.name = 'infile'
    extratag.value = infile_fullpath
    extratag.comment = ''
    extratags = (nextra eq 0) ? [extratag] : [extratags, extratag]
    nextra = nextra + 1L
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

  ; If requested, fix three tag values (xjcen, jsnr1, jsnr2) for spectra imported
  ; from tkplay, which counts these tags from 1 rather than from 0

  if isCW and keyword_set(tkplay) then begin
    tags_in[*].xjcen = tags_in[*].xjcen - 1L
    tags_in[*].jsnr1 = tags_in[*].jsnr1 - 1L
    tags_in[*].jsnr2 = tags_in[*].jsnr2 - 1L
    extratags = [extratags, extratags[0]]
    extratags[nextra].format = ''
    extratags[nextra].name = ''
    extratags[nextra].value = ''
    extratags[nextra].comment = $
           '# Conv from tkplay format ' + systime(/utc) + ' UT'
    nextra = nextra + 1L
  endif

  ; Turn each spectrum into a structure and add it to a temporary stack

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

  ; Back to the top of the main loop: Read in another rdf frame

endwhile

; All frames have been read in: close the input file

close,lun
free_lun,lun

; Check whether all CW spectra can be paired

if nreadCW eq 0 or nopairs or (nchan*(nreadCW/nchan) ne nreadCW) then begin

  ; Treat as single spectra, not pairs:
  ; no spectra, or an odd number of spectra, or else the "don't pair" flag was set

  pairSpectra = 0

endif else begin

  ; See if the file is broken into halves: first all the spectra from one
  ; polarization channel, then all the corresponding spectra from another

  nhalf = nreadCW/nchan
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
  ; corresponding spectrum from another channel

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
    n = n + nchan
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
      nvec = indgen(nchan) * nhalf + n
      if (*stackRead1[n]).pol eq 'OC' then begin
        n1 = n
        n2 = n + nhalf
      endif else begin
        n1 = n + nhalf
        n2 = n
      endelse
    endif else begin
      nvec = indgen(nchan) + n * nchan
      if (*stackRead1[nchan*n]).pol eq 'OC' then begin
        n1 = nchan*n
        n2 = nchan*n + 1
      endif else begin
        n1 = nchan*n + 1
        n2 = nchan*n
      endelse
    endelse
  ; OC and SC are special 
    nvec[0] = n1
    nvec[1] = n2
    pair_ntags =  (*stackRead1[n1]).ntags
    pair_tags = blanktags(npol=nchan)
  ; Need to loop in case they are out of order
    for j = 0, nchan-1 do begin
      pair_tags[j] = (*stackRead1[nvec[j]]).tags
    endfor
;    for k=0L,pair_ntags-1 do begin
;      pair_tags.(k) = (*stackRead1).tags.(k)
;    endfor
 ; For the extratags, just use OC and SC. 
    mergeExtra,(*stackRead1[n1]).extratags,(*stackRead1[n2]).extratags, $
               (*stackRead1[n1]).nextra,(*stackRead1[n2]).nextra, $
               pair_extratags,pair_nextra
    stackStruc = {group:groupCW, $
                  spec:fltarr(nchan,(*stackRead1[n1]).ndata), $
                  tags:pair_tags, $
                  extratags:pair_extratags, $
                  ndata:(*stackRead1[n1]).ndata, $
                  ntags:pair_ntags, $
                  nextra:pair_nextra, $
                  tname:(*stackRead1[n1]).tname}
 ; need to loop in case they are out of order
    for j = 0, nchan-1 do begin
      stackStruc.spec[j,*] = (*stackRead1[nvec[j]]).spec
    endfor
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
pro cw2csv,outfile,stack=n,chan=chan,help=help,_extra=_ext

; Write one channel of spectral data to a dat (ASCII text) file:
; -- Do this for the loaded pair, or else for stack pair n
;         if the stack keyword is set
; -- Do this for the OC spectrum, unless the /sc flag is set

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

if n_params() ne 1 or keyword_set(help) then begin
  print,'cw2csv,outfile[,/overwrite][,/append][,stack=n][,chan=chan][,/help]'
  print,' '
  print,'Write spectral data to a csv (ASCII text) file'
  print,' '
  print,'   Do this for the loaded pair, or else for stack pair n if the stack keyword is used'
  print,' '
  print,'   Do this for all spectrum, unless chan is set. Chan can be a list like [1,2,4]'
  print,' '
  return
endif else if size(outfile, /type) ne 7 then begin
  print,' '
  print,'cw2csv,outfile[,/overwrite][,/append][,stack=n][,/sc][,/help]'
  print,'Make sure that outfile is a quoted string!'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in cw2csv: The pair stack is empty, so nothing can be written to disk'
    return
  endif else if n le 0 then begin
    print,'ERROR in cw2csv: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in cw2csv: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in cw2csv: No pair is loaded, so nothing can be written to disk'
  return
endif

; Get the spectrum to be written to disk

if n_elements(n) gt 0 then begin
  pair = (*stack[n-1]).spec
  freq = (*stack[n-1]).freq
  npol = n_elements((*stack[n-1]).tags)
endif else begin
  pair = (*loaded).spec
  freq = (*loaded).freq
  npol = n_elements((*loaded).tags)
endelse
writepols = intarr(npol)
if n_elements(chan) eq 0 then writepols = writepols+1 else if n_elements(chan) eq 1 then begin
  if chan lt 0 || chan gt npol then begin
    print, "ERROR in cw2csv. chans must be 1 .. npol"
    return
  endif
  writepols[chan-1] = 1
endif else begin
  for i = 1, n_elements(chan) do begin
    if chan[i-1] lt 0 || chan[i-1] gt npol then begin
      print, "ERROR in cw2csv. chans must be 1 .. npol"
      return
    endif
    writepols[chan[i-1]-1] = 1
  endfor
endelse

; Write the data and close up shop

nwrite = total(writepols)
a = transpose(freq)
header=strarr(nwrite+1)
header[0] = "#Frequency"
for i = 0, npol-1 do begin
  if writepols[i] then begin
    a = [a, pair[i,*]]
    header[i+1] = chanstrings[i]
  endif
endfor
print, size(a)
write_csv, outfile, a, header=header

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cw2dat,outfile,sc=sc,stack=n,chan=chan,help=help,_extra=_ext

; Write one channel of spectral data to a dat (ASCII text) file:
; -- Do this for the loaded pair, or else for stack pair n
;         if the stack keyword is set
; -- Do this for the OC spectrum, unless the /sc flag is set

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

if n_params() ne 1 or keyword_set(help) then begin
  print,'cw2dat,outfile[,/overwrite][,/append][,stack=n][,/sc][,chan=chan][,/help]'
  print,' '
  print,'Write one channel of spectral data to a dat (ASCII text) file'
  print,' '
  print,'   Do this for the loaded pair, or else for stack pair n if the stack keyword is used'
  print,' '
  print,'   Do this for the OC spectrum, unless chan or /sc is set'
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
  npol = n_elements((*stack[n-1]).tags)
endif else begin
  pair = (*loaded).spec
  npol = n_elements((*loaded).tags)
endelse
if n_elements(chan) eq 0 then chan = 1
if chan lt 0 || chan gt npol then begin
  print, "ERROR in cw2dat. chan must be 1 .. npol"
  return
endif
if keyword_set(sc) then chan = 2
spec = reform(pair[chan,*])
chanstring = chanstrings[chan-1]

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
