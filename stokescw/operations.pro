;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function cubicconvpoly,x,cubic=cubic

; Evaluate the polynomial used to generate interpolation coefficients
; for cubic convolution

if n_elements(cubic) eq 0 then cubic = -1.0

absx = abs(x)
if absx lt 1.0 then begin
  interp_poly = (cubic + 2)*absx*absx*absx - (cubic + 3)*absx*absx + 1.0
endif else if absx lt 2.0 then begin
  interp_poly = cubic*(absx*absx*absx - 5*absx*absx + 8*absx - 4.0)
endif else begin
  interp_poly = 0.0
endelse

return,interp_poly

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function noisecorr,codemethod,spb,stride

; Compute a vector of noise correlation coefficients for image pixels that
; are in the same Doppler column and are a given number of delay rows apart;
; for example, corr[1] is the coefficient for pixels that are one row apart.
; The vector is just long enough to hold all nonzero coefficients.
;
; The expressions used here are asymptotically valid as fft length
; and (for short-code images) code length approach infinity.  Slight Doppler
; variations (for spb > 1) have been ignored, on the assumption that Doppler
; frequencies of interest are small compared to the unaliased bandwidth.

; First, compute a vector of noise correlation coefficients for pixels
; that are j samples apart (= j/stride rows apart)

n_corr = (codemethod eq 'long_orig') ? spb : 2*spb - 1
j = indgen(n_corr)
if codemethod eq 'long_orig' then begin
  corr = ( 1.0 - j/(1.0*spb) )^2
endif else begin
  indexmask = intarr(n_corr)
  indexmask[0:spb-1] = 1
  corr = ( ( (2*spb-j-1)*(2*spb-j)*(2*spb-j+1.0)            $
              - indexmask*4*(spb-j-1)*(spb-j)*(spb-j+1) )   $
           / ((2*spb)*(2*spb^2 + 1)) )^2
endelse

; Now collapse the vector so that, for example, corr[2] is the noise
; correlation coefficient for image pixels in the same Doppler column
; that are 2 rows apart (rather than 2 samples apart)

corr = corr[where(j mod stride eq 0)]

return, corr

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdevfactor_bilinear,rowshift,colshift,codemethod,spb,stride

; Compute the factor by which the r.m.s. noise is reduced when
; bilinear interpolation is used on an image

t = rowshift - floor(rowshift)    ; 0.0 <= t < 1.0
u = colshift - floor(colshift)    ; 0.0 <= u < 1.0
xpoly = [1.0 - t, t]
ypoly = [1.0 - u, u]
coeff = xpoly # ypoly

corr = noisecorr(codemethod,spb,stride)
if n_elements(corr) eq 1 then corr = [corr, 0.0]

sum = 0.0
for i=0L,1 do begin
  for j=0L,1 do begin
    corrcoeff = corr[abs(i-j)]
    for k=0L,1 do begin
      sum = sum + coeff[i,k]*coeff[j,k]*corrcoeff
    endfor
  endfor
endfor
sdevfactor = sqrt(sum)

return,sdevfactor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdevfactor_linear_cw,binshift

; Compute the factor by which the r.m.s. noise is reduced when
; linear interpolation is used on a spectrum

u = binshift - floor(binshift)    ; 0.0 <= u < 1.0
coeff = [1.0 - u, u]

sdevfactor = sqrt(total(coeff^2))

return,sdevfactor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdevfactor_cubicconv,rowshift,colshift,codemethod,spb,stride, $
                              cubic=cubic

; Compute the factor by which the r.m.s. noise is reduced when
; an image is interpolated using cubic convolution

if n_elements(cubic) eq 0 then cubic = -1.0

t = rowshift - floor(rowshift)    ; 0.0 <= t < 1.0
u = colshift - floor(colshift)    ; 0.0 <= u < 1.0
xpoly = fltarr(4)
ypoly = fltarr(4)
for k=0L,3 do begin
  xpoly[k] = cubicconvpoly(k-1-t, cubic=cubic)
  ypoly[k] = cubicconvpoly(k-1-u, cubic=cubic)
endfor
coeff = xpoly # ypoly

corr = noisecorr(codemethod,spb,stride)
ncorr = n_elements(corr)
while ncorr lt 4 do begin
  corr = [corr, 0.0]
  ncorr = ncorr + 1
endwhile

sum = 0.0
for i=0L,3 do begin
  for j=0L,3 do begin
    corrcoeff = corr[abs(i-j)]
    for k=0L,3 do begin
      sum = sum + coeff[i,k]*coeff[j,k]*corrcoeff
    endfor
  endfor
endfor
sdevfactor = sqrt(sum)

return,sdevfactor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdevfactor_cubicconv_cw,binshift,cubic=cubic

; Compute the factor by which the r.m.s. noise is reduced when
; a spectrum is interpolated using cubic convolution

if n_elements(cubic) eq 0 then cubic = -1.0

u = binshift - floor(binshift)    ; 0.0 <= u < 1.0
coeff = fltarr(4)
for k=0L,3 do coeff[k] = cubicconvpoly(k-1-u, cubic=cubic)

sdevfactor = sqrt(total(coeff^2))

return,sdevfactor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function gettag,tagname,chan=chan,stack=n,help=help

; Get a tag value for the loaded pair, or for stack pair #n
; if the stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'value = gettag(tagname[,chan=1 or 2][,stack=n][,/help])'
  print,'        tagname must be a quoted string'
  print,' '
  print,'Returns a cw tag value for the loaded pair, or else for stack pair #n'
  print,'if the stack keyword is set'
  print,' '
  print,'If chan=1 or chan=2, a scalar is returned;'
  print,'if chan is omitted, an array with all of the polarizations is returned'
  print,' '
  return, ''
endif else if size(tagname, /type) ne 7 then begin
  print,' '
  print,'value = gettag(tagname[,chan=1 or 2][,stack=n][,/help])'
  print,'Make sure that tagname is a quoted string!'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in gettag: The pair stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in gettag: Must have n >= 1'
    return, ''
  endif else if n gt nstack then begin
    print,'ERROR in gettag: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in gettag: No pair is loaded'
  return, ''
endif

; Get some elements of the loaded pair or of stack pair #n

tags = (n_elements(n) gt 0) ? (*stack[n-1]).tags : (*loaded).tags
tagnames = strlowcase(tag_names(tags[0]))

; Check that tagname is a valid tag name

tagname = strlowcase(tagname)
tagnum = where(tagnames eq tagname, count)
if (count eq 0) then begin
  print,tagname,' is not a valid tag name'
  return, ''
endif else begin
  tagnum = tagnum[0]
endelse

; Return the tag value(s) for the specified channel(s)

if (n_elements(chan) eq 0) then begin
  return, tags.(tagnum)
endif else begin
  if chan lt 1 or chan gt n_elements(tags) then begin
    print, 'Chan must be >= 1 and <= # of pols: ', n_elements(tags)
    return, ''
  endif
  return, tags[chan-1].(tagnum)
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function gettag1,tagname,stack=n,help=help

; Get a tag value for the loaded single-channel spectrum, or for
; stack1 spectrum #n if the stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'value = gettag1(tagname[,stack=n][,/help])'
  print,'        tagname must be a quoted string'
  print,' '
  print,'Returns a cw tag value for the loaded single-channel spectrum,'
  print,'or else for stack1 spectrum #n if the stack keyword is set'
  print,' '
  return, ''
endif else if size(tagname, /type) ne 7 then begin
  print,' '
  print,'value = gettag1(tagname[,stack=n][,/help])'
  print,'Make sure that tagname is a quoted string!'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in gettag1: The single-channel stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in gettag1: Must have n >= 1'
    return, ''
  endif else if n gt nstack1 then begin
    print,'ERROR in gettag1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in gettag1: No single-channel spectrum is loaded'
  return, ''
endif

; Get some elements of the loaded single-channel spectrum
; or of stack1 spectrum #n

tags1 = (n_elements(n) gt 0) ? (*stack1[n-1]).tags : (*loaded1).tags
tagnames = strlowcase(tag_names(tags1))

; Check that tagname is a valid tag name

tagname = strlowcase(tagname)
tagnum = where(tagnames eq tagname, count)
if (count eq 0) then begin
  print,tagname,' is not a valid tag name'
  return, ''
endif else begin
  tagnum = tagnum[0]
endelse

; Return the tag value

return, tags1.(tagnum)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getextra,name,stack=n,help=help

; Returns the value of an extra tag on the loaded pair, or else on
; stack pair n if the stack keyword is set.
;
; /comment returns the comment rather than the tag value
;
; Returns a null string if the tag isn't found.
;
; 2006 Jul 7: Add 'comment' keyword

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if size(name, /type) ne 7 then name = ''  ;  just so it's defined as a string

if keyword_set(help) or isnull(name) or n_params() ne 1 then begin
  print,' '
  print,'value = getextra(name[,stack=n][,/comment][,/help])'
  print,' '
  print,'Return the value of an extra tag for the loaded pair;'
  print,'    if the tag is not present, the null string is returned'
  print,' '
  print,'name must be a quoted string'
  print,' '
  print,'/comment returns the tag comment rather than the tag value'
  print,' '
  print,'stack=n returns a value for stack pair n rather than for the loaded pair'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in getextra: The pair stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getextra: Must have n >= 1'
    return, ''
  endif else if n gt nstack then begin
    print,'ERROR in getextra: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in getextra: No pair is loaded'
  return, ''
endif

; Get the array of extra tags

extratags = (n_elements(n) gt 0) ? (*stack[n-1]).extratags : (*loaded).extratags

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
              print,"ERROR: getextra isn't set up to handle extra tag format '", $
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getextra1,name,stack=n,help=help

; Returns the value of an extra tag on the loaded single-channel
; spectrum, or else on single-channel stack spectrum n if the 
; stack keyword is set.
;
; /comment returns the comment rather than the tag value
;
; Returns a null string if the tag isn't found.
;
; 2006 Jul 7: Add 'comment' keyword

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if size(name, /type) ne 7 then name = ''  ;  just so it's defined as a string

if keyword_set(help) or isnull(name) or n_params() ne 1 then begin
  print,' '
  print,'value = getextra1(name[,stack=n][,/comment][,/help])'
  print,' '
  print,'Return the value of an extra tag for the loaded single-channel spectrum;'
  print,'    if the tag is not present, the null string is returned'
  print,' '
  print,'name must be a quoted string'
  print,' '
  print,'/comment returns the tag comment rather than the tag value'
  print,' '
  print,'stack=n returns a value for stack1 spectrum n rather'
  print,'     than for the loaded single-channel spectrum'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in getextra1: The single-channel stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getextra1: Must have n >= 1'
    return, ''
  endif else if n gt nstack1 then begin
    print,'ERROR in getextra1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in getextra1: No single-channel spectrum is loaded'
  return, ''
endif

; Get the array of extra tags

extratags = (n_elements(n) gt 0) ? (*stack1[n-1]).extratags : (*loaded1).extratags

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
              print,"ERROR: getextra1 isn't set up to handle extra tag format '", $
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getextrai,name,stack=n,comment=comment,help=help

; Returns the value of an extra tag on the loaded image,
; or else on stacki image n if the stack keyword is set.
;
; /comment returns the comment rather than the tag value
;
; Returns a null string if the tag isn't found.
;
; 2006 Jul 7: Add 'comment' keyword

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if size(name, /type) ne 7 then name = ''  ;  just so it's defined as a string

if keyword_set(help) or isnull(name) or n_params() ne 1 then begin
  print,' '
  print,'value = getextrai(name[,stack=n][,/comment][,/help])'
  print,' '
  print,'Return the value of an extra tag for the loaded image;'
  print,'    if the tag is not present, the null string is returned'
  print,' '
  print,'name must be a quoted string'
  print,' '
  print,'/comment returns the tag comment rather than the tag value'
  print,' '
  print,'stack=n returns a value for stacki image n rather than for the loaded image'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in getextrai: The image stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getextrai: Must have n >= 1'
    return, ''
  endif else if n gt nstacki then begin
    print,'ERROR in getextrai: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in getextrai: No image is loaded'
  return, ''
endif

; Get the array of extra tags

extratags = (n_elements(n) gt 0) ? (*stacki[n-1]).extratags : (*loadedi).extratags

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
              print,"ERROR: getextrai isn't set up to handle extra tag format '", $
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getndata,stack=n,help=help

; Get the spectrum length for the loaded pair, or for stack pair #n
; if the stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'value = getndata([stack=n][,/help])'
  print,' '
  print,'Returns the spectrum length for the loaded pair, or else for'
  print,'stack pair #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in getndata: The pair stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getndata: Must have n >= 1'
    return, ''
  endif else if n gt nstack then begin
    print,'ERROR in getndata: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in getndata: No pair is loaded'
  return, ''
endif

; Return the spectrum length

ndata = (n_elements(n) gt 0) ? (*stack[n-1]).ndata : (*loaded).ndata
return, ndata

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getndata1,stack=n,help=help

; Get the spectrum length for the loaded single-channel spectrum,
; or for stack1 spectrum #n if the stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'value = getndata1([stack=n][,/help])'
  print,' '
  print,'Returns the spectrum length for the loaded single-channel spectrum,'
  print,'or else for stack pair #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in getndata1: The single-channel stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getndata1: Must have n >= 1'
    return, ''
  endif else if n gt nstack1 then begin
    print,'ERROR in getndata1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in getndata1: No single-channel spectrum is loaded'
  return, ''
endif

; Return the spectrum length

ndata = (n_elements(n) gt 0) ? (*stack1[n-1]).ndata : (*loaded1).ndata
return, ndata

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getwidth,stack=n,help=help

; Returns the value of the width of the loaded image,
; or else of stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() ne 0 then begin
  print,' '
  print,'value = getwidth([stack=n][,/help])'
  print,' '
  print,'Returns the width of the loaded image, or else of'
  print,'stacki image #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in getwidth: The image stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getwidth: Must have n >= 1'
    return, ''
  endif else if n gt nstacki then begin
    print,'ERROR in getwidth: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in getwidth: No image is loaded'
  return, ''
endif

; Return the image width

width = (n_elements(n) gt 0) ? (*stacki[n-1]).width : (*loadedi).width
return, width

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getheight,stack=n,help=help

; Returns the value of the height of the loaded image,
; or else of stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() ne 0 then begin
  print,' '
  print,'value = getheight([stack=n][,/help])'
  print,' '
  print,'Returns the height of the loaded image, or else of'
  print,'stacki image #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in getheight: The image stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getheight: Must have n >= 1'
    return, ''
  endif else if n gt nstacki then begin
    print,'ERROR in getheight: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loadedi).height le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in getheight: No image is loaded'
  return, ''
endif

; Return the image height

height = (n_elements(n) gt 0) ? (*stacki[n-1]).height : (*loadedi).height
return, height

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function gettarget,stack=n,help=help

; Get the target name for the loaded pair, or for stack pair #n
; if the stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'value = gettarget([stack=n][,/help])'
  print,' '
  print,'Returns the target name for the loaded pair, or else for'
  print,'stack pair #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in gettarget: The pair stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in gettarget: Must have n >= 1'
    return, ''
  endif else if n gt nstack then begin
    print,'ERROR in gettarget: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in getndata: No pair is loaded'
  return, ''
endif

; Return the target name

tname = (n_elements(n) gt 0) ? (*stack[n-1]).tname : (*loaded).tname
return, tname

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function gettarget1,stack=n,help=help

; Get the target name for the loaded single-channel spectrum,
; or for stack1 spectrum #n if the stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'value = gettarget1([stack=n][,/help])'
  print,' '
  print,'Returns the target name for the loaded single-channel spectrum,'
  print,'or else for stack pair #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in gettarget1: The single-channel stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in gettarget1: Must have n >= 1'
    return, ''
  endif else if n gt nstack1 then begin
    print,'ERROR in gettarget1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in gettarget1: No single-channel spectrum is loaded'
  return, ''
endif

; Return the target name

tname = (n_elements(n) gt 0) ? (*stack1[n-1]).tname : (*loaded1).tname
return, tname

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function gettargeti,stack=n,help=help

; Returns the target name for the loaded image,
; or else for stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() ne 0 then begin
  print,' '
  print,'value = gettargeti([stack=n][,/help])'
  print,' '
  print,'Returns the target name for the loaded image, or else for'
  print,'stacki image #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in gettargeti: The image stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in gettargeti: Must have n >= 1'
    return, ''
  endif else if n gt nstacki then begin
    print,'ERROR in gettargeti: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in gettargeti: No image is loaded'
  return, ''
endif

; Return the target name

tname = (n_elements(n) gt 0) ? (*stacki[n-1]).tname : (*loadedi).tname
return, tname

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getpol1,stack=n,help=help

; Get the polarization channel for the loaded single-channel spectrum,
; or for stack1 spectrum #n if the stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'value = getpol1([stack=n][,/help])'
  print,' '
  print,'Returns the polarization channel for the loaded single-channel spectrum,'
  print,'or else for stack pair #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in getpol1: The single-channel stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getpol1: Must have n >= 1'
    return, ''
  endif else if n gt nstack1 then begin
    print,'ERROR in getpol1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in getpol1: No single-channel spectrum is loaded'
  return, ''
endif

; Return the polarization channel

pol = (n_elements(n) gt 0) ? (*stack1[n-1]).pol : (*loaded1).pol
return, pol

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getpoli,stack=n,help=help

; Returns the polarization channel for the loaded image,
; or else for stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() ne 0 then begin
  print,' '
  print,'value = getpoli([stack=n][,/help])'
  print,' '
  print,'Returns the polarization channel for the loaded image, or else for'
  print,'stacki image #n if the stack keyword is set'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in getpoli: The image stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getpoli: Must have n >= 1'
    return, ''
  endif else if n gt nstacki then begin
    print,'ERROR in getpoli: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return, ''
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in getpoli: No image is loaded'
  return, ''
endif

; Return the polarization channel

pol = (n_elements(n) gt 0) ? (*stacki[n-1]).pol : (*loadedi).pol
return, pol

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getgroup,stack=n,help=help

; Get the group number for stack pair #n
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'value = getgroup(stack=n[,/help])'
  print,' '
  print,'Returns the group number for stack pair #n'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in getgroup: The pair stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getgroup: Must have n >= 1'
    return, ''
  endif else if n gt nstack then begin
    print,'ERROR in getgroup: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return, ''
  endif
endif else begin
  print,'ERROR in getgroup: Must use the stack keyword'
  return, ''
endelse

; Return the group number

return, (*stack[n-1]).group

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getgroup1,stack=n,help=help

; Get the group number for stack1 spectrum #n
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'value = getgroup1(stack=n[,/help])'
  print,' '
  print,'Returns the group number for stack1 spectrum #n'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in getgroup1: The single-channel stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getgroup1: Must have n >= 1'
    return, ''
  endif else if n gt nstack1 then begin
    print,'ERROR in getgroup1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return, ''
  endif
endif else begin
  print,'ERROR in getgroup1: Must use the stack keyword'
  return, ''
endelse

; Return the group number

return, (*stack1[n-1]).group

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getgroupi,stack=n,help=help

; Get the group number for stacki image #n
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'value = getgroupi(stack=n[,/help])'
  print,' '
  print,'Returns the group number for stacki image #n'
  print,' '
  return, ''
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in getgroupi: The image stack is empty'
    return, ''
  endif else if n le 0 then begin
    print,'ERROR in getgroupi: Must have n >= 1'
    return, ''
  endif else if n gt nstacki then begin
    print,'ERROR in getgroupi: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return, ''
  endif
endif else begin
  print,'ERROR in getgroupi: Must use the stack keyword'
  return, ''
endelse

; Return the group number

return, (*stacki[n-1]).group

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function beq_sf,n,B,help=help

; Computes the equivalent bandwidth of an "S(f)" function:
;   
;    S(f) = amplitude * [ 1 - (2f/B)^2 ]^(n/2)  ;  |f| < B/2
;
; Note that n needn't be an integer.
;
; If zero-crossing bandwidth B is given (in Hz), beq_sf returns
; the equivalent bandwidth (in Hz).  Otherwise it returns
; the equivalent bandwidth as a fraction of B.
;
; 2013 Jun 11: use !DPI instead of !PI so that the return value is double-precision

if n_params() eq 0 then n = -1
if n_params() le 1 then B = 1

if n lt 0 or B le 0 or n_params() gt 2 or keyword_set(help) then begin
  print,' '
  print,'beq = beq_sf(n[,B][,/help])'
  print,' '
  print,'For any real n >= 0 and B > 0, returns the equivalent bandwidth of the function'
  print,' '
  print,'     S(f) = amplitude * [ 1 - (2f/B)^2 ]^(n/2)  ;  |f| < B/2'
  print,' '
  print,'If zero-crossing bandwidth B is omitted, returned value is Beq/B'
  print,' '
  return,-1
endif

; Work with the logarithm of the gamma function to avoid overflow

beq_over_B = (0.5*sqrt(!DPI)) * exp(  lngamma(n + 1.5) + 2*lngamma(n/2.0 + 1.0)   $
                                    - lngamma(n + 1.0) - 2*lngamma(n/2.0 + 1.5) )

return, B*beq_over_B
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function makefilter,efb,n,dfreq,ndata,efb_use,help=help

; Create a frequency filter with a specified effective resolution efb
;
; The filter shape is the same as the shape of a radar echo from a
; spherical target with a (cos theta)^n scattering law:
;
;    S(f) = amplitude * [ 1 - (2f/B)^2 ]^(n/2)  ;  |f| < B/2
;
; Here B is the full width of the filter; the amplitude is chosen so as
; to normalize the filter.  B is chosen such that the filter's equivalent
; bandwidth is equal to efb; this must be done iteratively due to the
; finite frequency resolution (dfreq).
;
; Note that power-law exponent n needn't be an integer.
;
; Cribbed from the tkplay program filter.cc.  I've changed some
; variable names, have returned the entire filter rather than the
; positive-frequency portion (so that convolution can be done via the
; IDL "convol" function), and have normalized the filter.  I also get
; the equivalent bandwidth of a continuous "S(f)" function exactly
; (via the gamma function) rather than by Gaussian quadrature.  Finally,
; I let n be an input parameter rather than fixing it at 2.0.
;
; 2003 Nov 21: Modified to add the "efb_use" output parameter, the
;              actual value of efb used for the filter
;              (close but not quite equal to efb)
; 2013 Jun 11: Make frequency array double-precision so that calculations
;                  are done in double-precision when iterating to get
;                  filter efb

maxiter = 5000L
tolerance = 0.001

if n_params() ne 5 or keyword_set(help) then begin
  print,' '
  print,'filter = makefilter(efb,n,dfreq,ndata,efb_use[,/help])'
  print,' '
  print,'Create a filter with frequency resolution (bin width) dfreq and'
  print,'     effective frequency resolution (equivalent bandwidth) efb,'
  print,'     to be convolved with a spectrum in order to smooth it.'
  print,' '
  print,'     The filter shape is the same as the shape of echoes from a'
  print,'     spherical target with a (cos theta)^n scattering law:'
  print,' '
  print,'          S(f) = amplitude * [ 1 - (2f/B)^2 ]^(n/2)  ;  |f| < B/2'
  print,' '
  print,'     Here B is the full width of the filter; the amplitude is chosen'
  print,"     so as to normalize the filter to unit sum.  B is chosen such that'
  print,"     the filter's equivalent bandwidth is equal to efb."
  print,' '
  print,"     Note that power-law exponent n needn't be an integer."
  print,' '
  print,'     dfreq is the frequency resolution (bin width) of the filter,'
  print,'          which should be the same as that of the spectrum which is'
  print,'          to be smoothed using this filter'
  print,' '
  print,'     ndata is the length of the spectrum to be smoothed:'
  print,'          makefilter checks that the filter it returns is narrower than'
  print,'          the spectrum with which it will be convolved.'
  print,' '
  print,'     efb_use is an output parameter, the actual value of efb which is'
  print,'          used to create the filter (close but not quite equal to efb)'
  print,' '
  efb_use = dfreq
  return, [0]
endif else if ndata lt 4 then begin
  print,'ERROR in makefilter: Must have at least 4 spectral points'
  efb_use = dfreq
  return, [0]
endif else if efb le 0 then begin
  print,'ERROR in makefilter: Requested effective resolution must be positive'
  efb_use = dfreq
  return, [0]
endif else if dfreq le 0 then begin
  print,'ERROR in makefilter: Raw frequency resolution must be positive'
  efb_use = dfreq
  return, [0]
endif else if abs(efb/dfreq - 1.0) le tolerance then begin
  efb_use = dfreq
  return, [1.0]
endif else if efb lt dfreq then begin
  print,'ERROR in makefilter: Requested effective resolution (',efb, $
        ' Hz) is finer than the raw resolution (',dfreq,' Hz)', $
        format='(a,f7.3,a,f7.3,a)'
  efb_use = dfreq
  return, [0]
endif else if n lt 0 then begin
  print,"ERROR in makefilter: Power-law exponent n can't be negative"
  efb_use = dfreq
  return, [0]
endif

; Generate a starting approximation to the proper filter width.  For this
; we must first compute the equivalent bandwidth of an "S(f)" filter with
; continuous rather than discrete frequency.
;
; hwidth (B/2) is the half the zero-crossing filter width; we'll adjust it
; later for finite frequency resolution

beq_over_b = beq_sf(n)
hwidth = 0.5*efb/beq_over_b

; Compute the maximum allowed length for the positive-frequency portion of the filter

pos_maxsize = floor(ndata)/2 - 1

; Make sure to use a filter whose full bandwidth is narrower than the spectrum

width_check = long(hwidth/dfreq)
if width_check gt pos_maxsize then begin
  print,'ERROR in makefilter: Requested filter (',2*width_check,       $
        ' elements) has width >= spectrum width (',ndata,' elements)', $
        format='(a,i0,a,i0,a)'
  efb_use = dfreq
  return, [0]
endif

; Now iterate to find a filter efb approximately the same size as the efb requested

trial_efb = efb
iter = 0L
repeat begin
  iter = iter + 1L 
  npos = (ceil(hwidth/dfreq) - 1) < pos_maxsize
  npos = npos > 1L
  freq = (dindgen(npos) + 1)*dfreq              ; just the positive frequencies
  pos_filter = (1.0 - (freq/hwidth)^2)^(n/2.0)  ; the positive-frequency part of our filter
  fnorm = 1.0 + 2*total(pos_filter)
  fsq = 1.0 + 2*total(pos_filter^2)
  trial_efb = dfreq*(fnorm^2)/fsq
  hwidth = hwidth*(efb/trial_efb)^(0.1 - (0.05*iter)/maxiter)
endrep until ( (iter eq maxiter) or ( abs(trial_efb - efb)/dfreq le tolerance ) )

if iter eq maxiter then begin
  print,'WARNING in makefilter: Reached ',maxiter,'-iteration limit',format='(a,i0,a)'
  print,'                       Will use effective resolution ',trial_efb, $
        ' Hz instead of ',efb,' Hz',format='(a,d10.6,a,d10.6,a)'
endif

; Construct the normalized filter

nfilter = 2*npos + 1L
filter = dblarr(nfilter)
centerbin = nfilter/2
filter[centerbin] = 1.0/fnorm
filter[centerbin+1:nfilter-1] = pos_filter/fnorm
filter[0:centerbin-1] = reverse(filter[centerbin+1:nfilter-1])

efb_use = trial_efb
 
return, float(filter)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function gaussian, xi, parms, pderiv, DOUBLE=double

;**********************************
; This routine is copied from the IDL Astronomy Library
;     http://idlastro.gsfc.nasa.gov
;**********************************

;+
; NAME:
;       GAUSSIAN
; PURPOSE:
;       Compute the 1-d Gaussian function and optionally the derivative
; EXPLANATION:
;       Compute the 1-D Gaussian function and optionally the derivative 
;       at an array of points.
;
; CALLING SEQUENCE:
;       y = gaussian( xi, parms,[ pderiv ])
;
; INPUTS:
;       xi = array, independent variable of Gaussian function.
;
;       parms = parameters of Gaussian, 2, 3 or 4 element array:
;               parms[0] = maximum value (factor) of Gaussian,
;               parms[1] = mean value (center) of Gaussian,
;               parms[2] = standard deviation (sigma) of Gaussian.
;               (if parms has only 2 elements then sigma taken from previous
;               call to gaussian(), which is stored in a  common block).
;               parms[3] = optional, constant offset added to Gaussian.
; OUTPUT:
;       y -  Function returns array of Gaussian evaluated at xi.    Values will
;            be floating pt. (even if xi is double) unless the /DOUBLE keyword
;            is set.
;
; OPTIONAL INPUT:
;       /DOUBLE - set this keyword to return double precision for both
;             the function values and (optionally) the partial derivatives.
; OPTIONAL OUTPUT:
;       pderiv = [N,3] or [N,4] output array of partial derivatives,
;               computed only if parameter is present in call.
;
;               pderiv[*,i] = partial derivative at all xi absisca values
;               with respect to parms[i], i=0,1,2,[3].
;
;
; EXAMPLE:
;       Evaulate a Gaussian centered at x=0, with sigma=1, and a peak value
;       of 10 at the points 0.5 and 1.5.   Also compute the derivative
;
;       IDL> f = gaussian( [0.5,1.5], [10,0,1], DERIV )
;       ==> f= [8.825,3.25].   DERIV will be a 2 x 3 array containing the
;       numerical derivative at the two points with respect to the 3 parameters.
; 
; COMMON BLOCKS:
;       None
; HISTORY:
;       Written, Frank Varosi NASA/GSFC 1992.
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use machar() for machine precision, added /DOUBLE keyword,
;       add optional constant 4th parameter    W. Landsman   November 2001
;-
  On_error,2
  common gaussian, sigma

  if N_params() LT 2 then begin
        print,'Syntax - y = GAUSSIAN( xi, parms,[ pderiv, /DOUBLE ])'
        print,'         parms[0] = maximum value (factor) of Gaussian'
        print,'         parms[1] = mean value (center) of Gaussian'
        print,'         parms[2] = standard deviation (sigma) of Gaussian'
        print,'         parms[3] = optional constant to be added to Gaussian'
        return, -1
  endif

  common gaussian, sigma

        Nparmg = N_elements( parms )
        npts = N_elements(xi) 
        ptype = size(parms,/type)
        if (ptype LE 3) or (ptype GE 12) then parms = float(parms)
        if (Nparmg GE 3) then sigma = parms[2]

        double = keyword_set(DOUBLE)
        if double then $       ;Double precision?
            gauss = dblarr( npts ) else $
            gauss = fltarr( npts )
 
        z = ( xi - parms[1] )/sigma
        zz = z*z

; Get smallest value expressible on computer.   Set lower values to 0 to avoid
; floating underflow
        minexp = alog((machar(DOUBLE=double)).xmin)     
 
        w = where( zz LT -2*minexp, nw )
        if (nw GT 0) then gauss[w] = exp( -zz[w] / 2 )

        if N_params() GE 3 then begin

                if double then $ 
                pderiv = dblarr( npts, Nparmg ) else $
                pderiv = fltarr( npts, Nparmg )
                fsig = parms[0] / sigma

                pderiv[0,0] = gauss
                pderiv[0,1] = gauss * z * fsig

                if (Nparmg GE 3) then  pderiv[0,2] = gauss * zz * fsig
                if (Nparmg GE 4) then  pderiv[0,3] = replicate(1, npts)
           endif

 if Nparmg LT 4 then return, parms[0] * gauss else $
                     return, parms[0] * gauss + parms[3]
 end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function psf_gaussian, parameters, NPIXEL=npixel, NDIMENSION=ndim, FWHM=fwhm,  $
                        DOUBLE = double, CENTROID=cntrd, ST_DEV=st_dev,  $
                        XY_CORREL=xy_corr, NORMALIZE=normalize

;**********************************
; This routine is copied from the IDL Astronomy Library
;     http://idlastro.gsfc.nasa.gov
;**********************************

;+
; NAME:
;       PSF_GAUSSIAN
;
; PURPOSE:
;       Create a 1-d, 2-d, or 3-d Gaussian with specified FWHM, center 
; EXPLANATION:
;       Return a point spread function having Gaussian profiles,
;       as either a 1D vector, a 2D image, or 3D volumetric-data.
;
; CALLING SEQUENCE:
;       psf = psf_Gaussian( NPIXEL=, FWHM= , CENTROID = 
;                     [ /DOUBLE, /NORMALIZE, ST_DEV=,  NDIMEN= ] ) 
; or:
;       psf = psf_Gaussian( parameters, NPIXEL = ,NDIMEN = )
;
; REQUIRED INPUT KEYWORD:
;       NPIXEL = number pixels for each dimension, specify as an array,
;               or just one number to make all sizes equal.
;
; OPTIONAL KEYWORDS:
;       CENTROID = floating scalar or vector giving position of  PSF center.    
;               default is exact center of requested vector/image/volume.
;               The number of elements in CENTROID should equal the number of
;               dimensions.    **The definition of Centroid was changed in
;               March 2002, and now an integer defines the center of a pixel.**
;
;       /DOUBLE  = If set, then the output array is computed in double precision
;               the default is to return a floating point array.
;
;       FWHM = the desired Full-Width Half-Max (pixels) in each dimension,
;               specify as an array, or single number to make all the same.
;
;       NDIMEN = integer dimension of result: either 1 (vector), 2 (image), or 
;                3 (volume), default = 2 (an image result).
;
;       /NORMALIZE causes resulting PSF to be normalized so Total( psf ) = 1.
;
;       ST_DEV = optional way to specify width by standard deviation param.
;                Ignored if FWHM is specified.
;
;       XY_CORREL = scalar between 0 and 1 specifying correlation coefficient
;               Use this keyword, for example, to specify an elliptical 
;               Gaussian oriented at an angle to the X,Y axis.   Only valid
;               for 2-dimensional case.
;
;
; INPUTS (optional):
;
;       parameters = an NDIMEN by 3 array giving for each dimension:
;                       [ maxval, center, st_dev ],  overrides other keywords.
;
; EXAMPLE:
;       (1) Create a 31 x 31 array containing a normalized centered Gaussian 
;       with an X FWHM = 4.3 and a Y FWHM = 3.6
;
;       IDL> array = PSF_GAUSSIAN( Npixel=31, FWHM=[4.3,3.6], /NORMAL )
;
;       (2) Create a 50 pixel 1-d Gaussian vector with a maximum of 12, 
;          centered at  pixel 23 with a sigma of 19.2
;
;       IDL> psf = psf_gaussian([12,23,19.2],npixel=50)
; EXTERNAL CALLS:
;       function Gaussian()
; NOTES:
;       To improve speed, floating underflow exceptions are suppressed (using 
;       the MASK=32  keyword of CHECK_MATH() rather than being flagged.
;
; HISTORY:
;       Written, Frank Varosi NASA/GSFC 1991.
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Suppress underflow messages, add DOUBLE keyword. **Modified centroid
;       definition so integer position is pixel center** W. Landsman March 2002
;       Allow use of the ST_DEV (not STDEV) keyword W. Landsman Nov. 2002
;       Do not modify NPIXEL input keyword   W. Landsman  
;-
        On_error,2
	compile_opt idl2

        if (N_params() LT 1 ) and $
            not (keyword_set( FWHM) or keyword_set(ST_DEV)) then begin
                print,'Syntax - psf = PSF_GAUSSIAN( parameters, NPIXEL = )'
                print, $
       'or       psf = PSF_GAUSSIAN( FWHM = ,ST_DEV = ,NPIXEL = ,[CENTROID = ])'
                return, -1
        endif

        sp = size( parameters )
        if sp[0] EQ 1 then begin               ;Vector supplied?
                ndim = 1
                factor = parameters[0]
                cntrd = parameters[1]
                st_dev = parameters[2] 
         endif  else  if (sp[0] GE 1) then begin    ;Ndimen x 3 array supplied?
                 ndim = sp[1]
                 factor = total( parameters[*,0] )/float( ndim )
                cntrd = parameters[*,1]
                st_dev = parameters[*,2]
           endif

        double = keyword_set(double)
        if double then idltype = 5 else idltype = 4
        if N_elements( ndim ) NE 1 then ndim=2
        ndim = ndim>1

        if N_elements( npixel ) LE 0 then begin
                message,"must specify size of result with NPIX=",/INFO
                return,(-1)
          endif else begin 
	      npix = npixel
	      if N_elements( npix ) LT ndim then npix = replicate( npix[0], ndim )
         endelse

        if (N_elements( cntrd ) LT ndim) AND (N_elements( cntrd ) GT 0) then $
                        cntrd = replicate( cntrd[0], ndim )

        if N_elements( cntrd ) LE 0 then cntrd=(npix-1)/2. 
        if N_elements( fwhm ) GT 0 then begin 
               st_dev = fwhm/( 2.0d* sqrt( 2.0d* aLog(2.0d) ) )
               if not double then st_dev  = float(st_dev)
        endif 

        if N_elements( st_dev ) LE 0 then begin
                message,"must specify ST_DEV= or FWHM=",/INFO
                return,(-1)
          endif

        if N_elements( st_dev ) LT ndim then $
                        st_dev = replicate( st_dev[0], ndim )

        CASE ndim OF

        1: BEGIN
                x = findgen( npix[0] ) - cntrd[0]
                psf = gaussian( x, [1,0,st_dev] )
             END

        2: BEGIN
                psf = make_array( DIM=npix[0:ndim-1], TYPE = idltype )
                x = make_array( npix[0], /INDEX, TYPE=idltype ) - cntrd[0]
                y = make_array( npix[1], /INDEX, TYPE=idltype ) - cntrd[1]

                if N_elements( xy_corr ) EQ 1 then begin
                        sigfac = 1 / (2. * st_dev^2 )
                        y2 = sigfac[1] * y^2
                        x1 = sigfac[0] * x
                        yc = y * ( xy_corr/(st_dev[0]*st_dev[1]) )
                        for j=0,npix[1]-1 do begin
                                zz = x * (yc[j] + x1) + y2[j]
                                w = where( zz LT 86, nw )
                                if (nw GT 0) then psf[w,j] = exp( -zz[w] )
                          endfor
                  endif else begin
                        psfx = gaussian( x, [ 1, 0, st_dev[0] ], DOUBLE=double )
                        psfy = gaussian( y, [ 1, 0, st_dev[1] ], DOUBLE=double )
                        error = check_math(/print, MASK=32)
                        save_except = !EXCEPT & !EXCEPT = 0
                        for j=0,npix[1]-1 do psf[0,j] = psfx * psfy[j]
                        error = check_math(MASK=32)    ;Clear floating underflow
                        !EXCEPT = save_except  
                   endelse
             END

        3: BEGIN
                psf = make_array( DIM=npix[0:ndim-1], TYPE = idltype )
                x = make_array( npix[0], /INDEX, TYPE=idltype ) - cntrd[0]
                y = make_array( npix[1], /INDEX, TYPE=idltype ) - cntrd[1]
                z = make_array( npix[2], /INDEX, TYPE=idltype ) - cntrd[2]
                psfx = gaussian( x, [ 1, 0, st_dev[0] ], DOUBLE = double )
                psfy = gaussian( y, [ 1, 0, st_dev[1] ], DOUBLE = double)
                psfz = gaussian( z, [ 1, 0, st_dev[2] ], DOUBLE = double )
                error = check_math(MASK=32,/PRINT)
                save_except = !EXCEPT & !EXCEPT = 0
                for k=0,npix[2]-1 do begin
                    for j=0,npix[1]-1 do psf[0,j,k] = psfx * psfy[j] * psfz[k]
                 endfor
                 error = check_math(MASK=32)
                 !EXCEPT = save_except  
             END

        ENDCASE

        if keyword_set( normalize ) then return, psf/total( psf )

        if N_elements( factor ) EQ 1 then begin
                if (factor NE 1) then return,factor*psf else return,psf
           endif else return, psf
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function gaussfold_cm, lam, flux, fwhm, lammin=lammin, lammax=lammax

;**********************************
; This routine was adapted from
;     http://astro.uni-tuebingen.de/software/idl/aitlib/misc/gaussfold.pro
;     by suppressing floating-point underflow messages and by making two changes
;     relevant for very coarse smoothing: ensuring that the resolution of the
;     flux interpolation grid is at least as fine as the input resolution; and
;     zero-filling the interpolated flux array if it would otherwise be shorter
;     than the smoothing kernel
;**********************************

;+
; NAME:     
;           gaussfold
;
;
; PURPOSE:
;           Smoothes a plot by convolving with a Gaussian profile.
;           Main purpose is to convolve a spectrum (flux against
;           wavelength) with a given instrument resolution.
;           Also applicable e.g. to smooth ligthcurves or in
;           time-series analysis. 
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;           smoothedFlux = gaussfold(lambda, flux, sigma,
;                                    LAMMIN=lammin, LAMMAX=lammax) 
;
;
; INPUTS:
;           lambda = In ascending order sorted array (double or float)
;                    containing the values of the x-axis (e.g. the
;                    wavelength- or frequencygrid). Irregularly spaced
;                    grids are acepted. 
;           flux   = Array (double or float) same size as lambda
;                    containing the values of the y-axis (e.g. the
;                    flux).  
;           fwhm   = FWHM in units of the x-axis (e.g. Angstroem) of
;                    the Gaussian profile.
;
;
; OPTIONAL INPUTS:
;           NONE
;
;
; KEYWORD PARAMETERS:
;           lammin &
;           lammax   = Defines the x-axis range. All y-data within
;                      this range will be smoothed.  
;                      CAUTION: improtant in case of large arrays (memory!)
;                      DEFAULTS: lammin = MIN(lambda) i.e. lambda[0]
;                                lammax = MAX(lambda) i.e. lambda[N_ELEMENTS(lambda)-1]
;
;
; OUTPUTS:
;           smoothedFlux = Array (double of float) same size as lambda
;                          containing the smoothed flux.  
;
;
; OPTIONAL OUTPUTS:
;           NONE
;
;
; COMMON BLOCKS:
;           NONE
;
;
; SIDE EFFECTS:
;           NONE
;
;
; RESTRICTIONS:
;           NONE
;
;
; PROCEDURE:
;           1. Interpolation of the flux on a fine spaced grid.
;              Oversamplingfaktor: 17
;           2. Convolution of the flux with a gaussian profile 
;           3. Interpolation of the smoothed flux on the original grid.
;
;
; EXAMPLE:
;           fluxS = GAUSSFOLD( lambda, flux, 0.5)
;
;
; MODIFICATION HISTORY:
;           Version 1.0  by Katja Pottschmidt  1999
;           Version 2.0  by Jochen Deetjen     1999
;                        - Documentation added
;                        - Keyword defaults defined
;                        - convolution with the gaussian kernel is
;                          restricted within a window 
;                          -> by far less time consuming 
;
;-

   IF (NOT KEYWORD_SET(lammin)) THEN $
     lammin = MIN(lam)

   IF (NOT KEYWORD_SET(lammax)) THEN $
     lammax = MAX(lam)

   ; Prepare to interpolate the flux onto a fine grid, 17 times
   ; finer than the FWHM of the Gaussian smoothing filter --
   ; unless the input resolution is even finer than that

   nlam = n_elements(lam)
   lamincr = lam[1:nlam-1] - lam[0:nlam-2]
   inputres = min(abs(lamincr))
   dlambda   = (fwhm / 17D0) < inputres

   ;; get a 1D gaussian profile (and suppress underflow messages)

   fwhm_pix = fwhm / dlambda
   window   = FIX( 17 * fwhm_pix )

   error = check_math(/print, MASK=32)
   save_except = !EXCEPT & !EXCEPT = 0
   gauss = PSF_GAUSSIAN( NP=window, FWHM=fwhm_pix, /NORMALIZE, NDIMEN=1 )
   error = check_math(MASK=32)    ;Clear floating underflow
   !EXCEPT = save_except  

   ; Now interpolate the flux, zero-filling if necessary to make sure
   ; that there are at least as many flux values as smoothing kernel elements
   ; (or else the convol function will choke).  Note that the interpolated
   ; spectrum will have the x-axis in increasing order even if the input
   ; spectrum is in decreasing order.

   interlam  = lammin + dlambda * DINDGEN( LONG( (lammax-lammin)/dlambda+1 ))
   interflux = INTERPOL( flux, lam, interlam )   
   n_interp = n_elements(interlam)
   n_extra = (n_elements(gauss) - n_interp) > 0L
   if n_extra gt 0 then begin
     n_half = (n_extra + 1)/2
     interlam  = (lammin - dlambda*n_half) + dlambda*dindgen(n_interp + 2*n_half)
     interflux = [fltarr(n_half),interflux,fltarr(n_half)]
   endif

   ;; convolve input spectrum with the gauss profile
   fold  = CONVOL( interflux, gauss, /CENTER, /EDGE_TRUNCATE)

   ; Interpolate in order to get back to the input resolution
   ; (and, if the input spectrum had the x-axis in decreasing order,
   ; to go back to decreasing order)

   fluxfold  = INTERPOL( fold, interlam, lam )

   RETURN, fluxfold

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function bandwidth1,threshold,leftedge,rightedge,mode=mode,center=center,help=help

; Get the spectral bandwidth of the loaded single-channel spectrum,
; measured across a contiguous set of points all of which are at least
; threshold standard deviations strong.  The default level is threshold = 0,
; that is, the zero-crossing bandwidth.  Linear interpolation is used to 
; estimate the "exact" spectral edges.
;
; The mode parameter defines the starting point for the search:
;   0: xjcen
;   1: peak signal within the defined signal bins, jsnr1-jsnr2 (default)
;   2: frequency specified via the center keyword
;   3: frequency bin specified via the center keyword
;   4: try all points within the signal bins and use the maximum width
;
; The width is returned, and the parameters leftedge and rightedge
; (if included in the procedure call) are set to the frequencies of the
; left and right edges, respectively.  Note that bandwidth1 is positive
; even if frequency increases from right to left.

common loadedBlock,loadedi,loaded1,loaded

if n_params() eq 0 then threshold = 0.0
if n_elements(mode) eq 0 then mode = 1

if keyword_set(help) or n_params() gt 3 or mode lt 0 or mode gt 4 then begin
  print,' '
  print,'bw = bandwidth1([,threshold[,leftedge,rightedge]][,mode=mode,center=center][,/help])'
  print,'  threshold is the number of sigmas for the crossing points (default = 0)'
  print,'  leftedge, rightedge optionally return the frequencies of the left and right edges'
  print,'  mode defines the starting point for the search:'
  print,'       0: xjcen'
  print,'       1: peak signal within the defined signal bins, jsnr1-jsnr2 (default)'
  print,'       2: frequency specified via the center keyword'
  print,'       3: frequency bin (0-based) specified via the center keyword'
  print,'       4: try all points within the signal bins and use the maximum width'
  print,' '
  return, -999.9
endif

if (*loaded1).ndata le 2 then begin
  print,"ERROR in bandwidth1: No single-channel spectrum is loaded, so bandwidth can't be computed"
  return, -999.9
endif

; Get some elements of the loaded spectrum

spec = (*loaded1).spec
f = (*loaded1).freq
tags1 = (*loaded1).tags
lim1 = tags1.jsnr1
lim2 = tags1.jsnr2
deltaf = f[1] - f[0]  ; < 0 if frequency increases leftward
posfr = tags1.posfr   ;  -1 if frequency increases leftward, else +1

; Use the mode parameter to determine how to define the bandwidth

if mode ne 4 then begin

  case mode of
    0: centerbin = tags1.xjcen
    1: begin
         peak = max(spec[lim1:lim2], centerbin)
         centerbin = lim1 + centerbin
       endcase
    2: centerbin = round( (center - f[0])/deltaf )
    3: centerbin = center
  endcase

  findedges,threshold,centerbin,leftedge,rightedge

endif else begin

  ; Use the maximum-width contiguous set of bins which all
  ; are above the threshold.  Insist that at least one bin
  ; lies within the signal limits [jsnr1, jsnr2], but allow
  ; the set to extend beyond these limits.

  centerbin = lim1
  maxwidth = -999.9
  repeat begin
    findedges,threshold,centerbin,l_edge,r_edge
    width = posfr*(r_edge - l_edge)
    if width gt maxwidth then begin
      leftedge = l_edge
      rightedge = r_edge
      maxwidth = width
    endif
    centerbin = round( (r_edge - f[0])/deltaf ) + 1L
  endrep until centerbin gt lim2

endelse

return, posfr*(rightedge - leftedge) > 0L
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro findedges,threshold,centerbin,leftedge,rightedge,help=help

; Use linear interpolation to estimate the frequencies of the left and right
; spectral edges: start from bin centerbin and work outward in both
; directions until the innermost crossing points (where the signal first
; drops below threshold standard deviations) are reached.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 4 or keyword_set(help) then begin
  print,' '
  print,'findedges,threshold,centerbin,leftedge,rightedge[,/help]'
  print,' '
  return
endif

if (*loaded1).ndata le 2 then begin
  print,"ERROR in findedges: No single-channel spectrum is loaded,", $
        " so signal edges can't be determined",format='(2a)'
  return
endif

; Get some elements of the loaded spectrum

f = (*loaded1).freq
spec = (*loaded1).spec
ndata = (*loaded1).ndata

df = f[1] - f[0]  ;  could be positive or negative

; Check that the starting point is within the spectrum

if (centerbin lt 0 or centerbin ge ndata) then begin
  center = f[0] + centerbin*df
  print,' '
  print,'ERROR in findedges: Center bin #',centerbin,' (f = ',center, $
        ' Hz) is out of bounds',format='(a,i0,a,f0,a)'
  print,' '
  leftedge = f[ndata-1] + 0.5*df
  rightedge = f[0] - 0.5*df
  return
endif

; Make sure that the starting point is above the threshold

if spec[centerbin] lt threshold then begin
  leftedge = f[ndata-1] + 0.5*df
  rightedge = f[0] - 0.5*df
  return
endif

; Starting from centerbin, find the first crossing point to the left,
; then get its frequency

leftbin = centerbin
while (leftbin gt 1 and spec[leftbin-1] ge threshold) do begin
  leftbin = leftbin - 1L
endwhile
if (leftbin eq 1 and spec[0] ge threshold) then begin
  print,'WARNING: Signal extends past left edge of spectrum'
  leftedge = f[0] - 0.5*df
endif else begin
  fracbin = (spec[leftbin] - threshold) / (spec[leftbin] - spec[leftbin-1])
  leftedge = f[leftbin] - fracbin*df
endelse

; Starting from centerbin, find the first crossing point to the right,
; then get its frequency

rightbin = centerbin
while (rightbin lt ndata-2 and spec[rightbin+1] ge threshold) do begin
  rightbin = rightbin + 1L
endwhile
if (rightbin eq ndata-2 and spec[ndata-1] ge threshold) then begin
  print,'WARNING: Signal extends past right edge of spectrum'
  rightedge = f[ndata-1] + 0.5*df
endif else begin
  fracbin = (spec[rightbin] - threshold) / (spec[rightbin] - spec[rightbin+1])
  rightedge = f[rightbin] + fracbin*df
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showbw,threshold,mode=mode,center=center,chan=chan,help=help

; Display the signal bandwidth at a given threshold for the loaded pair
 
common loadedBlock,loadedi,loaded1,loaded
common channelBlock, chanstrings, maxchan


if n_params() gt 1 or keyword_set(help) then begin
  print,' '
  print,'showbw[,threshold][,mode=mode[,center=center]][,chan=1 or 2][,/help]'
  print,'  threshold is the number of sigmas for the crossing points (default = 0)'
  print,'  mode defines the starting point for the search:'
  print,'       0: xjcen'
  print,'       1: peak signal within the defined signal bins, jsnr1-jsnr2 (default)'
  print,'       2: frequency specified via the center keyword'
  print,'       3: frequency bin (0-based) specified via the center keyword'
  print,'       4: try all points within the signal bins and use the maximum width'
  print,' '
  return
endif

if (*loaded).ndata le 2 then begin
  print,"ERROR in showbw: No pair is loaded, so signal bandwidth can't be determined"
  return
endif

npol = n_elements((*loaded).tags)

; Decide which channel to use

showpol = intarr(npol)

if n_elements(chan) eq 0 then begin
  showpol = showpol + 1
endif else begin
  if chan-1 gt npol then begin
    print,'Chan may not be larger than the number of polarizations'
    return
  endif
  showpol[chan-1] = 1
endelse

; Set the defaults for parameters which haven't been specified

if n_params() eq 0 then threshold = 0
if n_elements(mode) eq 0 then begin
  mode = 1
  center = 0   ; any defined value for center will do here
endif

; Store the loaded single-channel spectrum so we can use that "space" and
; reload the spectrum later on

storeloaded1 = ptr_new(*loaded1)

; Display the bandwidth(s) and the left and right spectral edges

printstring = ''
for ch=1,npol do begin
  if showpol[ch-1] then begin
    two2one,ch
    bw1 = bandwidth1(threshold,leftedge,rightedge,mode=mode,center=center)
    printstring = printstring + chanstrings[ch-1] +": " $
                  + string(bw1,format='(f7.2)') $
                  + ' Hz  (' + string(leftedge,format='(f8.2)') $
                  + ' Hz to  ' + string(rightedge,format='(f8.2)') $
                  + ' Hz)     '
  endif
endfor
print,' '
print,printstring
print,' '

; Reload the single-channel spectrum that was there at the start

*loaded1 = *storeloaded1
ptr_free,storeloaded1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showbw1,threshold,mode=mode,center=center,help=help
 
; Display the signal bandwidth at a given threshold for the loaded 
; single-channel spectrum
 
common loadedBlock,loadedi,loaded1,loaded

if n_params() gt 1 or keyword_set(help) then begin
  print,' '
  print,'showbw1[,threshold][,mode=mode[,center=center]][,/help]'
  print,'  threshold is the number of sigmas for the crossing points (default = 0)'
  print,'  mode defines the starting point for the search:'
  print,'       0: xjcen'
  print,'       1: peak signal within the defined signal bins, jsnr1-jsnr2 (default)'
  print,'       2: frequency specified via the center keyword'
  print,'       3: frequency bin (0-based) specified via the center keyword'
  print,'       4: try all points within the signal bins and use the maximum width'
  print,' '
  return
endif

if (*loaded1).ndata le 2 then begin
  print,"ERROR in showbw1: No spectrum is loaded, so signal bandwidth can't be determined"
  return
endif

; Set the defaults for parameters which haven't been specified

if n_params() eq 0 then threshold = 0
if n_elements(mode) eq 0 then begin
  mode = 1
  center = 0   ; any defined value for center will do here
endif

; Display the bandwidth(s) and the left and right spectral edges

print,' '
bw1 = bandwidth1(threshold,leftedge,rightedge,mode=mode,center=center)
print,bw1,' Hz  (',leftedge,' Hz to ',rightedge,' Hz)', $
      format='(f7.2,a,f8.2,a,f8.2,a)'
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro xsec1,xsec,xsecErr,cumulative=cumulative,maxf=maxf,xrange=xrange,        $
          noplot=noplot,noerase=noerase,wide_err=wide_err,base_err=base_err, $
          bwsigma=bwsigma,sdev=sdev,set=set,xtitle=xtitle,                   $
          chanstring=chanstring,nopos=nopos,help=help

; Compute the radar cross section of the loaded single-channel spectrum, along
; with its rms uncertainty and several associated parameters
;
; -- see the long help message below for usage details
;
; 2004 Jul 19: Added /wide_err keyword
;
; 2004 Dec 13: Added /cumulative keyword and associated keywords
;
; 2006 Jun 25: frequency refers to center of bin, not left edge
;
; 2022 Jan 12 Added "nopos" if you don't want to truncate at zero

common loadedBlock,loadedi,loaded1,loaded

; Check that frequency range was correctly input (if input at all)

if not keyword_set(cumulative) then begin
  xrangeOK = 1
endif else if n_elements(xrange) gt 0 then begin
  if n_elements(xrange) ne 2 then begin
    xrangeOK = 0
  endif else begin
    xrangeOK = (xrange[0] lt xrange[1])
  endelse
endif else if n_elements(maxf) gt 0 then begin
  if n_elements(maxf) ne 1 then begin
    xrangeOK = 0
  endif else begin
    xrangeOK = (maxf gt 0)
  endelse
endif else begin
  xrangeOK = 1
endelse

if keyword_set(help) or n_params() gt 2 or not xrangeOK then begin
  print,' '
  print,'xsec1[,xsec][,xsecErr][,/cumulative][,xrange=[minf,maxf]][,maxf=maxf]         $'
  print,'     [,/noplot][,/noerase][,/wide_err][,base_err=base_err][,bwsigma=bwsigma]  $'
  print,'     [,sdev=sdev][,xtitle=xtitle][,chanstring=chanstring][,/set][,/help]'
  print,' '
  print,'Estimate the radar cross section of the loaded single-channel spectrum'
  print,'    in one of two ways:'
  print,' '
  print,'    If /cumulative is set, compute the cumulative signal as a function of'
  print,'         frequency and then do a zeroth-order fit to the portion of this curve'
  print,'         beyond the maximum signal frequency (determined by tag jsnr2 or,'
  print,'         if frequency increases leftward, by jsnr1).'
  print,'    Otherwise, just sum the signal between the signal limits (jsnr1-jsnr2).'
  print,' '
  print,'Also get the rms cross-section uncertainty, the equivalent bandwidth, the'
  print,'    bandwidth (zero-crossing or otherwise), and the radar albedo (crudely'
  print,'    assuming a spherical target) for the loaded single-channel spectrum'
  print,' '
  print,'    The cross-section uncertainty is the sum in quadrature of 2-3 contributions,'
  print,'       one or two for noise fluctuations and one for baseline uncertainties:'
  print,' '
  print,'       1) sdev times the square root of the number of channels within'
  print,'              the equivalent bandwidth; if /wide_err is set, use the full'
  print,'              signal width (jsnr1-jsnr2) instead of the equivalent bandwidth'
  print,' '
  print,'       2) the prediction error on the zeroth-order fit (if /cumulative is set)'
  print,' '
  print,'       3) the signal within a box whose width is the signal width (jsnr1-jsnr2)'
  print,'              and whose height is base_err noise standard deviations'
  print,' '
  print,'xsec and xsecErr return the cross section and uncertainty if included'
  print,' '
  print,'/cumulative changes how the cross section is computed (see above)'
  print,' '
  print,'The next four keywords apply only if /cumulative is set:'
  print,' '
  print,'    xrange is the frequency range over which the spectrum is integrated'
  print,'        -- must have maxf > minf even if frequency increases leftward'
  print,'        (default = full frequency range)'
  print,'    maxf=100 is the same as xrange=[-100,100]'
  print,'    /noplot prevents the cumulative curve from being plotted'
  print,'    /noerase prevents the window from being erased before plotting'
  print,' '
  print,'/wide_err causes the cross section uncertainty due to noise fluctuations'
  print,'          (see item 1 above) to be computed using the full signal width'
  print,'          rather than the equivalent bandwidth; this increases the'
  print,'          uncertainty, although it probably is still much smaller'
  print,'          than the systematic calibration uncertainty'
  print,'base_err is how many noise standard deviations to use for the'
  print,'   rms baseline uncertainty when computing xsecErr (see above)'
  print,'   (default = 0)'
  print,'bwsigma is used for bandwidth estimation:'
  print,'   the full bandwidth displayed is the width between the two innermost'
  print,'   points which are bwsigma noise standard deviations above the baseline'
  print,'   (default = 0  -->  zero-crossing baseline)'
  print,'sdev= to input sdev instead of using the stored tag value'
  print,'xtitle is the title used for the x-axis'
  print,"chanstring='OC' or 'SC' is just prepended to screen output"
  print,'/set  to set the cross section and cross section error tags'
  print,' '
  return
endif

if (*loaded1).ndata le 2 then begin
  print,"ERROR in xsec1: No spectrum is loaded, so cross section can't be computed"
  print,' '
  return
endif

if not keyword_set(wide_err) then wide_err = 0
if n_elements(base_err) eq 0 then base_err = 0
if n_elements(bwsigma) eq 0 then bwsigma = 0
if not keyword_set(set) then set = 0
if n_elements(xtitle) eq 0 then xtitle = 'Doppler frequency  (Hz)'

if n_elements(chanstring) eq 0 then begin
  string1 = ''
  string2 = ''
  string3 = ''
endif else begin
  string1 = chanstring + ': '
  string2 = '    '
  string3 = chanstring + ' '
endelse

; Get some elements of the loaded spectrum

ndata = (*loaded1).ndata
f = (*loaded1).freq
s = (*loaded1).spec
tags1 = (*loaded1).tags
dfreq = tags1.dfreq
posfr = tags1.posfr
bin0 = tags1.xjcen
lim1 = tags1.jsnr1
lim2 = tags1.jsnr2
if keyword_set(cumulative) and $
       ((posfr eq  1 and lim2 eq ndata-1) or $
        (posfr eq -1 and lim1 eq 0      )      ) then begin
  print,'ERROR in xsec1: /cumulative is set but there is no noise baseline'
  print,'                beyond the maximum signal frequency'
  print,' '
  return
endif

; If /cumulative is set, get the range of frequency bins
; over which the spectrum will be integrated

if keyword_set(cumulative) then begin
  if n_elements(xrange) gt 0 then begin
    xr = (posfr eq 1) ? xrange : reverse(xrange)
  endif else if n_elements(maxf) gt 0 then begin
    xr = (posfr eq 1) ? [-maxf, maxf] : [maxf, -maxf]
  endif else begin
    xr = [(f[0] - 0.5*posfr*dfreq), (f[ndata-1] + 0.5*posfr*dfreq)]
  endelse
  bin1 = 0L > (bin0 + posfr*round(xr[0]/dfreq))
  bin2 = (ndata - 1) < (bin0 + posfr*round(xr[1]/dfreq))
  if bin1 gt lim2 or bin2 lt lim1 then begin
    print,'ERROR in xsec1: /cumulative is set but the requested frequency range'
    print,'                lies entirely outside the signal limits'
    print,' '
    return
  endif else if (posfr eq  1 and bin2 le lim2) or $
                (posfr eq -1 and bin1 ge lim1) then begin
    print,'ERROR in xsec1: /cumulative is set but the requested frequency range'
    print,'                does not extend beyond the maximum signal frequency'
    print,' '
    return
  endif else if (posfr eq  1 and bin1 gt lim1) or $
                (posfr eq -1 and bin2 lt lim2) then begin
    print,'WARNING in xsec1: /cumulative is set but part of the signal range'
    print,'                  lies outside of the requested frequency range'
    print,' '
  endif
endif

; Allow a temporary resetting of sdev (cross section uncertainty per bin)

if n_elements(sdev) gt 0 then _sdev = sdev else _sdev = tags1.sdev

; Compute the equivalent bandwidth (and demand that it not be negative)

sumSignal = total(s[lim1:lim2])
sumSignalSq = total(s[lim1:lim2]^2)
ibeq = (sumSignal^2)/sumSignalSq  ;  # of frequency bins in equivalent bandwidth
beq = ibeq*dfreq > 0.0            ;  equivalent bandwidth in Hz

; Compute the cross section by using the cumulative curve or else by summing
; the signal withing the signal limits

if keyword_set(cumulative) then begin
  cum_sumSignal = fltarr(ndata)
  in_mask = intarr(ndata)
  if posfr eq 1 then begin
    cum_sumSignal[bin1] = s[bin1]
    for n=bin1+1L,bin2 do cum_sumSignal[n] = cum_sumSignal[n-1] + s[n]
    in_mask[lim2+1:bin2] = 1
  endif else begin
    cum_sumSignal[bin2] = s[bin2]
    for n=bin2-1L,bin1,-1 do cum_sumSignal[n] = cum_sumSignal[n+1] + s[n]
    in_mask[bin1:lim1-1] = 1
  endelse
  cum_sumSignal = cum_sumSignal*_sdev
  xsec = polymaskfit(f,cum_sumSignal,0,in_mask=in_mask,yerror=xsecErr_fit)
endif else begin
  xsec = sumSignal*_sdev
endelse

; If /cumulative is set, plot the cumulative curve (with a horizontal line
; indicating the cross section) unless otherwise requested

if keyword_set(cumulative) and not keyword_set(noplot) then begin

  ; Get the maximum and minimum baseline frequencies and the
  ; minimum frequency of the fitting range.

  if posfr eq 1 then begin
    fmin = f[bin1] - 0.5*dfreq
    fmax = f[bin2] + 0.5*dfreq
    fmin_fit = f[lim2] + 0.5*dfreq
  endif else begin
    fmin = f[bin2] - 0.5*dfreq
    fmax = f[bin1] + 0.5*dfreq
    fmin_fit = f[lim1] + 0.5*dfreq
  endelse

  maxy = max(cum_sumSignal[bin1:bin2])
  miny = min(cum_sumSignal[bin1:bin2])
  yspan = maxy - miny
  maxy = maxy + 0.08*yspan
  miny = miny - 0.08*yspan
  plottitle = 'Cumulative ' + string3 + 'signal'
  if notnull((*loaded1).tname) then plottitle = plottitle + ' for ' + (*loaded1).tname
  if not keyword_set(noerase) then erase
  plot,f,cum_sumSignal,xrange=xr,yrange=[miny,maxy], $
       xtitle=xtitle,ytitle='cross section  (km^2)',title=plottitle,linestyle=1
  oplot,[fmin_fit,fmax],[xsec,xsec]
  oplot,[fmin,fmax],[0,0]
endif

; Make sure that the cross section isn't negative

if not keyword_set(nopos) then xsec = xsec > 0.0

; Compute the cross section uncertainty, accounting both for random statistical
; fluctuations and also for possible systematic (baselining) errors.
;
; These uncertainties will be useful for computing errors on the polarization 
; ratio SC/OC, but they can't be taken seriously on their own without adding
; in the (typically much larger) absolute calibration uncertainties.

if keyword_set(wide_err) then begin
  xsecErr_stat = sqrt(lim2-lim1+1)*_sdev      ; cross-section equivalent of the noise power
                                              ;      within the full signal range
endif else begin
  xsecErr_stat = sqrt(ibeq)*_sdev             ; cross-section equivalent of the noise power
                                              ;      within the equivalent bandwidth
endelse

xsecErr_bline = (lim2-lim1+1)*base_err*_sdev  ; cross-section error if the baseline is off
                                              ;     by base_err noise standard deviations

if keyword_set(cumulative) then begin
  xsecErr = sqrt(xsecErr_stat^2 + xsecErr_fit^2 + xsecErr_bline^2)
endif else begin
  xsecErr = sqrt(xsecErr_stat^2 + xsecErr_bline^2)
endelse

; Get the bandwidth at the signal level specified by bwsigma

bw = bandwidth1(bwsigma,mode=1)

; Compute the radar albedo. Do NOT compute its uncertainty, which is dominated
; by calibration uncertainties, target elongation, and incompleteness of rotation
; phase coverage.  The formal diameter uncertainty might contribute for NEAs;
; the formal uncertainty on cross section will always be insignificant.

diameter = getextra1('diameter')
if isnull(diameter) then diameter = 1.0
area = 0.25 * !pi * diameter^2
albedo = xsec/area

; Set the cross section and cross section error tags if requested

if set then begin
  (*loaded1).tags.cross = xsec
  (*loaded1).tags.crerr = xsecErr
endif

; Output the results

print,string1,'dfreq (Hz)    sdev (km^2)    Beq (Hz)    B (Hz)    ', $
      'xsec (km^2) +/- err (km^2)    D (km)    albedo',format='(3a)'
sdev_format = (_sdev ge 1.0) ? 'f12.4' : 'f12.10'
xsec_format = (xsec ge 1.0 or xsecErr ge 1.0) ? 'f12.3' : 'f12.10'
if xsec lt 0 then xsec_format='f12.9'
formatstring = '(a,1x,f7.3,5x,' + sdev_format + ',3x,f7.2,4x,f7.2,4x,'      $
               + xsec_format + ',4x,' + xsec_format + ',2x,f6.2,5x,f5.3)'
print,string2,dfreq,_sdev,beq,bw,xsec,xsecErr,diameter,albedo, $
      format=formatstring
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sigmamu,cumulative=cumulative,chan=chan,sdev=sdev,set=set,help=help,_extra=_ext

; Print the two cross sections (sigmas) and the polarization ratio (mu)
;
; -- see the long help message below for usage details
;
; 2002 Nov 13: Realized that the "conservative error" algorithm (for assigning
;              cross-section errors for purposes of getting SC/OC) was coded
;              incorrectly.  Switch to using the actual SC and OC cross-section
;              errors, but revise procedure xsec1 to assign those errors more
;              conservatively if desired
;
; 2004 Jul 19: Added /wide_err keyword
;
; 2004 Dec 13: Added /cumulative keyword and associated keywords
;
; 2005 Jan 05: When making the plot required by the /cumulative keyword,
;              for some reason the "linestyle" keyword to oplot is ignored
;              unless the window is at least partly hidden while being drawn,
;              so I hide it, draw it, and then show it.

common loadedBlock,loadedi,loaded1,loaded
common channelBlock, chanstrings, maxchan

if keyword_set(help) or n_params() ne 0 then begin
  print,' '
  print,'sigmamu[,/cumulative][,xrange=[minf,maxf]][,maxf=maxf][,/noplot]       $'
  print,'       [,/wide_err][,base_err=base_err][,bwsigma=bwsigma][,sdev=sdev]  $'
  print,'       [,chan=1 or 2][,/set][,/help]'
  print,' '
  print,'For each specified channel of the loaded pair, estimate the cross section'
  print,'    in one of two ways:'
  print,' '
  print,'    If /cumulative is set, compute the cumulative signal as a function of'
  print,'         frequency and then do a zeroth-order fit to the portion of this curve'
  print,'         beyond the maximum signal frequency (determined by tag jsnr2 or,'
  print,'         if frequency increases leftward, by jsnr1).'
  print,'    Otherwise, just sum the signal between the signal limits (jsnr1-jsnr2).'
  print,' '
  print,'Also get the rms cross-section uncertainty, the equivalent bandwidth, the'
  print,'    bandwidth (zero-crossing or otherwise), and the radar albedo (crudely'
  print,'    assuming a spherical target) for each specified channel.'
  print,' '
  print,'    Cross-section uncertainties are sums in quadrature of 2-3 contributions,'
  print,'       one or two for noise fluctuations and one for baseline uncertainties:'
  print,' '
  print,'       1) sdev times the square root of the number of channels within'
  print,'              the equivalent bandwidth; if /wide_err is set, use the full'
  print,'              signal width (jsnr1-jsnr2) instead of the equivalent bandwidth'
  print,' '
  print,'       2) the prediction error on the zeroth-order fit (if /cumulative is set)'
  print,' '
  print,'       3) the signal within a box whose width is the signal width (jsnr1-jsnr2)'
  print,'              and whose height is base_err noise standard deviations'
  print,' '
  print,'If both channels are included, also compute circular polarization ratio SC/OC'
  print,'    and its uncertainty'
  print,' '
  print,'/cumulative changes how cross sections are computed (see above)'
  print,' '
  print,'The next three keywords apply only if /cumulative is set:'
  print,' '
  print,'    xrange is the frequency range over which the spectrum is integrated'
  print,'        -- must have maxf > minf even if frequency increases leftward'
  print,'        (default = full frequency range)'
  print,'    maxf=100 is the same as xrange=[-100,100]'
  print,'    /noplot prevents the cumulative curve from being plotted'
  print,' '
  print,'/wide_err causes the cross section uncertainties due to noise fluctuations'
  print,'    (see item 1 above) to be computed using the full signal width'
  print,'    rather than the equivalent bandwidth; this increases the'
  print,'    uncertainties, although they probably are still much smaller'
  print,'    than the systematic calibration uncertainty'
  print,'base_err is how many noise standard deviations to use for the rms baseline'
  print,'    uncertainty when computing cross-section uncertainties (see above)'
  print,'    (default = 0)'
  print,'bwsigma is used for bandwidth estimation:'
  print,'    the full bandwidth displayed is the width between the two innermost'
  print,'    points which are bwsigma noise standard deviations above the baseline'
  print,'    (default = 0  -->  zero-crossing bandwidth)'
  print,'sdev= to input sdev for both channels instead of using stored tag values'
  print,'/set  to set the cross section and cross section error tags'
  print,' '
  return
endif

if (*loaded).ndata le 2 then begin
  print,"ERROR in sigmamu: No pair is loaded, so can't compute cross sections and SC/OC"
  print,' '
  return
endif

npol = n_elements((*loaded).spec[*,0])

; Check which channels to process

  usepol = intarr(npol)
if n_elements(chan) eq 0 then begin
  usepol[0:npol-1] = 1
endif else if chan ge 1 and chan le npol then begin
  usepol[chan-1] = 1
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
  return
endelse
xsecs = fltarr(npol)
xsecerrs=fltarr(npol)

; Store the loaded single-channel spectrum so we can use that "space" and
; reload the spectrum later on

storeloaded1 = ptr_new(*loaded1)

; Get tags from the loaded pair

tags = (*loaded).tags

; Use the user-specified sdev or else the tag value.

if n_elements(sdev) gt 0 then begin
  sdevs = fltarr(npol) + sdev
endif else begin
  sdevs = tags[*].sdev
endelse

if not keyword_set(cumulative) then cumulative = 0
if not keyword_set(set) then set = 0

; Set up the plots (if /cumulative is set)

if keyword_set(cumulative) then begin
  pmultisave = !p.multi
  nplots = total(usepol)
  !p.multi = [0,1,nplots]
  erase
  if !d.window ge 0 then wshow,!d.window,0  ; temporary fix: hide the window while drawing
endif

; Compute and display the cross sections

didone = 0
aremore = total(usepol)-1
print,' '
for i = 0, npol-1 do begin
  if usepol[i] then begin
    two2one,i+1
    noerase = didone
    xtitlestring = (aremore) ? '' : 'Doppler frequency  (Hz)'
    nopos = i ge 2 ; negative xsec is fine for RE,IM
    xsec1,MYxsec,MYxsecErr,cumulative=cumulative,sdev=sdevs[i],set=set, $
          xtitle=xtitlestring,chanstring=chanstrings[i],noerase=noerase,nopos=nopos,_extra=_ext
    if n_elements(MYxsec) eq 0 then return
    xsecs[i] = MYxsec
    xsecerrs[i] = MYxsecErr
    if set then one2two
    didone = 1
    aremore = aremore - 1
  endif
endfor

; Compute the circular polarization ratio SC/OC
;
; When using Fieller's theorem (function ratioerr) to get the error on the
; ratio, be somewhat conservative by boosting the smaller of the two
; statistical errors computed by xsec1 to equal the larger error
; (taking into account the different sdev values)

if usepol[0] and usepol[1] then begin

  sdev_ratio = sdevs[1]/sdevs[0]
  OCxsec = xsecs[0]
  SCxsec = xsecs[1]
  OCxsecErr = xsecerrs[0]
  SCxsecErr = xsecerrs[1]
  if SCxsecErr/OCxsecErr lt sdev_ratio then begin
    SCxsecErr = OCxsecErr*sdev_ratio
  endif else begin
    OCxsecErr = SCxsecErr/sdev_ratio
  endelse
  ratioerr,SCxsec,SCxsecErr,OCxsec,OCxsecErr,0,1,ratio,ratlo,rathi
  ratlo = (finite(ratlo)) ? (ratlo > 0.0) : 0.0

  ; Display the polarization ratio

  print,'SC/OC = ',ratio,'  (68% confidence interval = ',ratlo,' to ',rathi,')', $
        format='(a,f5.3,a,f5.3,a,f5.3,a)'
  print,' '
endif

; Reload the single-channel spectrum that was there at the start

*loaded1 = *storeloaded1
ptr_free,storeloaded1

; Reset plotting parameters

if keyword_set(cumulative) then begin
  if !d.window ge 0 then wshow  ; temporary fix: show the window now that the plot is completed
  !p.multi = pmultisave
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ratioerr,x1,err1,x2,err2,corr,t,ratio,ratlo,rathi,help=help

; Uses Fieller's theorem to compute ratio x1/x2 
; and its t-sigma upper and lower bounds
;
; x1 = numerator (err1 is its standard deviation)
; x2 = denominator (err2 is its standard deviation)
; corr = correlation between numerator and denominator
; t = number of standard deviations for the output errors

if n_params() ne 9 or keyword_set(help) then begin
  print,' '
  print,'ratioerr,numerator,numer_err,denominator,denom_err,correlation,t, $'
  print,'         ratio,ratio_lowerlimit,ratio_upperlimit[,/help]'
  print,' '
  print,"Use Fieller's theorem to compute a ratio and its t-sigma ", $
             'confidence interval',format='(2a)'
  print,' '
  return
endif

if x2 eq 0 then begin
  print,'WARNING in ratioerr: denominator = 0, so ratio set to infinity'
  print,' '
  ratio = !VALUES.F_INFINITY
  ratlo = -!VALUES.F_INFINITY
  rathi = !VALUES.F_INFINITY
  return
endif

ratio = x1/x2

v11 = err1^2          ; variance of numerator
v22 = err2^2          ; variance of denominator
v12 = corr*err1*err2  ; covariance of numerator and denominator

g = (t^2)*v22/(x2^2)
sqrt_arg = v11 - 2*ratio*v12 + (ratio^2)*v22 - g*v11*(1 - corr^2)
if sqrt_arg ge 0 then begin
  y = (t/x2)*sqrt(sqrt_arg)
  ratlo = (ratio - g*v12/v22 - y)/(1.0 - g)   ; lower bound
  if g lt 1 then begin
    rathi = (ratio - g*v12/v22 + y)/(1.0 - g) ; upper bound
  endif else begin
    rathi = !VALUES.F_INFINITY
  endelse
endif else begin
  ratlo = -!VALUES.F_INFINITY
  rathi = !VALUES.F_INFINITY
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro xseci,xsec,xsecErr,base_err=base_err,set=set,help=help

common loadedBlock,loadedi,loaded1,loaded

if keyword_set(help) or n_params() gt 2 then begin
  print,' '
  print,'xseci[,xsec][,xsecErr][,base_err=base_err][,/set][,/help]'
  print,' '
  print,'Compute the cross section, its rms uncertainty, the equivalent bandwidth,'
  print,'    and the full bandwidth (zero-crossing or otherwise) for the loaded image'
  print,' '
  print,'    The cross-section uncertainty is the sum in quadrature of two contributions,'
  print,'       one for noise fluctuations and one for baseline uncertainties:'
  print,' '
  print,'       1) sdev times the square root of the number of channels within'
  print,'              the equivalent bandwidth;'
  print,'       2) the signal within a box whose width is the signal width (dopbin1-dopbin2)'
  print,'              and whose height is base_err noise standard deviations'
  print,' '
  print,'    xsec and xsecErr return the cross section and uncertainty if included'
  print,'    base_err is how many noise standard deviations to use for the'
  print,'       rms baseline uncertainty when computing xsecErr (see above)'
  print,'       (default = 0)'
  print,'    /set  to set the cross section and cross section error tags'
  print,' '
  return
endif

if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,"ERROR in xseci: No image is loaded, so cross section can't be computed"
  print,' '
  return
endif

if n_elements(base_err) eq 0 then base_err = 0
if not keyword_set(set) then set = 0

; Get some elements of the loaded image

image = (*loadedi).image
width = (*loadedi).width
height = (*loadedi).height
dopbin1 = getextrai('dopbin1')
dopbin2 = getextrai('dopbin2')
if (isnull(dopbin1) or isnull(dopbin2)) then begin
  print,'WARNING in xseci: Using full Doppler range'
  dopbin1 = 0L
  dopbin2 = width - 1
endif
delbin1 = getextrai('delbin1')
delbin2 = getextrai('delbin2')
if (isnull(delbin1) or isnull(delbin2)) then begin
  print,'WARNING in xseci: Using full delay range'
  delbin1 = 0L
  delbin2 = height - 1
endif
ndop = dopbin2 - dopbin1 + 1
ndel = delbin2 - delbin1 + 1
dfreq = getextrai('fres')
if isnull(dfreq) then begin
  print,'WARNING in xseci: Assuming dfreq = 1.0 Hz'
  dfreq = 1.0
endif
delayunit = getextrai('delayunit')
baud = getextrai('baudlen')
spb = getextrai('samples_per_baud')
if isnull(spb) then spb = 1L
rowsPerBaud = getextrai('rows_per_baud')
if isnull(rowsPerBaud) then rowsPerBaud = spb
if isnull(delayunit) then begin
  delayunit = (notnull(baud)) ? baud/rowsPerBaud : -9.99
endif
if isnull(baud) then baud = -9.99
stride = spb/rowsPerBaud
codemethod = getextrai('codemethod')
if isnull(codemethod) then begin
  print,"WARNING in xseci: Assuming codemethod = 'short'"
  codemethod = 'short'
endif
diameter = getextrai('diameter')
if isnull(diameter) then begin
  print,'WARNING in xseci: Assuming diameter = 1.0 km'
  diameter = 1.0
endif

; Get sdev

sdev = getextrai('sdev')
if isnull(sdev) then begin
  print,'WARNING in xseci: Assuming sdev = 1.0 km^2'
  sdev = 1.0
endif else if sdev le 0 then begin
  print,'WARNING in xseci: Assuming sdev = 1.0 km^2'
  sdev = 1.0
endif

; Compute noise covariances for image rows that are j samples apart in delay
; (= k rows apart where j = k*stride).
;
; The "covar" vector holds the factors that multiply (k Tsys df)^2 / N_looks.
; The expressions used here are asymptotically valid as fft length
; and (for short-code images) code length approach infinity.  Slight Doppler
; variations (for spb > 1) have been ignored, on the assumption that Doppler
; frequencies of interest are small compared to the unaliased bandwidth.

n_covar = (codemethod eq 'long_orig') ? spb : 2*spb - 1
j = indgen(n_covar)
if codemethod eq 'long_orig' then begin
  covar = ( 1 - j/(1.0*spb) )^2
endif else begin
  indexmask = intarr(n_covar)
  indexmask[0:spb-1] = 1
  covar = ( ( (2*spb-j-1)*(2*spb-j)*(2*spb-j+1)              $
               - indexmask*4*(spb-j-1)*(spb-j)*(spb-j+1) )   $
            / (6.0*spb^3) )^2
endelse

; Divide the covariances (covar) by the variance (covar[0]) to
; get the correlation coefficients (corr) for image rows that
; are j samples apart in delay

corr = covar/covar[0]

; In case stride > 1, "collapse" the correlation coefficient vector so that
; it is indexed by the difference in row number rather than by the
; difference in sample number

w = where(j mod stride eq 0, rowcount)
corr_row = corr[w]

; The rms summed noise for a rectangular region m columns x n rows is
;
;   (rms noise for individual pixels)
;   * sqrt( m * [sum over k of (n - |k|)*corr(k*stride)] )
;
; where k counts how many rows separate two pixels within a given colummn
; while j = k*stride counts how many samples separate them.
;
; The individual-pixel rms noise is just sdev, which was already multiplied
; by a "noise reduction" factor = sqrt(covar[0]) when it was computed in
; procedure imagecal.
;
; For the typical case where n >> spb, the rms summed noise is larger than
; the "standard" result (sdev * sqrt(# pixels)) by the factor
; sqrt(sum over k of corr_row(k)).

delfactor_row = (ndel - indgen(rowcount)) > 0L
corrsum = 2*total(delfactor_row*corr_row) - ndel   ; includes k < 0

; Compute the cross section and equivalent bandwidth
; (and demand that neither one is negative)

sumSignal = total(image[dopbin1:dopbin2,delbin1:delbin2])
sumSignalSq = total( total(image[dopbin1:dopbin2,delbin1:delbin2], 2)^2 )
ibeq = (sumSignal^2)/sumSignalSq  ;  # of frequency bins in equivalent bandwidth
beq = ibeq*dfreq > 0.0            ;  equivalent bandwidth in Hz
xsec = sumSignal*sdev > 0.0       ;  radar cross section in km^2

; Compute the cross section uncertainty, accounting both for random statistical
; fluctuations and also for possible systematic (baselining) errors.
;
; These uncertainties can't be taken seriously without adding
; in the (typically much larger) absolute calibration uncertainties.

xsecErr_stat = sqrt(ibeq*corrsum)*sdev               ; cross-section equivalent of the noise power
                                                     ;     within the equivalent bandwidth
xsecErr_bline = ndop*ndel*base_err*sdev              ; cross-section error if the baseline is off
                                                     ;     by base_err noise standard deviations
xsecErr = sqrt(xsecErr_stat^2 + xsecErr_bline^2)

; Compute the radar albedo.  Do NOT compute its uncertainty, which is dominated
; by calibration uncertainties, target elongation, and incompleteness of rotation
; phase coverage.  The formal diameter uncertainty might contribute for NEAs;
; the formal uncertainty on cross section will always be insignificant.

area = 0.25 * !pi * diameter^2
albedo = xsec/area

; Set the cross section and cross section error extra tags if requested

xsec_format2 = (xsec ge 1.0 or xsecErr ge 1.0) ? 'f12.3' : 'f12.10'
xsec_format1 = '(' + xsec_format2 + ')'
if set then begin
  setextrai,'f','cross_sec',string(xsec,format=xsec_format1),         $
            comment='radar cross section in km^2'
  setextrai,'f','cross_sec_err',string(xsecErr,format=xsec_format1),  $
            comment='radar cross section error in km^2'
endif

; Output the results

print,' '
print,'dfreq (Hz)  baud,ddel (us)    sdev (km^2)    Beq (Hz)   ', $
      'xsec (km^2) +/- err (km^2)    D (km)    albedo',format='(2a)'
sdev_format = (sdev ge 1.0) ? 'f12.4' : 'f12.10'
formatstring = '(1x,f7.3,4x,f6.2,1x,f6.2,4x,' + sdev_format + ',3x,f7.2,4x,'      $
               + xsec_format2 + ',4x,' + xsec_format2 + ',2x,f6.2,5x,f5.3)'
print,dfreq,baud,delayunit,sdev,beq,xsec,xsecErr,diameter,albedo,format=formatstring
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showfreq,stack=n,help=help

; Displays frequency information about the loaded pair,
; or else about stack pair n if the stack keyword is set.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showfreq[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in showfreq: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showfreq: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in showfreq: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in showfreq: No pair is loaded'
  return
endif

; Get the structure whose parameters are to be displayed

dispStruc = (n_elements(n) gt 0) ? *stack[n-1] : *loaded
tags = dispStruc.tags
dfreq = tags[0].dfreq
posfr = tags[0].posfr
xjcen = tags[0].xjcen
ndata = dispStruc.ndata
f = (n_elements(n) gt 0) ? posfr*dfreq*(findgen(ndata) - xjcen) : dispStruc.freq

; Display the frequency-related parameters

print,'Frequency: ',ndata,' points at ',dfreq,' Hz resolution  (', $
      f[0]-0.5*posfr*dfreq,' Hz to ',f[ndata-1]+0.5*posfr*dfreq,' Hz)', $
      format='(a,i0,a,f7.3,a,f9.2,a,f9.2,a)'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showfreq1,stack=n,help=help

; Displays frequency information about the loaded single-channel spectrum,
; or else about stack1 spectrum n if the stack keyword is set.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showfreq1[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in showfreq1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showfreq1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in showfreq1: There are only ',nstack1,' spectra in stack1', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in showfreq1: No single-channel spectrum is loaded'
  return
endif

; Get the structure whose parameters are to be displayed

dispStruc = (n_elements(n) gt 0) ? *stack1[n-1] : *loaded1
tags1 = dispStruc.tags
dfreq = tags1.dfreq
posfr = tags1.posfr
xjcen = tags1.xjcen
ndata = dispStruc.ndata
f = (n_elements(n) gt 0) ? posfr*dfreq*(findgen(ndata) - xjcen) : dispStruc.freq

; Display the frequency-related parameters

print,'Frequency: ',ndata,' points at ',dfreq,' Hz resolution  (', $
      f[0]-0.5*posfr*dfreq,' Hz to ',f[ndata-1]+0.5*posfr*dfreq,' Hz)', $
      format='(a,i0,a,f7.3,a,f9.2,a,f9.2,a)'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showdimsi,stack=n,help=help

; Displays delay-Doppler dimensions for the loaded image,
; or else for stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showdimsi[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in showdimsi: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showdimsi: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in showdimsi: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in showdimsi: No image is loaded'
  return
endif

; Get the Doppler and delay vectors

if n_elements(n) gt 0 then begin
  width = (*stacki[n-1]).width
  height = (*stacki[n-1]).height
  dfreq = getextrai('fres',stack=n)
  eph_col = getextrai('eph_col',stack=n)
  eph_row = getextrai('eph_row',stack=n)
  delayunit = getextrai('delayunit',stack=n)
  spb = getextrai('samples_per_baud',stack=n)
  rows_per_baud = getextrai('rows_per_baud',stack=n)
  baud = getextrai('baudlen',stack=n)
endif else begin
  width = (*loadedi).width
  height = (*loadedi).height
  dfreq = getextrai('fres')
  eph_col = getextrai('eph_col')
  eph_row = getextrai('eph_row')
  delayunit = getextrai('delayunit')
  spb = getextrai('samples_per_baud')
  rows_per_baud = getextrai('rows_per_baud')
  baud = getextrai('baudlen')
endelse
if isnull(delayunit) and notnull(baud) then begin
  if isnull(spb) then spb = 1L
  if isnull(rows_per_baud) then rows_per_baud = spb
  delayunit = baud/rows_per_baud
endif
if notnull(eph_col) and notnull(eph_row) and notnull(dfreq) and notnull(delayunit) then begin
  f = dfreq*(findgen(width) - eph_col)       ; Doppler (Hz)
  d = delayunit*(findgen(height) - eph_row)  ; delay (usec)
endif else begin
  print,"ERROR in showdimsi: Need the 'eph_col' and 'eph_row' and 'fres' tags,"
  print,"                    plus either the 'delayunit' or 'baudlen' tag"
  return
endelse

; Display the delay-Doppler parameters

nchars = strlen(string((width > height), format='(i0)'))
formatstring = '(i' + string(nchars,format='(i0)') + ')'
printstring1 = 'Image:     '                                              $
               + string(width,format=formatstring) + ' Doppler cols at '  $
               + string(dfreq,format='(f7.3)') + ' Hz resolution  ('      $
               + string(f[0]-0.5*dfreq,format='(f8.2)') + ' Hz to '       $
               + string(f[width-1]+0.5*dfreq,format='(f8.2)') + ' Hz)'
printstring2 = '           '                                                $
               + string(height,format=formatstring) + ' delay   rows at '   $
               + string(delayunit,format='(f7.3)') + ' us resolution  ('    $
               + string(d[0]-0.5*delayunit,format='(f8.2)') + ' us to '     $
               + string(d[height-1]+0.5*delayunit,format='(f8.2)')          $
               + ' us)'
print,' '
print,printstring1
print,printstring2
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showtags,tagname,chan=chan,stack=n,help=help

; Display tags of the loaded pair, or else of stack pair #n if the
; stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline

if n_params() eq 0 then tagname = 'all'

if n_params() gt 1 or keyword_set(help) then begin
  print,'showtags[,tagname][,chan=1 or 2][,stack=n][,/help]'
  print,' '
  print,'If no CW tag name is given, all CW tags are displayed'
  return
endif else if size(tagname, /type) ne 7 then begin
  print,' '
  print,'showtags,tagname[,chan=1 .. npol][,stack=n][,/help]'
  print,'Make sure that tagname is a quoted string!'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in showtags: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showtags: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in showtags: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in showtags: No pair is loaded'
  return
endif

; Get some elements of the loaded pair or else of stack pair #n

if n_elements(n) eq 0 then begin
  tags = (*loaded).tags
  ntags = (*loaded).ntags
endif else begin
  tags = (*stack[n-1]).tags
  ntags = (*stack[n-1]).ntags
endelse
tagnames = strlowcase(tag_names(tags[0]))
npol = n_elements(tags)

; Check that tagname is a valid tag name

tagname = strlowcase(tagname)
if (tagname ne 'all') then begin
  tagnum = where(tagnames eq tagname, count)
  if (count eq 0) then begin
    print,tagname,' is not a valid tag name'
    return
  endif else begin
    tagnum = tagnum[0]
  endelse
endif

; Check which channels to display

  printpol = intarr(4)
if (n_elements(chan) eq 0) then begin
  printpol[0:npol-1] = 1
endif else if (chan gt 1 and chan le npol) then begin
  printpol[chan] = 1
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or 3 or 4 or omit (both)'
  return
endelse
polstring = ['OC: ','SC: ','RE: ','IM: ']

; Display the tag(s)

if (tagname eq 'all') then begin

  ; print all tag values

  for k0 = 0L,ntags-1,iomaxline do begin
    lastOnLine = (k0 + iomaxline - 1) < (ntags - 1)
    printstring = 'TAG: '
    for k=k0,lastOnLine do begin
      printstring = printstring + string(tagnames[k],format=nameformat)
    endfor
    print,' '
    print,printstring
    for ch=1,npol do begin
      if printpol[ch-1] then begin
        printstring = polstring[ch-1] + ' '
        for k=k0,lastOnLine do begin
          printstring = printstring + string(tags[ch-1].(k),format=tagformats[k])
        endfor
        print,printstring
      endif
    endfor
  endfor
  print,' '

endif else begin

  ; print one tag value

    printstring = tagname + ':'
    for ch=1,npol do begin
      if printpol[ch-1] then begin
        printstring = printstring + '     ' + polstring[ch-1] $
                      + string(tags[ch-1].(tagnum),format=tagformats[tagnum])
      endif
    endfor
    print,' '
    print,printstring
    print,' '

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showtags1,tagname,stack=n,help=help
 
; Display tags of the loaded single-channel spectrum, or else of
; stack1 spectrum #n if the stack keyword is used
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline

if n_params() eq 0 then tagname = 'all'

if n_params() gt 1 or keyword_set(help) then begin
  print,'showtags1[,tagname][,stack=n][,/help]'
  print,' '
  print,'If no CW tag name is given, all CW tags are displayed'
  return
endif else if size(tagname, /type) ne 7 then begin
  print,' '
  print,'showtags1,tagname[,stack=n][,/help]'
  print,'Make sure that tagname is a quoted string!'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in showtags1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showtags1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in showtags1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in showtags1: No single-channel spectrum is loaded'
  return
endif

; Get some elements of the loaded single-channel spectrum or else of
; stack1 spectrum #n

if n_elements(n) eq 0 then begin
  tags1 = (*loaded1).tags
  ntags = (*loaded1).ntags
endif else begin
  tags1 = (*stack1[n-1]).tags
  ntags = (*stack1[n-1]).ntags
endelse
tagnames = strlowcase(tag_names(tags1))

; Check that tagname is a valid tag name

tagname = strlowcase(tagname)
if (tagname ne 'all') then begin
  tagnum = where(tagnames eq tagname, count)
  if (count eq 0) then begin
    print,tagname,' is not a valid tag name'
    return
  endif else begin
    tagnum = tagnum[0]
  endelse
endif

; Display the tag(s)

if (tagname eq 'all') then begin

  ; print all tag values

  for k0 = 0L,ntags-1,iomaxline do begin
    lastOnLine = (k0 + iomaxline - 1) < (ntags - 1)
    printstring = 'TAG: '
    for k=k0,lastOnLine do begin
      printstring = printstring + string(tagnames[k],format=nameformat)
    endfor
    print,' '
    print,printstring
    printstring = '     '
    for k=k0,lastOnLine do begin
      printstring = printstring + string(tags1.(k),format=tagformats[k])
    endfor
    print,printstring
  endfor
  print,' '

endif else begin

  ; print one tag value

  printstring = tagname + ':' $
                + string(tags1.(tagnum),format=tagformats[tagnum])
  print,' '
  print,printstring
  print,' '

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro settag,tagname,newvalue,chan=chan,stack=n,silent=silent,help=help

; Set a tag value for the loaded pair
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline
common channelBlock, chanstrings, maxchan

if n_params() ne 2 or keyword_set(help) then begin
  print,' '
  print,'settag,tagname,newvalue[,chan=1 or 2][,stack=n][,/silent][,/help]'
  print,' '
  print,'Set a cw tag value for the loaded pair, or else for stack pair #n'
  print,'    if the stack keyword is set'
  print,' '
  print,'newvalue can be a scalar which will be assigned to both the OC and SC'
  print,'    channels (unless the chan keyword is used); otherwise it can be a'
  print,'    two-element array [OCvalue, SCvalue], even if the chan keyword is used.'
  print,' '
  return
endif else if size(tagname, /type) ne 7 then begin
  print,' '
  print,'settag,tagname,newvalue[,chan=1 or 2][,stack=n][,/silent][,/help]'
  print,'Make sure that tagname is a quoted string!'
  print,' '
  return
endif else if n_elements(newvalue) gt 2 then begin
  print,' '
  print,'settag,tagname,newvalue[,chan=1 or 2][,stack=n][,/silent][,/help]'
  print,'newvalue must be a scalar or a two-element [OCvalue, SCvalue] array'
  print,' '
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in settag: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in settag: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in settag: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in settag: No pair is loaded'
  return
endif

; Get some elements of the loaded pair or of stack pair #n

tags = (n_elements(n) gt 0) ? (*stack[n-1]).tags : (*loaded).tags
tagnames = strlowcase(tag_names(tags[0]))
npol = n_elements(tags)

; Check that tagname is a valid tag name

tagname = strlowcase(tagname)
tagnum = where(tagnames eq tagname, count)
if (count eq 0) then begin
  print,tagname,' is not a valid tag name'
  return
endif else begin
  tagnum = tagnum[0]
endelse

; Check which values to assign to which channel

if n_elements(newvalue) gt 1 then begin
  if n_elements(chanvalue) eq npol then chanvalue = newvalue else begin
    print, "ERROR: settag: number of values must be 1 or the same as number or pols", n_elements(chanvalue), npol
    return
  endelse
endif else begin
  chanvalue = replicate(newvalue, npol)
endelse

; Check which channels to change

  setpol = intarr(npol)
if (n_elements(chan) eq 0) then begin
  setpol = setpol + 1
endif else if (chan gt 0 && chan le npol) then begin
  setpol[chan-1] = 1
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
  return
endelse

; Change the tag

printstring = ''
for ch=1,npol do begin
  if setpol[ch-1] then begin
    oldvalue = tags[ch-1].(tagnum)
    if n_elements(n) eq 0 then begin
      (*loaded).tags[ch-1].(tagnum) = chanvalue[ch-1]
    endif else begin
      (*stack[n-1]).tags[ch-1].(tagnum) = chanvalue[ch-1]
    endelse
    printstring = printstring + 'Changing ' + tagname + '(' + chanstrings[ch-1]     $
                  + ') from ' + string(oldvalue,format=tagformats[tagnum])        $
                  + ' to ' + string(chanvalue[ch-1],format=tagformats[tagnum])    $
                  + '     '
  endif
endfor
if not keyword_set(silent) then begin
  print,' '
  print,printstring
  print,' '
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro settag1,tagname,newvalue,stack=n,silent=silent,help=help

; Set a tag value for the loaded single-channel spectrum
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline

if n_params() ne 2 or keyword_set(help) then begin
  print,' '
  print,'settag1,tagname,newvalue[,stack=n][,/silent][,/help]'
  print,' '
  print,'Set a cw tag value for the loaded single-channel spectrum,'
  print,'or else for stack1 spectrum #n if the stack keyword is set'
  print,' '
  return
endif else if size(tagname, /type) ne 7 then begin
  print,' '
  print,'settag1,tagname,newvalue[,stack=n][,/silent][,/help]'
  print,'Make sure that tagname is a quoted string!'
  print,' '
  return
endif else if n_elements(newvalue) ne 1 then begin
  print,' '
  print,'settag1,tagname,newvalue[,chan=1 or 2][,stack=n][,/silent][,/help]'
  print,'newvalue must be a scalar'
  print,' '
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in settag1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in settag1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in settag1: There are only ',nstack1, $
          ' spectra in the single-channel stack', format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in settag1: No single-channel spectrum is loaded'
  return
endif

; Get some elements of the loaded spectrum

tags1 = (n_elements(n) gt 0) ? (*stack1[n-1]).tags : (*loaded1).tags
tagnames = strlowcase(tag_names(tags1))

; Check that tagname is a valid tag name

tagname = strlowcase(tagname)
tagnum = where(tagnames eq tagname, count)
if (count eq 0) then begin
  print,tagname,' is not a valid tag name'
  return
endif else begin
  tagnum = tagnum[0]
endelse

; Change the tag

oldvalue = tags1.(tagnum)
if n_elements(n) eq 0 then begin
  (*loaded1).tags.(tagnum) = newvalue
endif else begin
  (*stack1[n-1]).tags.(tagnum) = newvalue
endelse
printstring = 'Changing ' + tagname + ' from '                       $
              + string(oldvalue,format=tagformats[tagnum]) + ' to '  $
              + string(newvalue,format=tagformats[tagnum])
if not keyword_set(silent) then begin
  print,' '
  print,printstring
  print,' '
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro siglim,mode,lim1,lim2,chan=chan,stack=st,silent=silent,help=help

; Set the signal limits (tags jsnr1, jsnr2) of the loaded pair
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

if n_params() eq 0 then mode = -1  ; just so it's defined

if keyword_set(help) || n_params() eq 0 || n_params() gt 3 || $
     mode lt 0 || mode gt 2 || (mode eq 1 && n_params() eq 1) || $
     (mode eq 2 && n_params() lt 3) then begin
  print,' '
  print,'siglim,mode[,lim1[,lim2]][,chan=1 or 2][,/stack][,/silent][,/help]'
  print,'set signal limits based on'
  print,'  predicted bandwidth                  (mode = 0)'
  print,'  frequencies lim1, lim2, or +/- lim1  (mode = 1)'
  print,'  frequency bins lim1, lim2 (0-based)  (mode = 2)'
  print,' '
  print,'/stack sets signal limits for every stack pair'
  print,'       rather than for the loaded pair'
  print,' '
  print,'/silent suppresses screen output'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Set signal limits for the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in siglim: No pair is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Set signal limits for every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in siglim: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; Determine which polarization channels were specified

; Assign the left and right signal limits according to the
; method specified via the mode parameter

if not keyword_set(silent) then print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack[n] : *loaded
  tags = useStruc.tags
  df = tags[0].dfreq
  posfr = tags[0].posfr
  bin0 = tags[0].xjcen
  ndata = useStruc.ndata
  npol = n_elements(tags)
  f = (keyword_set(st)) ? posfr*df*(findgen(ndata) - bin0) : useStruc.freq

  setpol = intarr(npol)
  if n_elements(chan) eq 0 then begin
    setpol = setpol + 1
    chanstring = ''
  endif else if chan ge 1 && chan le npol then begin
    setpol[chan-1] = 1
    chanstring = chanstrings[chan-1]
  endif else begin
    print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
    return
  endelse

  ; Assign the left and right signal limits according to the
  ; method specified via the mode parameter

  if mode eq 0 then begin

    ; Use predicted bandwidth (and be a bit generous)

    if keyword_set(st) then begin
      diameter = getextra('diameter', stack=(n+1))
      period = getextra('period', stack=(n+1))
      lambda = getextra('lambda', stack=(n+1))
    endif else begin
      diameter = getextra('diameter')
      period = getextra('period')
      lambda = getextra('lambda')
    endelse
    if isnull(diameter) then begin
      print,'WARNING in siglim: Assuming diameter = 1.0 km'
      diameter = 1.0
    endif
    if isnull(period) then begin
      print,'WARNING in siglim: Assuming period = 1.0 hr'
      period = 1.0
    endif
    if isnull(lambda) then begin
      print,'WARNING in siglim: Assuming lambda = 0.1259632 m'
      lambda = 0.1259632
    endif

    bw = 4*!pi*diameter*1000.0/(lambda*period*3600.0)
    bin1 = bin0 - round(0.55*bw/df)
    bin2 = bin0 + round(0.55*bw/df)

  endif else if mode eq 1 then begin

    ; Specify two frequencies or a symmetric +/- interval

    if n_params() eq 3 then begin
      if posfr eq 1 then begin
        bin1 = bin0 + round(lim1/df)
        bin2 = bin0 + round(lim2/df)
      endif else begin
        bin1 = bin0 - round(lim2/df)
        bin2 = bin0 - round(lim1/df)
      endelse
    endif else begin
      bin1 = bin0 - round(lim1/df)
      bin2 = bin0 + round(lim1/df)
    endelse
  endif else begin

    ; Directly specify frequency bins

    bin1 = round(1.0*lim1)
    bin2 = round(1.0*lim2)
  endelse

  ; Make sure that the limits fall within the spectrum

  if bin1 lt 0 then begin
    bin1 = 0L
    if not keyword_set(silent) then print,'Resetting left signal limit to frequency bin 0'
  endif
  if bin2 ge ndata then begin
    bin2 = ndata - 1L
    if not keyword_set(silent) then print,'Resetting right signal limit to frequency bin ',bin2
  endif

  ; Display the results and set the tags

  if keyword_set(st) then begin
    startstring = 'Stack #' + string(n+1, format='(i0)') + ': Setting '
  endif else begin
    startstring = 'Setting '
  endelse
  if not keyword_set(silent) then begin
    print,startstring,chanstring,'signal limits to freq. bins ',  $
          bin1,'-',bin2,'  (',f[bin1]-0.5*posfr*df,' Hz to ',f[bin2]+0.5*posfr*df,' Hz)',  $
          format='(3a, i0, a, i0, a, f8.2, a, f8.2, a)'
    if (mode le 1 && bin0 ne ndata/2) then $
        print,'(Full range is bins 0-',ndata-1,', 0 Hz is bin ',bin0,')', $
              format='(a,i0,a,i0,a)'
  endif

  for ch=1,npol do begin
    if setpol[ch-1] then begin

      if keyword_set(st) then begin
        (*stack[n]).tags[ch-1].jsnr1 = bin1
        (*stack[n]).tags[ch-1].jsnr2 = bin2
      endif else begin
        (*loaded).tags[ch-1].jsnr1 = bin1
        (*loaded).tags[ch-1].jsnr2 = bin2
      endelse
    endif
  endfor

endfor

if not keyword_set(silent) then print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro siglim1,mode,lim1,lim2,chanstring=chanstring,stack=st,silent=silent,help=help

; Set the signal limits (tags jsnr1, jsnr2) of the loaded single-channel spectrum

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

if n_params() eq 0 then mode = -1 ; just so it's defined

if keyword_set(help) || n_params() eq 0 || n_params() gt 3 || $
     mode lt 0 || mode gt 2 || (mode eq 1 && n_params() eq 1) || $
     (mode eq 2 && n_params() lt 3) then begin
  print,' '
  print,'siglim1,mode[,lim1[,lim2]][,chanstring=chanstring][,/stack][,/silent][,/help]'
  print,'set signal limits based on'
  print,'  predicted bandwidth                  (mode = 0)'
  print,'  frequencies lim1, lim2, or +/- lim1  (mode = 1)'
  print,'  frequency bins lim1, lim2 (0-based)  (mode = 2)'
  print,' '
  print,'/stack sets signal limits for every stack1 spectrum'
  print,'       rather than for the loaded single-channel spectrum'
  print,' '
  print,'/silent suppresses screen output'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Set signal limits for the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in siglim1: No single-channel spectrum is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Set signal limits for every spectrum in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in siglim1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
endelse

if n_elements(chanstring) eq 0 then begin
  chanstring = ''
endif else begin
  chanstring = chanstring + ' '
endelse
if not keyword_set(silent) then print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack1[n] : *loaded1
  tags1 = useStruc.tags
  df = tags1.dfreq
  posfr = tags1.posfr
  bin0 = tags1.xjcen
  ndata = useStruc.ndata
  f = (keyword_set(st)) ? posfr*df*(findgen(ndata) - bin0) : useStruc.freq

  ; Assign the left and right signal limits according to the
  ; method specified via the mode parameter

  if mode eq 0 then begin

    ; Use predicted bandwidth (and be a bit generous)

    if keyword_set(st) then begin
      diameter = getextra1('diameter', stack=(n+1))
      period = getextra1('period', stack=(n+1))
      lambda = getextra1('lambda', stack=(n+1))
    endif else begin
      diameter = getextra1('diameter')
      period = getextra1('period')
      lambda = getextra1('lambda')
    endelse
    if isnull(diameter) then begin
      print,'WARNING in siglim1: Assuming diameter = 1.0 km'
      diameter = 1.0
    endif
    if isnull(period) then begin
      print,'WARNING in siglim1: Assuming period = 1.0 hr'
      period = 1.0
    endif
    if isnull(lambda) then begin
      print,'WARNING in siglim1: Assuming lambda = 0.1259632 m'
      lambda = 0.1259632
    endif

    bw = 4*!pi*diameter*1000.0/(lambda*period*3600.0)
    bin1 = bin0 - round(0.55*bw/df)
    bin2 = bin0 + round(0.55*bw/df)

  endif else if mode eq 1 then begin

    ; Specify two frequencies or a symmetric +/- interval

    if n_params() eq 3 then begin
      if posfr eq 1 then begin
        bin1 = bin0 + round(lim1/df)
        bin2 = bin0 + round(lim2/df)
      endif else begin
        bin1 = bin0 - round(lim2/df)
        bin2 = bin0 - round(lim1/df)
      endelse
    endif else begin
      bin1 = bin0 - round(lim1/df)
      bin2 = bin0 + round(lim1/df)
    endelse
  endif else begin

    ; Directly specify frequency bins

    bin1 = round(1.0*lim1)
    bin2 = round(1.0*lim2)
  endelse

  ; Make sure that the limits fall within the spectrum

  if bin1 lt 0 then begin
    bin1 = 0L
    if not keyword_set(silent) then print,'Resetting left ',chanstring,'signal limit to frequency bin 0'
  endif
  if bin2 ge ndata then begin
    bin2 = ndata - 1L
    if not keyword_set(silent) then print,'Resetting right ',chanstring,'signal limit to frequency bin ',bin2
  endif

  ; Display the results and set the tags

  if keyword_set(st) then begin
    startstring = 'Stack1 #' + string(n+1, format='(i0)') + ': Setting '
  endif else begin
    startstring = 'Setting '
  endelse
  if not keyword_set(silent) then begin
    print,startstring,chanstring,'signal limits to frequency bins ',  $
          bin1,'-',bin2,'  (',f[bin1]-0.5*posfr*df,' Hz to ',f[bin2]+0.5*posfr*df,' Hz)',   $
          format='(3a, i0, a, i0, a, f8.2, a, f8.2, a)'
    if (mode le 1 && bin0 ne ndata/2) then $
        print,'(Full range is bins 0-',ndata-1,', 0 Hz is bin ',bin0,')', $
              format='(a,i0,a,i0,a)'
  endif

  if keyword_set(st) then begin
    (*stack1[n]).tags.jsnr1 = bin1
    (*stack1[n]).tags.jsnr2 = bin2
  endif else begin
    (*loaded1).tags.jsnr1 = bin1
    (*loaded1).tags.jsnr2 = bin2
  endelse

endfor

if not keyword_set(silent) then print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro siglimi,mode,doprange=doprange,delrange=delrange,stack=st,silent=silent,help=help

; Set the signal limits (extra tags doplim1, doplim2, dellim1, dellim2)
; of the loaded image

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then mode = -1 ; just so it's defined

if keyword_set(help) or n_params() ne 1 or mode lt 0 or mode gt 2                    $
        or (mode eq 0 && (n_elements(delrange) gt 0 or n_elements(doprange) gt 0))  $
        or (mode ne 0 && n_elements(delrange) eq 0 && n_elements(doprange) eq 0)   $
        then begin
  print,' '
  print,'siglimi,mode[,doprange=doprange][,delrange=delrange][,/stack][,/silent][,/help]'
  print,'set signal limits based on'
  print,'  mode = 0: full delay range, predicted Doppler bandwidth'
  print,'            (do not use delrange or doprange in this mode)'
  print,'  mode = 1: delay range in usec and/or Doppler range in Hz'
  print,'            (doprange=100 is the same as doprange=[-100,100])'
  print,'  mode = 2: delay bins and/or Doppler bins (0-based)'
  print,' '
  print,'/stack sets signal limits for every stacki image'
  print,'       rather than for the loaded image'
  print,' '
  print,'/silent suppresses screen output'
  print,' '
  return
endif

; Get delay and/or Doppler limits from the keywords

if n_elements(doprange) gt 0 then begin
  if n_elements(doprange) eq 2 then begin
    doplim1 = doprange[0]
    doplim2 = doprange[1]
    if doplim2 lt doplim1 then begin
      print,'ERROR in siglimi: Must have doprange[1] >= doprange[0]'
      return
    endif
  endif else if n_elements(doprange) eq 1 then begin
    doplim1 = -doprange[0]
    doplim2 = doprange[0]
    if doplim2 lt 0 then begin
      print,'ERROR in siglimi: Must have scalar doprange >= 0'
      return
    endif
  endif else begin
    print,'ERROR in siglimi: doprange must be a scalar or a 2-element vector'
    return
  endelse
endif

if n_elements(delrange) gt 0 then begin
  if n_elements(delrange) eq 2 then begin
    dellim1 = delrange[0]
    dellim2 = delrange[1]
    if dellim2 lt dellim1 then begin
      print,'ERROR in siglimi: Must have delrange[1] >= delrange[0]'
      return
    endif
  endif else begin
    print,'ERROR in siglimi: delrange must be a 2-element vector'
    return
  endelse
endif

if not keyword_set(st) then begin

  ; Set signal limits for the loaded image

  if (*loadedi).width le 2 && (*loadedi).height le 2 then begin
    print,'ERROR in siglimi: No image is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Set signal limits for every image in the stack

  if nstacki eq 0 then begin
    print,'ERROR in siglimi: The image stack is empty'
    return
  endif
  nuse = nstacki
endelse

if not keyword_set(silent) then print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this image

  if keyword_set(st) then begin
    width = (*stacki[n]).width
    height = (*stacki[n]).height
    dfreq = getextrai('fres',stack=(n+1))
    eph_col = getextrai('eph_col',stack=(n+1))
    eph_row = getextrai('eph_row',stack=(n+1))
    delayunit = getextrai('delayunit',stack=(n+1))
    spb = getextrai('samples_per_baud',stack=(n+1))
    rows_per_baud = getextrai('rows_per_baud',stack=(n+1))
    baud = getextrai('baudlen',stack=(n+1))
    diameter = getextrai('diameter',stack=(n+1))
    period = getextrai('period',stack=(n+1))
    lambda = getextrai('lambda',stack=(n+1))
  endif else begin
    width = (*loadedi).width
    height = (*loadedi).height
    dfreq = getextrai('fres')
    eph_col = getextrai('eph_col')
    eph_row = getextrai('eph_row')
    delayunit = getextrai('delayunit')
    spb = getextrai('samples_per_baud')
    rows_per_baud = getextrai('rows_per_baud')
    baud = getextrai('baudlen')
    diameter = getextrai('diameter')
    period = getextrai('period')
    lambda = getextrai('lambda')
  endelse
  if isnull(delayunit) && notnull(baud) then begin
    if isnull(spb) then spb = 1L
    if isnull(rows_per_baud) then rows_per_baud = spb
    delayunit = baud/rows_per_baud
  endif
  if notnull(eph_col) && notnull(eph_row) && notnull(dfreq) && notnull(delayunit) then begin
    f = dfreq*(findgen(width) - eph_col)       ; Doppler (Hz)
    d = delayunit*(findgen(height) - eph_row)  ; delay (usec)
  endif else begin
    print,"ERROR in siglimi: Need the 'eph_col' and 'eph_row' and 'fres' tags,"
    print,"                  plus either the 'delayunit' or 'baudlen' tag"
    return
  endelse

  ; Assign the left and right signal limits according to the
  ; method specified via the mode parameter

  change_delay = mode eq 0 or n_elements(delrange) gt 0
  change_doppler = mode eq 0 or n_elements(doprange) gt 0

  if mode eq 0 then begin

    ; Use full delay range;
    ; use predicted Doppler bandwidth (and be a bit generous)

    delbin1 = 0L
    delbin2 = height - 1

    if isnull(diameter) then begin
      print,'WARNING in siglimi: Assuming diameter = 1.0 km'
      diameter = 1.0
    endif
    if isnull(period) then begin
      print,'WARNING in siglimi: Assuming period = 1.0 hr'
      period = 1.0
    endif
    if isnull(lambda) then begin
      print,'WARNING in siglimi: Assuming lambda = 0.1259632 m'
      lambda = 0.1259632
    endif

    bw = 4*!pi*diameter*1000.0/(lambda*period*3600.0)
    dopbin1 = round(eph_col - 0.55*bw/dfreq)
    dopbin2 = round(eph_col + 0.55*bw/dfreq)

  endif else if mode eq 1 then begin

    ; Specify delay range (in usec) and/or Doppler range (in Hz)

    if change_delay then begin
      delbin1 = round(eph_row + dellim1/delayunit)
      delbin2 = round(eph_row + dellim2/delayunit)
    endif

    if change_doppler then begin
      dopbin1 = round(eph_col + doplim1/dfreq)
      dopbin2 = round(eph_col + doplim2/dfreq)
    endif

  endif else begin

    ; Directly specify delay bins and/or Doppler bins

    if change_delay then begin
      delbin1 = round(1.0*dellim1)
      delbin2 = round(1.0*dellim2)
    endif

    if change_doppler then begin
      dopbin1 = round(1.0*doplim1)
      dopbin2 = round(1.0*doplim2)
    endif

  endelse

  ; Make sure that the limits fall within the spectrum

  if change_delay then begin
    if delbin1 lt 0 then begin
      delbin1 = 0L
      if not keyword_set(silent) then print,'Resetting low-delay signal limit to bin 0'
    endif
    if delbin2 ge height then begin
      delbin2 = height - 1
      if not keyword_set(silent) then $
          print,'Resetting high-delay signal limit to bin ',delbin2,format='(a,i0)'
    endif
  endif

  if change_doppler then begin
    if dopbin1 lt 0 then begin
      dopbin1 = 0L
      if not keyword_set(silent) then print,'Resetting low-Doppler signal limit to bin 0'
    endif
    if dopbin2 ge width then begin
      dopbin2 = width - 1
      if not keyword_set(silent) then $
          print,'Resetting high-Doppler signal limit to bin ',dopbin2,format='(a,i0)'
    endif
  endif

  ; Display the results

  if not keyword_set(silent) then begin
    if keyword_set(st) then begin
      printstring1 = 'Stacki #' + string(n+1, format='(i0)') $
                               + ': Setting '
    endif else begin
      printstring1 = 'Setting '
    endelse
    nskip1 = strlen(printstring1)
    if change_doppler then begin
      printstring1 = printstring1 + 'Doppler limits to bins '  $
                     + string(dopbin1,format='(i0)')  + '-'    $
                     + string(dopbin2,format='(i0)')
      if change_delay then begin
        formatstring = '(' + string(nskip1,format='(i0)') + 'x,a)'
        printstring2 = string('delay   limits to bins ',format=formatstring)  $
                       + string(delbin1,format='(i0)')  + '-'                 $
                       + string(delbin2,format='(i0)')
        nskip2 = strlen(printstring1) - strlen(printstring2)
        if nskip2 gt 0 then begin
          for i=0L,nskip2-1 do printstring2 = printstring2 + ' '
        endif else if nskip2 lt 0 then begin
          for i=0L,-nskip2-1 do printstring1 = printstring1 + ' '
        endif
      endif
      printstring1 = printstring1 + '  ('                                        $
                     + string(f[dopbin1]-0.5*dfreq,format='(f8.2)') + ' Hz to '  $
                     + string(f[dopbin2]+0.5*dfreq,format='(f8.2)') + ' Hz)'
      print,printstring1
      if change_delay then begin
        printstring2 = printstring2 + '  ('                                            $
                       + string(d[delbin1]-0.5*delayunit,format='(f8.2)') + ' us to '  $
                       + string(d[delbin2]+0.5*delayunit,format='(f8.2)') + ' us)'
        print,printstring2
      endif
    endif else begin
      printstring1 = printstring1 + 'delay limits to bins '                          $
                     + string(delbin1,format='(i0)')  + '-'                          $
                     + string(delbin2,format='(i0)') + '  ('                         $
                     + string(d[delbin1]-0.5*delayunit,format='(f8.2)') + ' us to '  $
                     + string(d[delbin2]+0.5*delayunit,format='(f8.2)') + ' us)'
      print,printstring1
    endelse
  endif

  ; Write the extra tags

  dop1comment = 'min Doppler bin containing signal, 0-based'
  dop2comment = 'max Doppler bin containing signal, 0-based'
  del1comment = 'min delay bin containing signal, 0-based'
  del2comment = 'max delay bin containing signal, 0-based'

  if keyword_set(st) then begin
    if change_doppler then begin
      setextrai,'i','dopbin1',dopbin1,comment=dop1comment,stack=(n+1)
      setextrai,'i','dopbin2',dopbin2,comment=dop2comment,stack=(n+1)
    endif
    if change_delay then begin
      setextrai,'i','delbin1',delbin1,comment=del1comment,stack=(n+1)
      setextrai,'i','delbin2',delbin2,comment=del2comment,stack=(n+1)
    endif
  endif else begin
    if change_doppler then begin
      setextrai,'i','dopbin1',dopbin1,comment=dop1comment
      setextrai,'i','dopbin2',dopbin2,comment=dop2comment
    endif
    if change_delay then begin
      setextrai,'i','delbin1',delbin1,comment=del1comment
      setextrai,'i','delbin2',delbin2,comment=del2comment
    endif
  endelse

endfor

if not keyword_set(silent) then print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro vignette,mode,lim1,lim2,stack=st,help=help

; Vignette the loaded pair
;
; If keyword st is set, vignette every pair in the stack instead
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then mode = -1 ; just so it's defined

if keyword_set(help) or n_params() eq 0 or n_params() gt 3 || $
     mode lt 0 or mode gt 2 || $
     (mode eq 0 && (n_params() ne 1 || not keyword_set(st))) || $
     (mode eq 1 && n_params() eq 1) || $
     (mode eq 2 && n_params() lt 3) then begin
  print,' '
  print,'vignette,mode[,lim1[,lim2]][,/stack][,/help]'
  print,' '
  print,'vignette pair based on'
  print,'  max freq range common to all stack pairs (mode = 0)', $
        ' (must set /stack)',format='(2a)'
  print,'  frequencies lim1, lim2, or +/- lim1      (mode = 1)'
  print,'  frequency bins lim1, lim2 (0-based)      (mode = 2)'
  print,' '
  print,'/stack vignettes every stack pair instead of the loaded pair'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Vignette the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in vignette: No pair is loaded'
    return
  endif
  nuse = 1
endif else begin

  ; Vignette every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in vignette: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; Set the min and max frequencies (or frequency bins)

if n_params() eq 1 then begin

  ; mode = 0: Find the maximum frequency range common to all stack pairs
  
  lim1_use = -9.99e20
  lim2_use = 9.99e20
  for n=0L,nuse-1 do begin
    tags = (*stack[n]).tags
    dfreq = double(tags[0].dfreq) ; need precision for large ndata
    posfr = tags[0].posfr
    xjcen = tags[0].xjcen
    ndata = (*stack[n]).ndata
    if posfr eq 1 then begin
      fmin = -dfreq*(xjcen + 0.5)
      fmax = dfreq*(ndata - xjcen - 0.5)
    endif else begin
      fmin = -dfreq*(ndata - xjcen - 0.5)
      fmax = dfreq*(xjcen + 0.5)
    endelse
    lim1_use = lim1_use > fmin
    lim2_use = lim2_use < fmax
  endfor
  if lim1_use ge lim2_use then begin
    print,'ERROR in vignette: The pair stack has no common frequency range'
    return
  endif
endif else if n_params() eq 2 then begin
  lim1_use = -lim1
  lim2_use = lim1
endif else begin
  lim1_use = lim1
  lim2_use = lim2
endelse

; Now process the pair(s)

for n=0L,nuse-1 do begin

  ; Get a few elements of this pair

  useStruc = (keyword_set(st)) ? *stack[n] : *loaded
  pair = useStruc.spec
  tags = useStruc.tags
  dfreq = double(tags[0].dfreq)
  posfr = tags[0].posfr
  xjcen = tags[0].xjcen
  jsnr1 = tags.jsnr1
  jsnr2 = tags.jsnr2
  ndata = useStruc.ndata
  npol = n_elements(pair[*,0])
  f = (keyword_set(st)) ? posfr*dfreq*(findgen(ndata) - xjcen) : useStruc.freq
  
  ; Assign the left and right limits according to the
  ; method specified via the mode parameter;
  ; be careful about roundoff error influencing the limits

  if mode le 1 then begin

    ; Convert frequency limits to frequency bins

    if posfr eq 1 then begin
      bin1 = xjcen + round(lim1_use/dfreq + 1e-7)
      bin2 = xjcen + round(lim2_use/dfreq - 1e-7)
    endif else begin
      bin1 = xjcen - round(lim2_use/dfreq - 1e-7)
      bin2 = xjcen - round(lim1_use/dfreq + 1e-7)
    endelse
  endif else begin

    ; Directly specify frequency bins

    bin1 = round(1.0*lim1_use + 1e-7)
    bin2 = round(1.0*lim2_use - 1e-7)
  endelse

  ; Make sure that the limits fall within the spectrum
  ; and that more than two frequency bins remain after vignetting

  if bin1 lt 0 or bin2 ge ndata then begin
    print,'ERROR on pair #',n+1,': Requested frequency bins [', $
          bin1,'-',bin2,'] extend beyond spectrum [0-',ndata-1,']', $
          format='(4(a,i0),a)'
    print,'       Pair left unchanged'
  endif else if bin2 le bin1+1 then begin
    print,'ERROR on pair #',n+1,': Requested frequency bins [', $
          bin1,'-',bin2,'] would leave too few points in vignetted spectrum', $
          format='(3(a,i0),a)'
    print,'       Pair left unchanged'
  endif else begin
    ndata = bin2 - bin1 + 1
    tags.xjcen = xjcen - bin1
    newjsnr1 = jsnr1 - bin1
    newjsnr2 = jsnr2 - bin1
    for ch=1,npol do begin
      if newjsnr1[ch-1] lt 0 && newjsnr2[ch-1] lt 0 then begin
        print,'WARNING on channel ',ch,' of pair #',n+1, $
              ': Signal range reset to leftmost frequency bin', $
              format='(a,i0,a)'
        newjsnr1[ch-1] = 0L
        newjsnr2[ch-1] = 0L
      endif else if newjsnr1[ch-1] ge ndata && newjsnr2[ch-1] ge ndata then begin
        print,'WARNING on channel ',ch,' of pair #',n+1, $
              ': Signal range reset to rightmost frequency bin', $
              format='(a,i0,a,i0,a)'
        newjsnr1[ch-1] = ndata - 1
        newjsnr2[ch-1] = ndata - 1
      endif else begin
        newjsnr1[ch-1] = newjsnr1[ch-1] > 0L
        newjsnr2[ch-1] = newjsnr2[ch-1] < (ndata - 1)
      endelse
    endfor
    tags.jsnr1 = newjsnr1
    tags.jsnr2 = newjsnr2
    print,'Vignetting pair #',n+1,' to ',ndata,' points, frequency bins ',  $
          bin1,'-',bin2,'  (',f[bin1]-0.5*posfr*dfreq,' Hz to ',f[bin2]+0.5*posfr*dfreq,' Hz)', $
          format='(a,i0,a,i0,a,i0,a,i0,a,f9.2,a,f9.2,a)'
    if keyword_set(st) then begin
      newStruc = {group:useStruc.group, spec:pair[*,bin1:bin2], tags:tags, $
                  extratags:useStruc.extratags, ndata:ndata,  $
                  ntags:useStruc.ntags, nextra:useStruc.nextra, tname:useStruc.tname}
    endif else begin
      newStruc = {freq:f[bin1:bin2], spec:pair[*,bin1:bin2], tags:tags, $
                  extratags:useStruc.extratags, ndata:ndata, ntags:useStruc.ntags, $
                  nextra:useStruc.nextra, tname:useStruc.tname}
    endelse

    if keyword_set(st) then *stack[n] = newStruc else *loaded = newStruc

  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro vignette1,mode,lim1,lim2,stack=st,help=help

; Vignette the loaded single-channel spectrum
;
; If keyword st is set, vignette every spectrum in the
; single-channel stack instead
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then mode = -1 ; just so it's defined

if keyword_set(help) or n_params() eq 0 or n_params() gt 3 or $
     mode lt 0 or mode gt 2 or $
     (mode eq 0 and (n_params() ne 1 or not keyword_set(st))) or $
     (mode eq 1 and n_params() eq 1) or $
     (mode eq 2 and n_params() lt 3) then begin
  print,' '
  print,'vignette1,mode[,lim1[,lim2]][,/stack][,/help]'
  print,' '
  print,'vignette single-channel spectrum based on'
  print,'  max freq range common to all stack1 spectra (mode = 0)', $
        ' (must set /stack)',format='(2a)'
  print,'  frequencies lim1, lim2, or +/- lim1         (mode = 1)'
  print,'  frequency bins lim1, lim2 (0-based)         (mode = 2)'
  print,' '
  print,'/stack vignettes every stack1 spectrum instead of the ', $
           'loaded single-channel spectrum',format='(2a)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Vignette the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in vignette1: No single-channel spectrum is loaded'
    return
  endif
  nuse = 1
endif else begin

  ; Vignette every spectrum in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in vignette1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
endelse

; Set the min and max frequencies (or frequency bins)

if n_params() eq 1 then begin

  ; mode = 0: Find the maximum frequency range common to all stack1 spectra
  
  lim1_use = -9.99e20
  lim2_use = 9.99e20
  for n=0L,nuse-1 do begin
    tags1 = (*stack1[n]).tags
    dfreq = double(tags1.dfreq) ; need precision for large ndata
    posfr = tags1.posfr
    xjcen = tags1.xjcen
    ndata = (*stack1[n]).ndata
    if posfr eq 1 then begin
      fmin = -dfreq*(xjcen + 0.5)
      fmax = dfreq*(ndata - xjcen - 0.5)
    endif else begin
      fmin = -dfreq*(ndata - xjcen - 0.5)
      fmax = dfreq*(xjcen + 0.5)
    endelse
    lim1_use = lim1_use > fmin
    lim2_use = lim2_use < fmax
  endfor
  if lim1_use ge lim2_use then begin
    print,'ERROR in vignette1: The single-channel stack has ', $
          'no common frequency range',format='(2a)'
    return
  endif
endif else if n_params() eq 2 then begin
  lim1_use = -lim1
  lim2_use = lim1
endif else begin
  lim1_use = lim1
  lim2_use = lim2
endelse

; Now process the spectrum (or spectra)

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack1[n] : *loaded1
  spec = useStruc.spec
  tags1 = useStruc.tags
  dfreq = tags1.dfreq
  posfr = tags1.posfr
  xjcen = tags1.xjcen
  jsnr1 = tags1.jsnr1
  jsnr2 = tags1.jsnr2
  ndata = useStruc.ndata
  f = (keyword_set(st)) ? posfr*dfreq*(findgen(ndata) - xjcen) : useStruc.freq
  
  ; Assign the left and right limits according to the
  ; method specified via the mode parameter;
  ; be careful about roundoff error influencing the limits

  if mode le 1 then begin

    ; Convert frequency limits to frequency bins

    if posfr eq 1 then begin
      bin1 = xjcen + round(lim1_use/dfreq + 1e-7)
      bin2 = xjcen + round(lim2_use/dfreq - 1e-7)
    endif else begin
      bin1 = xjcen - round(lim2_use/dfreq - 1e-7)
      bin2 = xjcen - round(lim1_use/dfreq + 1e-7)
    endelse
  endif else begin

    ; Directly specify frequency bins

    bin1 = round(1.0*lim1_use + 1e-7)
    bin2 = round(1.0*lim2_use - 1e-7)
  endelse

  ; Make sure that the limits fall within the spectrum
  ; and that more than two frequency bins remain after vignetting

  if bin1 lt 0 or bin2 ge ndata then begin
    print,'ERROR on spectrum #',n+1,': Requested frequency bins [', $
          bin1,'-',bin2,'] extend beyond spectrum [0-',ndata-1,']', $
          format='(4(a,i0),a)'
    print,'       Spectrum left unchanged'
  endif else if bin2 le bin1+1 then begin
    print,'ERROR on spectrum #',n+1,': Requested frequency bins [', $
          bin1,'-',bin2,'] would leave too few points in vignetted spectrum', $
          format='(3(a,i0),a)'
    print,'       Spectrum left unchanged'
  endif else begin
    ndata = bin2 - bin1 + 1
    tags1.xjcen = xjcen - bin1
    newjsnr1 = jsnr1 - bin1
    newjsnr2 = jsnr2 - bin1
    if newjsnr1 lt 0 and newjsnr2 lt 0 then begin
      print,'WARNING on spectrum #',n+1,': Signal range reset to leftmost frequency bin', $
            format='(a,i0,a)'
      tags1.jsnr1 = 0L
      tags1.jsnr2 = 0L
    endif else if newjsnr1 ge ndata and newjsnr2 ge ndata then begin
      print,'WARNING on spectrum #',n+1,': Signal range reset to rightmost frequency bin', $
            format='(a,i0,a)'
      tags1.jsnr1 = ndata - 1
      tags1.jsnr2 = ndata - 1
    endif else begin
      tags1.jsnr1 = newjsnr1 > 0L
      tags1.jsnr2 = newjsnr2 < (ndata - 1)
    endelse
    print,'Vignetting spectrum #',n+1,' to ',ndata,' points, frequency bins ',  $
          bin1,'-',bin2,'  (',f[bin1]-0.5*posfr*dfreq,' Hz to ',f[bin2]+0.5*posfr*dfreq,' Hz)',     $
          format='(a,i0,a,i0,a,i0,a,i0,a,f9.2,a,f9.2,a)'
    if keyword_set(st) then begin
      newStruc = {group:useStruc.group, spec:spec[bin1:bin2], tags:tags1, $
                  extratags:useStruc.extratags, ndata:ndata, ntags:useStruc.ntags, $
                  nextra:useStruc.nextra, tname:useStruc.tname, pol:useStruc.pol}
    endif else begin
      newStruc = {freq:f[bin1:bin2], spec:spec[bin1:bin2], tags:tags1, $
                  extratags:useStruc.extratags, ndata:ndata, ntags:useStruc.ntags, $
                  nextra:useStruc.nextra, tname:useStruc.tname, pol:useStruc.pol}
    endelse

    if keyword_set(st) then *stack1[n] = newStruc else *loaded1 = newStruc

  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro vignettei,mode,doprange=doprange,delrange=delrange,stack=st,help=help

; Vignette the loaded image
;
; If keyword st is set, vignette every image in the image stack instead

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then mode = -1 ; just so it's defined

if keyword_set(help) or n_params() ne 1 or mode lt 0 or mode gt 2                   $
        or (mode eq 0 and (n_elements(delrange) gt 0 or n_elements(doprange) gt 0   $
                                                 or not keyword_set(st)  ))         $
        or (mode ne 0 and n_elements(delrange) eq 0 and n_elements(doprange) eq 0)  $
        then begin
  print,' '
  print,'vignettei,mode[,doprange=doprange][,delrange=delrange][,/stack][,/help]'
  print,' '
  print,'vignette image based on'
  print,'  mode = 0: max delay-Doppler range common to all stacki images'
  print,'            (must set /stack and not use delrange or doprange in this mode)'
  print,'  mode = 1: delay range in usec and/or Doppler range in Hz'
  print,'            (doprange=100 is the same as doprange=[-100,100])'
  print,'  mode = 2: delay bins and/or Doppler bins (0-based)'
  print,' '
  print,'/stack vignettes every stacki image instead of the loaded image'
  print,' '
  return
endif

; Decide which images we're vignetting, check that all necessary extra tags
; are present, and store some information for later use

if not keyword_set(st) then begin

  ; Vignette the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,'ERROR in vignettei: No image is loaded'
    return
  endif
  nuse = 1L

  width = (*loadedi).width
  height = (*loadedi).height
  dfreq = getextrai('fres')
  eph_col = getextrai('eph_col')
  eph_row = getextrai('eph_row')
  dc_col = getextrai('dc_col')
  delayunit = getextrai('delayunit')
  spb = getextrai('samples_per_baud')
  rows_per_baud = getextrai('rows_per_baud')
  baud = getextrai('baudlen')
  if isnull(delayunit) and notnull(baud) then begin
    if isnull(spb) then spb = 1L
    if isnull(rows_per_baud) then rows_per_baud = spb
    delayunit = baud/rows_per_baud
  endif
  if isnull(eph_col) or isnull(eph_row) or isnull(dfreq) or isnull(delayunit) then begin
    print,"ERROR in vignettei: The loaded image", $
          " needs the 'eph_col' and 'eph_row' and 'fres' tags,",format='(2a)'
    print,"                    plus either the 'delayunit' or 'baudlen' tag"
    return
  endif

endif else begin

  ; Vignette every image in the image stack

  if nstacki eq 0 then begin
    print,'ERROR in vignettei: The image stack is empty'
    return
  endif
  nuse = nstacki
  width = fltarr(nuse)
  height = fltarr(nuse)
  dfreq = dblarr(nuse) ; need precision for large width
  delayunit = fltarr(nuse)
  eph_col = fltarr(nuse)
  eph_row = fltarr(nuse)
  dc_col = fltarr(nuse)

  for n=0L,nuse-1 do begin
    width[n] = (*stacki[n]).width
    height[n] = (*stacki[n]).height
    dfreq[n] = getextrai('fres',stack=(n+1))
    eph_col[n] = getextrai('eph_col',stack=(n+1))
    eph_row[n] = getextrai('eph_row',stack=(n+1))
    dc_col[n] = getextrai('dc_col',stack=(n+1))
    delayunit[n] = getextrai('delayunit',stack=(n+1))
    spb = getextrai('samples_per_baud',stack=(n+1))
    rows_per_baud = getextrai('rows_per_baud',stack=(n+1))
    baud = getextrai('baudlen',stack=(n+1))
    if isnull(delayunit[n]) and notnull(baud) then begin
      if isnull(spb) then spb = 1L
      if isnull(rows_per_baud) then rows_per_baud = spb
      delayunit[n] = baud/rows_per_baud
    endif
    if isnull(eph_col) or isnull(eph_row) or isnull(dfreq) or isnull(delayunit) then begin
      print,'ERROR in vignettei: stacki image #',n+1, $
            " needs the 'eph_col' and 'eph_row' and 'fres' tags,",format='(a,i0,a)'
      print,"                    plus either the 'delayunit' or 'baudlen' tag"
      return
    endif
  endfor

endelse

; Set the min and max delay-Doppler values (or image bins)

if mode eq 0 then begin

  ; Find the maximum delay-Doppler range common to all stacki images
  
  doplim1 = -9.99e20
  doplim2 = 9.99e20
  dellim1 = -9.99e20
  dellim2 = 9.99e20
  for n=0L,nuse-1 do begin
    fmin = -dfreq[n]*(eph_col[n] + 0.5)                 ; Doppler (Hz)
    fmax = dfreq[n]*(width[n] - eph_col[n] - 0.5)
    dmin = -delayunit[n]*(eph_row[n] + 0.5)             ; delay (usec)
    dmax = delayunit[n]*(height[n] - eph_row[n] - 0.5)
    doplim1 = doplim1 > fmin
    doplim2 = doplim2 < fmax
    dellim1 = dellim1 > dmin
    dellim2 = dellim2 < dmax
  endfor
  if doplim1 ge doplim2 then begin
    print,'ERROR in vignettei: The image stack has no common Doppler range'
    return
  endif else if dellim1 ge dellim2 then begin
    print,'ERROR in vignettei: The image stack has no common delay range'
    return
  endif

endif else begin

  ; Get delay and/or Doppler limits from the keywords

  if n_elements(doprange) gt 0 then begin
    if n_elements(doprange) eq 2 then begin
      doplim1 = doprange[0]
      doplim2 = doprange[1]
      if doplim2 lt doplim1 then begin
        print,'ERROR in vignettei: Must have doprange[1] >= doprange[0]'
        return
      endif
    endif else if n_elements(doprange) eq 1 then begin
      doplim1 = -doprange[0]
      doplim2 = doprange[0]
      if doplim2 lt 0 then begin
        print,'ERROR in vignettei: Must have scalar doprange >= 0'
        return
      endif
    endif else begin
      print,'ERROR in vignettei: doprange must be a scalar or a 2-element vector'
      return
    endelse
  endif

  if n_elements(delrange) gt 0 then begin
    if n_elements(delrange) eq 2 then begin
      dellim1 = delrange[0]
      dellim2 = delrange[1]
      if dellim2 lt dellim1 then begin
        print,'ERROR in vignettei: Must have delrange[1] >= delrange[0]'
        return
      endif
    endif else begin
      print,'ERROR in vignettei: delrange must be a 2-element vector'
      return
    endelse
  endif

endelse

; Now process the image(s)

change_delay = mode eq 0 or n_elements(delrange) gt 0
change_doppler = mode eq 0 or n_elements(doprange) gt 0

print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this image

  if keyword_set(st) then begin
    image = (*stacki[n]).image
    dopbin1 = getextrai('dopbin1',stack=(n+1))
    dopbin2 = getextrai('dopbin2',stack=(n+1))
    delbin1 = getextrai('delbin1',stack=(n+1))
    delbin2 = getextrai('delbin2',stack=(n+1))
    crpix1 = getextrai('crpix1',stack=(n+1))
    crpix2 = getextrai('crpix2',stack=(n+1))
  endif else begin
    image = (*loadedi).image
    dopbin1 = getextrai('dopbin1')
    dopbin2 = getextrai('dopbin2')
    delbin1 = getextrai('delbin1')
    delbin2 = getextrai('delbin2')
    crpix1 = getextrai('crpix1')
    crpix2 = getextrai('crpix2')
  endelse
  f = dfreq[n]*(findgen(width[n]) - eph_col[n])       ; Doppler (Hz)
  d = delayunit[n]*(findgen(height[n]) - eph_row[n])  ; delay (usec)
  
  ; Assign the left and right limits according to the
  ; method specified via the mode parameter;
  ; be careful about roundoff error influencing the limits

  if mode le 1 then begin

    ; Convert delay-Doppler limits to image bins

    if change_delay then begin
      deledgebin1 = round(eph_row[n] + dellim1/delayunit[n] + 1e-7)
      deledgebin2 = round(eph_row[n] + dellim2/delayunit[n] - 1e-7)
    endif else begin
      deledgebin1 = 0L
      deledgebin2 = height[n] - 1
    endelse

    if change_doppler then begin
      dopedgebin1 = round(eph_col[n] + doplim1/dfreq[n] + 1e-7)
      dopedgebin2 = round(eph_col[n] + doplim2/dfreq[n] - 1e-7)
    endif else begin
      dopedgebin1 = 0L
      dopedgebin2 = width[n] - 1
    endelse

  endif else begin

    ; Directly specify image bins

    if change_delay then begin
      deledgebin1 = round(1.0*dellim1 + 1e-7)
      deledgebin2 = round(1.0*dellim2 - 1e-7)
    endif else begin
      deledgebin1 = 0L
      deledgebin2 = height[n] - 1
    endelse

    if change_doppler then begin
      dopedgebin1 = round(1.0*doplim1 + 1e-7)
      dopedgebin2 = round(1.0*doplim2 - 1e-7)
    endif else begin
      dopedgebin1 = 0L
      dopedgebin2 = width[n] - 1
    endelse

  endelse

  ; Make sure that the limits fall within the image and that
  ; more than two bins remain in each dimension after vignetting

  if dopedgebin1 lt 0 or dopedgebin2 ge width[n] then begin
    print,'ERROR on image #',n+1,': Requested Doppler bins [', $
          dopedgebin1,'-',dopedgebin2,'] extend beyond image [0-', $
          width[n]-1,']',format='(4(a,i0),a)'
    print,'       Image left unchanged'
  endif else if dopedgebin2 le dopedgebin1+1 then begin
    print,'ERROR on image #',n+1,': Requested Doppler bins [', $
          dopedgebin1,'-',dopedgebin2, $
          '] would leave too few points in vignetted image', $
          format='(3(a,i0),a)'
    print,'       Image left unchanged'
  endif else if deledgebin1 lt 0 or deledgebin2 ge height[n] then begin
    print,'ERROR on image #',n+1,': Requested delay bins [', $
          deledgebin1,'-',deledgebin2,'] extend beyond image [0-', $
          height[n]-1,']',format='(4(a,i0),a)'
    print,'       Image left unchanged'
  endif else if deledgebin2 le deledgebin1+1 then begin
    print,'ERROR on image #',n+1,': Requested delay bins [', $
          deledgebin1,'-',deledgebin2, $
          '] would leave too few points in vignetted image', $
          format='(3(a,i0),a)'
    print,'       Image left unchanged'
  endif else begin

    ; Everything's OK, so change the image and display the results
    ;
    ; Note that you have to replace the entire image structure: IDL won't
    ; let you change the dimensions of the existing image array.

    image = image[dopedgebin1:dopedgebin2,deledgebin1:deledgebin2]
    width[n] = dopedgebin2 - dopedgebin1 + 1
    height[n] = deledgebin2 - deledgebin1 + 1
    eph_col[n] = eph_col[n] - dopedgebin1
    eph_row[n] = eph_row[n] - deledgebin1
    if notnull(dc_col[n]) then dc_col[n] = dc_col[n] - dopedgebin1
    if notnull(dopbin1) and notnull(dopbin2) then begin
      newdopbin1 = dopbin1 - dopedgebin1
      newdopbin2 = dopbin2 - dopedgebin1
      if newdopbin1 lt 0 and newdopbin2 lt 0 then begin
        print,'WARNING on image #',n+1,': Signal range reset to leftmost column', $
              format='(a,i0,a)'
        dopbin1 = 0L
        dopbin2 = 0L
      endif else if newdopbin1 ge width[n] and newdopbin2 ge width[n] then begin
        print,'WARNING on spectrum #',n+1,': Signal range reset to rightmost column', $
              format='(a,i0,a)'
        dopbin1 = width[n] - 1
        dopbin2 = width[n] - 1
      endif else begin
        dopbin1 = newdopbin1 > 0L
        dopbin2 = newdopbin2 < (width[n] - 1)
      endelse
    endif else if notnull(dopbin1) then begin
      dopbin1 = (dopbin1 - dopedgebin1) > 0L
    endif else if notnull(dopbin2) then begin
      dopbin2 = (dopbin2 - dopedgebin1) < (width[n] - 1)
    endif
    if notnull(delbin1) and notnull(delbin2) then begin
      newdelbin1 = delbin1 - deledgebin1
      newdelbin2 = delbin2 - deledgebin1
      if newdelbin1 lt 0 and newdelbin2 lt 0 then begin
        print,'WARNING on image #',n+1,': Signal range reset to bottom row', $
              format='(a,i0,a)'
        delbin1 = 0L
        delbin2 = 0L
      endif else if newdelbin1 ge height[n] and newdelbin2 ge height[n] then begin
        print,'WARNING on spectrum #',n+1,': Signal range reset to top row', $
              format='(a,i0,a)'
        delbin1 = height[n] - 1
        delbin2 = height[n] - 1
      endif else begin
        delbin1 = newdelbin1 > 0L
        delbin2 = newdelbin2 < (height[n] - 1)
      endelse
    endif else if notnull(delbin1) then begin
      delbin1 = (delbin1 - deledgebin1) > 0L
    endif else if notnull(delbin2) then begin
      delbin2 = (delbin2 - deledgebin1) < (height[n] - 1)
    endif
    if notnull(crpix1) then crpix1 = crpix1 - dopedgebin1
    if notnull(crpix2) then crpix2 = crpix2 - deledgebin1

    if keyword_set(st) then begin
      setextrai,'f','eph_col',eph_col[n],stack=(n+1)
      setextrai,'f','eph_row',eph_row[n],stack=(n+1)
      if notnull(dc_col[n]) then setextrai,'f','dc_col',dc_col[n],stack=(n+1)
      if notnull(dopbin1) then setextrai,'i','dopbin1',dopbin1,stack=(n+1)
      if notnull(dopbin2) then setextrai,'i','dopbin2',dopbin2,stack=(n+1)
      if notnull(delbin1) then setextrai,'i','delbin1',delbin1,stack=(n+1)
      if notnull(delbin2) then setextrai,'i','delbin2',delbin2,stack=(n+1)
      if notnull(crpix1) then setextrai,'f','crpix1',crpix1,stack=(n+1)
      if notnull(crpix2) then setextrai,'f','crpix2',crpix2,stack=(n+1)
      oldStruc = *stacki[n]
      newStruc = {group:oldStruc.group,                       $
                  image:image, extratags:oldStruc.extratags,  $
                  width:width[n], height:height[n],           $
                  nextra:oldStruc.nextra, tname:oldStruc.tname, pol:oldStruc.pol}
      *stacki[n] = newStruc
    endif else begin
      setextrai,'f','eph_col',eph_col[n]
      setextrai,'f','eph_row',eph_row[n]
      if notnull(dc_col[n]) then setextrai,'f','dc_col',dc_col[n]
      if notnull(dopbin1) then setextrai,'i','dopbin1',dopbin1
      if notnull(dopbin2) then setextrai,'i','dopbin2',dopbin2
      if notnull(delbin1) then setextrai,'i','delbin1',delbin1
      if notnull(delbin2) then setextrai,'i','delbin2',delbin2
      if notnull(crpix1) then setextrai,'f','crpix1',crpix1
      if notnull(crpix2) then setextrai,'f','crpix2',crpix2
      oldStruc = *loadedi
      newStruc = {image:image, extratags:oldStruc.extratags,  $
                  width:width[n], height:height[n],           $
                  nextra:oldStruc.nextra, tname:oldStruc.tname, pol:oldStruc.pol}
      *loadedi = newStruc
    endelse

    printstring1 = 'Vignetting image #' + string(n+1,format='(i0)') + ' to '  $
                   + string(width[n],format='(i0)') + 'x'                     $
                   + string(height[n],format='(i0)') + ' pixels, '
    nskip1 = strlen(printstring1)
    printstring1 = printstring1 + 'Doppler bins '             $
                   + string(dopedgebin1,format='(i0)') + '-'  $
                   + string(dopedgebin2,format='(i0)')
    formatstring = '(' + string(nskip1,format='(i0)') + 'x,a)'
    printstring2 = string('delay   bins ',format=formatstring)  $
                   + string(deledgebin1,format='(i0)')  + '-'   $
                   + string(deledgebin2,format='(i0)')
    nskip2 = strlen(printstring1) - strlen(printstring2)
    if nskip2 gt 0 then begin
      for i=0L,nskip2-1 do printstring2 = printstring2 + ' '
    endif else if nskip2 lt 0 then begin
      for i=0L,-nskip2-1 do printstring1 = printstring1 + ' '
    endif
    printstring1 = printstring1 + '  ('                                               $
                   + string(f[dopedgebin1]-0.5*dfreq[n],format='(f8.2)') + ' Hz to '  $
                   + string(f[dopedgebin2]+0.5*dfreq[n],format='(f8.2)') + ' Hz)'
    printstring2 = printstring2 + '  ('                                       $
                   + string(d[deledgebin1]-0.5*delayunit[n],format='(f8.2)')  $
                   + ' us to '                                                $
                   + string(d[deledgebin2]+0.5*delayunit[n],format='(f8.2)')  $
                   + ' us)'
    print,printstring1
    print,printstring2

  endelse

endfor

print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fold,chan=chan,stack=st,help=help

; Fold the loaded pair, that is, make it symmetric and increase SNR by
; averaging it with a version of itself that's been flipped about zero Doppler
;
; If stack keyword is set, do this instead for all pairs in the stack
;
; Each flipped spectrum is shifted in order to get the 0-Hz bin to match onto
; itself, and this leads one edge to get wrapped around to the other edge.
; This should only be a problem if the signal takes up a large fraction of the
; spectrum or if xjcen is WAY off center.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

sqrt_2 = sqrt(2.0)

if n_params() ne 0 or keyword_set(help) then begin
  print,'fold[,chan=1 .. npol2][,/stack][,/help]'
  return
endif

if not keyword_set(st) then begin

  ; Fold the loaded pair

  if (*loaded).ndata le 2 then begin
    print,"ERROR: There's no loaded pair to fold"
    return
  endif

; Check which channels to fold

  npol = n_elements((*loaded).tags)
  foldpol = intarr(npol) 

  if n_elements(chan) eq 0 then begin
    foldpol = foldpol + 1
  endif else if chan ge 1 and chan le npol then begin
    foldpol[chan-1] = 1
  endif else begin
    print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
    return
  endelse

  ndata = (*loaded).ndata

  for ch=1,npol do begin
    if foldpol[ch-1] then begin
      xjcen = (*loaded).tags[ch-1].xjcen
      sdev = (*loaded).tags[ch-1].sdev
      orig = reform((*loaded).spec[ch-1,*])
      nshift = 2*xjcen - ndata + 1L
      flip = shift(reverse(orig), nshift)

      ; Increase SNR by root 2; decrease noise by root 2

      (*loaded).spec[ch-1,*] = (orig + flip)/sqrt_2  ;  that's sqrt_2 * (orig + flip)/2
      (*loaded).tags[ch-1].sdev = sdev/sqrt_2

      ; Get new signal limits by using the wider of the original and flipped limits

      orig_jsnr1 = (*loaded).tags[ch-1].jsnr1
      orig_jsnr2 = (*loaded).tags[ch-1].jsnr2
      orig_offset1 = xjcen - orig_jsnr1
      orig_offset2 = orig_jsnr2 - xjcen
      flip_jsnr1 = (xjcen - orig_offset2) > 0L
      flip_jsnr2 = (xjcen + orig_offset1) < (ndata - 1)
      (*loaded).tags[ch-1].jsnr1 = orig_jsnr1 < flip_jsnr1
      (*loaded).tags[ch-1].jsnr2 = orig_jsnr2 > flip_jsnr2
    endif
  endfor

  ; Add an extra tag to show that folding has been performed

  addextra,'s','folded','true',foldpol

endif else begin

  ; Fold all pairs in the stack

  if nstack eq 0 then begin
    print,'ERROR in fold: The pair stack is empty'
    return
  endif

  for n=0L,nstack-1 do begin
    ndata = (*stack[n]).ndata
; Check which channels to fold

    npol = n_elements((*stack[n]).tags)
    foldpol = intarr(npol) 
  
    if n_elements(chan) eq 0 then begin
      foldpol = foldpol + 1
    endif else if chan ge 1 and chan le npol then begin
      foldpol[chan-1] = 1
    endif else begin
      print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
      return
    endelse

    for ch=1,npol do begin
      if foldpol[ch-1] then begin
        xjcen = (*stack[n]).tags[ch-1].xjcen
        sdev = (*stack[n]).tags[ch-1].sdev
        orig = reform((*stack[n]).spec[ch-1,*])
        nshift = 2*xjcen - ndata + 1L
        flip = shift(reverse(orig), nshift)
        (*stack[n]).spec[ch-1,*] = (orig + flip)/sqrt_2
        (*stack[n]).tags[ch-1].sdev = sdev/sqrt_2
        orig_jsnr1 = (*stack[n]).tags[ch-1].jsnr1
        orig_jsnr2 = (*stack[n]).tags[ch-1].jsnr2
        orig_offset1 = xjcen - orig_jsnr1
        orig_offset2 = orig_jsnr2 - xjcen
        flip_jsnr1 = (xjcen - orig_offset2) > 0L
        flip_jsnr2 = (xjcen + orig_offset1) < (ndata - 1)
        (*stack[n]).tags[ch-1].jsnr1 = orig_jsnr1 < flip_jsnr1
        (*stack[n]).tags[ch-1].jsnr2 = orig_jsnr2 > flip_jsnr2
      endif
    endfor
    addextra,'s','folded','true',foldpol,stack=(n+1)
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fold1,stack=st,help=help

; Fold the loaded single-channel spectrum, that is, make it symmetric and 
; increase SNR by averaging it with a version of itself that's been flipped
; about zero Doppler
;
; If stack keyword is set, do this instead for all spectra in stack1
;
; The flipped spectrum is shifted in order to get the 0-Hz bin to match onto
; itself, and this leads one edge to get wrapped around to the other edge.
; This should only be a problem if the signal takes up a large fraction of the
; spectrum or if xjcen is WAY off center.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

sqrt_2 = sqrt(2.0)

if n_params() ne 0 or keyword_set(help) then begin
  print,'fold1[,/stack][,/help]'
  return
endif

if not keyword_set(st) then begin

  ; Fold the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,"ERROR: There's no loaded single-channel spectrum to fold"
    return
  endif

  ndata = (*loaded1).ndata
  xjcen = (*loaded1).tags.xjcen
  sdev = (*loaded1).tags.sdev
  orig = (*loaded1).spec
  nshift = 2*xjcen - ndata + 1L
  flip = shift(reverse(orig), nshift)

  ; Increase SNR by root 2; decrease noise by root 2

  (*loaded1).spec = (orig + flip)/sqrt_2  ;  that's sqrt_2 * (orig + flip)/2
  (*loaded1).tags.sdev = sdev/sqrt_2

  ; Get new signal limits by using the wider of the original and flipped limits

  orig_jsnr1 = (*loaded1).tags.jsnr1
  orig_jsnr2 = (*loaded1).tags.jsnr2
  orig_offset1 = xjcen - orig_jsnr1
  orig_offset2 = orig_jsnr2 - xjcen
  flip_jsnr1 = (xjcen - orig_offset2) > 0L
  flip_jsnr2 = (xjcen + orig_offset1) < (ndata - 1)
  (*loaded1).tags.jsnr1 = orig_jsnr1 < flip_jsnr1
  (*loaded1).tags.jsnr2 = orig_jsnr2 > flip_jsnr2

  ; Add an extra tag to show that folding has been performed

  addextra1,'s','folded','true'

endif else begin

  ; Fold all spectra in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in fold1: The single-channel stack is empty'
    return
  endif

  for n=0L,nstack1-1 do begin
    ndata = (*stack1[n]).ndata
    xjcen = (*stack1[n]).tags.xjcen
    sdev = (*stack1[n]).tags.sdev
    orig = (*stack1[n]).spec
    nshift = 2*xjcen - ndata + 1L
    flip = shift(reverse(orig), nshift)
    (*stack1[n]).spec = (orig + flip)/sqrt_2
    (*stack1[n]).tags.sdev = sdev/sqrt_2
    orig_jsnr1 = (*stack1[n]).tags.jsnr1
    orig_jsnr2 = (*stack1[n]).tags.jsnr2
    orig_offset1 = xjcen - orig_jsnr1
    orig_offset2 = orig_jsnr2 - xjcen
    flip_jsnr1 = (xjcen - orig_offset2) > 0L
    flip_jsnr2 = (xjcen + orig_offset1) < (ndata - 1)
    (*stack1[n]).tags.jsnr1 = orig_jsnr1 < flip_jsnr1
    (*stack1[n]).tags.jsnr2 = orig_jsnr2 > flip_jsnr2
    addextra1,'s','folded','true',stack=(n+1)
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro foldi,stack=st,help=help

; Fold the loaded image, that is, make it symmetric and increase SNR
; by averaging it with a version of itself that's been flipped about
; zero Doppler (actually about the center of the 0-Hz column, which
; might not be quite equal to the "eph_col" extra tag value).
;
; If stack keyword is set, do this instead for all images in stacki

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

sqrt_2 = sqrt(2.0)

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'foldi[,/stack][,/help]'
  print,' '
  print,'Fold the loaded image -- that is, average it with a version of'
  print,"     itself that's been flipped about zero Doppler (actually"
  print,'     about the center of the 0-Hz column, which might not'
  print,"     exactly coincide with the 'eph_col' extra tag value)"
  print,' '
  print,'/stack does this instead for all images in the image stack'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Fold the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,' '
    print,'ERROR in foldi: No image is loaded'
    print,' '
    return
  endif

  ; Get new 0-Hz column

  orig_eph_col = getextrai('eph_col')
  new_eph_col = 1.0*round(orig_eph_col)
  setextrai,'f','eph_col',new_eph_col

  ; Read the image, reverse it in Doppler, and shift it to keep
  ; the center of the ephemeris column in the same place

  width = (*loadedi).width
  nshift = 2*round(new_eph_col) - width + 1L
  origimage = (*loadedi).image
  flipimage = shift(reverse(origimage, 1), [nshift,0])

  ; Take the mean of the two, increase SNR by root 2, and decrease noise by root 2

  (*loadedi).image = (origimage + flipimage)/sqrt_2   ;  that's sqrt_2 * (orig + flip)/2
  sdev = getextrai('sdev')
  if notnull(sdev) then setextrai,'f','sdev',sdev/sqrt_2

  ; Get new signal limits by using the wider of the original and flipped limits

  orig_dopbin1 = getextrai('dopbin1')
  orig_dopbin2 = getextrai('dopbin2')
  if notnull(dopbin1) and notnull(dopbin2) then begin
    orig_offset1 = round(new_eph_col - orig_dopbin1)  ; offset from the column CENTER
    orig_offset2 = round(orig_dopbin2 - new_eph_col)
    flip_dopbin1 = round(new_eph_col - orig_offset2) > 0L
    flip_dopbin2 = round(new_eph_col + orig_offset1) < (width - 1)
    new_dopbin1 = orig_dopbin1 < flip_dopbin1
    new_dopbin2 = orig_dopbin2 > flip_dopbin2
    setextrai,'i','dopbin1',new_dopbin1
    setextrai,'i','dopbin2',new_dopbin2
  endif

  ; DC column is now meaningless (unless it coincided with eph_col)

  orig_dc_col = getextrai('dc_col')
  if notnull(orig_dc_col) then begin
    if orig_dc_col ne orig_eph_col then deleteextrai,'dc_col'
  endif

  ; Add an extra tag to show that folding has been performed

  addextrai,'s','folded','true'

endif else begin

  ; Fold all images in the stack

  if nstacki eq 0 then begin
    print,'ERROR in foldi: The image stack is empty'
    return
  endif

  for n=0L,nstacki-1 do begin
    orig_eph_col = getextrai('eph_col',stack=(n+1))
    new_eph_col = 1.0*round(orig_eph_col)
    setextrai,'f','eph_col',new_eph_col,stack=(n+1)
    width = (*stacki[n]).width
    nshift = 2*round(new_eph_col) - width + 1L
    origimage = (*stacki[n]).image
    flipimage = shift(reverse(origimage, 1), [nshift,0])
    (*stacki[n]).image = (origimage + flipimage)/sqrt_2
    sdev = getextrai('sdev',stack=(n+1))
    if notnull(sdev) then setextrai,'f','sdev',sdev/sqrt_2,stack=(n+1)
    orig_dopbin1 = getextrai('dopbin1',stack=(n+1))
    orig_dopbin2 = getextrai('dopbin2',stack=(n+1))
    if notnull(dopbin1) and notnull(dopbin2) then begin
      orig_offset1 = round(new_eph_col - orig_dopbin1)  ; offset from the column CENTER
      orig_offset2 = round(orig_dopbin2 - new_eph_col)
      flip_dopbin1 = round(new_eph_col - orig_offset2) > 0L
      flip_dopbin2 = round(new_eph_col + orig_offset1) < (width - 1)
      new_dopbin1 = orig_dopbin1 < flip_dopbin1
      new_dopbin2 = orig_dopbin2 > flip_dopbin2
      setextrai,'i','dopbin1',new_dopbin1,stack=(n+1)
      setextrai,'i','dopbin2',new_dopbin2,stack=(n+1)
    endif
    orig_dc_col = getextrai('dc_col',stack=(n+1))
    if notnull(orig_dc_col) then begin
      if orig_dc_col ne orig_eph_col then deleteextrai,'dc_col',stack=(n+1)
    endif
    addextrai,'s','folded','true',stack=(n+1)
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro bline,degree,omitf=farray,chan=chan,stack=st,nosubtract=nosub, $
    divide=div,plotbase=plotbase,help=help,_extra=_ext

; Subtract a polynomial baseline from the loaded pair
;
; If stack keyword is set, do this instead for all pairs in the stack
;
; If nosubtract is set, just compute and display the fit coefficients 
; but don't subtract anything from the spectrum.
;
; If divide is set, subtract the baseline and then normalize the result
; by dividing by the baseline
;
; If plotbase is set, overplot the fitted baseline on the current plot.
; The default is to use long dashes for OC and short for SC; this can
; be changed by specifying oplot's 'linestyle' keyword, color can
; be added by specifying oplot's 'color' keyword, and the line thickness
; can be changed by specifying oplot's 'thick' keyword.
;
; This routine uses polymaskfit, a modified version of IDL's poly_fit; this
; is a lot faster than svdfit, but allegedly is less stable.  Consider using svdfit
; (with the double keyword set, or else it's not very precise!) if this turns out
; to be a persistent problem.
;
; Modified 4/12/02 to replace polyfitw (obsolete starting with IDL 5.4)
; with poly_fit, and to shift/scale the x-axis to the range [-1,1] before
; fitting so as to reduce roundoff error.  Before this I was using
; frequency bin number [0,ndata-1] as the x-coordinate; the change results
; in changed coefficients being written as extra tags.
;
; Modified 8/9/02 to add the omitf=farray keyword option for omitting
; frequency ranges from the baseline fit (in addition to the signal
; region); this was spurred by the need to omit the imperfectly
; interpolated DC region when baselining some unhopped spectra.
;
; Modified 11/20/02 to switch to using polymaskfit rather than the
; built-in poly_fit, and to use the frequency vector corresponding to
; bin *centers* as the x-variable for fits.
;
; Modified 11/17/05 to add the plotbase keyword
;
; Modified 6/25/05 to recognize that frequency *already* refers to the center
; of each bin, not to the left edge

; Modified 01/11/22 For stokes: There is no offset to compute, we have to
; use the OC and SC ones. This is an unfortunate non-generalization that
; might be fixed with some more flags. However, this routine already knows
; that dividing is in the calibration step, so OK, I guess.
; So for stokes, it does the fit for the subtraction, but the division
; is sqrt(base_OC * base_SC)
;

 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that the sets of frequency ranges to be omitted from the baseline fit
; have been correctly input (if input at all).  Define the number of such
; regions to be the number input plus one: one extra for the signal region
; specified via tags jsnr1 - jsnr2, which is always omitted from the fit.

if n_elements(farray) eq 0 then begin
  n_omitregions = 1
  omitOK = 1
endif else begin
  omitsize = size(farray)
  omitdims = omitsize[0]
  dim1 = omitsize[1]
  if omitdims eq 0 or omitdims gt 2 or dim1 ne 2 then begin
    omitOK = 0
  endif else begin
    omitOK = 1
    n_omitregions = (omitdims eq 1) ? 2 : (1 + omitsize[2])
    for k=0L,n_omitregions-2 do begin
      if farray[0,k] gt farray[1,k] then omitOK = 0
    endfor
  endelse
endelse

if n_params() eq 0 then degree = 0  ;  just so it's defined

if n_params() ne 1 or degree lt 0 or keyword_set(help) $
                   or (keyword_set(nosub) and keyword_set(div)) $
                   or not omitOK then begin
  print,' '
  print,'bline,degree[,omitf=farray][,chan=1 .. npo][,/stack][,/nosubtract] $'
  print,'            [,/divide][,/plotbase[,oplot keywords]][,/help]'
  print,' '
  print,'        farray is an array of frequency pairs defining regions'
  print,'             (IN ADDITION to the signal region) which should be omitted'
  print,'             from the baseline fit.  Values are in Hz and must be'
  print,'             in increasing order even if the spectrum has frequency'
  print,'             increasing leftward.  Examples of farray are [-150,-120]'
  print,'             or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  print,'        /stack operates on all stack pairs instead of on the loaded pair'
  print,'        /nosubtract computes and displays the fit without subtracting'
  print,'        /divide subtracts the fit and then divides by the fit'
  print,' '
  print,'        /plotbase overplots the fitted baseline on the current plot.'
  print,'             The default is to use long dashes for OC and short for SC;'
  print,'             this can be changed, line thickness can be increased, and'
  print,"             color can be added by specifying oplot keywords 'linestyle'"
  print,"             'thick' and 'color'"
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Subtract a baseline from the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in bline: No pair is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Subtract a baseline from each pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in bline: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; Define the format code for displaying the fit coefficients

ndegree = round(degree)
coeff_format = strtrim(ndegree+1L, 2) + 'g'

; Process the loaded pair, or else loop through the pair stack

for n=0L,nuse-1 do begin

  ; Get a few elements of this pair, including the frequency vector
  ; corresponding to bin *centers*

  if keyword_set(st) then begin
    useStruc = *stack[n]
    startstring = 'Stack pair #' + string(n+1, format='(i0)') + ': '
  endif else begin
    useStruc = *loaded
    startstring = ''
  endelse
  ndata = useStruc.ndata
  tags = useStruc.tags
  npol = n_elements(useStruc.spec[*,0])
  df = tags[0].dfreq
  posfr = tags[0].posfr
  bin0 = tags[0].xjcen
  if keyword_set(st) then begin
    freq = posfr*df*(findgen(ndata) - bin0)
  endif else begin
    freq = useStruc.freq
  endelse
  f_leftmost = -posfr*df*(bin0 + 0.5)
  f_rightmost = posfr*df*(ndata - bin0 - 0.5)
  fmin = f_leftmost < f_rightmost
  fmax = f_leftmost > f_rightmost

; Check which channels to fit

  blinepol = intarr(npol)
  if (n_elements(chan) eq 0) then begin
    blinepol = blinepol + 1
  endif else if chan ge 1 and chan le npol then begin
    blinepol[chan-1] = 1
  endif else begin
    print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
    return
  endelse

  ; Mask out the frequency bin ranges to be omitted from the fit
  ; (other than the signal range, which could vary by polarization)

  mask_partial = intarr(ndata) + 1
  for k=0L,n_omitregions-2 do begin
    if posfr eq 1 then begin
      omit_left = bin0 + round(farray[0,k]/df)
      omit_right = bin0 + round(farray[1,k]/df)
    endif else begin
      omit_left = bin0 - round(farray[1,k]/df)
      omit_right = bin0 - round(farray[0,k]/df)
    endelse
    if omit_left lt 0 or omit_right ge ndata then begin 
      print,' '
      print,'WARNING for pair #',n+1,': Omitted frequency range (', $
            farray[0,k],' Hz to ',farray[1,k], $
            ' Hz) extends beyond the spectrum (',fmin,' Hz to ', $
            fmax,' Hz)',format='(a,i0,a,4(f9.2,a))'
    endif
    if omit_left lt ndata and omit_right ge 0 then begin
      omit_left = omit_left > 0
      omit_right = omit_right < (ndata - 1)
      mask_partial[omit_left:omit_right] = 0
    endif
  endfor

  ; Loop through the polarization channels and do the baseline fit(s)
  ; and (if specified) subtraction and division

  for ch=1,npol do begin
    if blinepol[ch-1] then begin

      ; Get the spectrum and the signal range for this channel,
      ; then mask out the signal range

      spec = reform(useStruc.spec[ch-1,*])
      in_mask = mask_partial
      omit_left = tags[ch-1].jsnr1
      omit_right = tags[ch-1].jsnr2
      in_mask[omit_left:omit_right] = 0
      n_noise = total(in_mask)

      ; If there's any baseline to fit, go ahead and fit it, then
      ; (if specified) subtract it and (if also specified) divide by it

      if n_noise eq 0 then begin
        print,startstring,'Channel ',ch, $
              ' has no baseline points, so no baseline was computed', $
              format='(2a,i0,a)'
      endif else begin
        fit_coeffs = polymaskfit(freq,spec,ndegree,in_mask=in_mask,yfit=baseline)
        if keyword_set(plotbase) then oplot,freq,baseline,linestyle=(8-3*ch),_extra=_ext
        print,startstring,'Channel ',ch,' fit coefficients (',n_noise,' points) =', $
              float(fit_coeffs),format='(2a,i0,a,i0,a,'+coeff_format+')'
        if not keyword_set(nosub) then begin
          newspec = spec - baseline
          if not keyword_set(div) then begin
            print,startstring,'Channel ',ch,' baseline subtracted',format='(2a,i0,a)'

            ; If a baseline was subtracted, add an extra tag for each fit coefficient
            ; for each polarization channel processed.
            ;
            ; (Don't do this if normalization was performed, since this is a part of
            ; initial data reduction and no one wants a record of it.)

            fit_coeffs_str = string(float(fit_coeffs), format='(g)')
            chanflags = intarr[npol]
            chanflags[ch-1] = 1
            for i=0L,ndegree do begin
              etagstem = 'bline_c' + strtrim(i, 2)
              if keyword_set(st) then begin
                addextra,'f',etagstem,fit_coeffs_str[i],chanflags,stack=(n+1)
              endif else begin
                addextra,'f',etagstem,fit_coeffs_str[i],chanflags
              endelse
            endfor
          endif else begin
; Stokes stuff
            if ch eq 1 then xxbase = double(baseline)
            if ch eq 2 then begin
              xxbase = float(sqrt(xxbase * double(baseline)))
              xxbad = where(~finite(xxbase),count)
              if count gt 0 then xxbase[xxbad] = -1
            endif
            if ch gt 2 and total(blinepol) gt 1 then baseline = xxbase
; end of Stokes stuff
            not_pos = where(baseline le 0.0, count)
            if count gt 0 then begin
              print,startstring,"Channel ",ch, $
                    " baseline has non-positive values, so can't divide", $
                    format='(2a,i0,a)'
            endif else begin
              newspec = newspec/baseline
              print,startstring,'Channel ',ch, $
                    ' -- subtracted baseline, then divided by baseline', $
                    format='(2a,i0,a)'
            endelse
          endelse
          if keyword_set(st) then begin
            (*stack[n]).spec[ch-1,*] = float(newspec)
          endif else begin
            (*loaded).spec[ch-1,*] = float(newspec)
          endelse
        endif
      endelse
    endif
  endfor
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro bline1,degree,omitf=farray,stack=st,nosubtract=nosub,divide=div, $
           plotbase=plotbase,help=help,_extra=_ext

; Subtract a polynomial baseline from the loaded single-channel spectrum 
;
; If stack keyword is set, do this instead for all spectra in the
; single-channel stack
;
; If nosubtract is set, just compute and display the fit coefficients 
; but don't subtract anything from the spectrum.
;
; If divide is set, subtract the baseline and then normalize the result
; by dividing by the baseline
;
; If plotbase is set, overplot the fitted baseline on the current plot.
; The default is to use long dashes; this can be changed by specifying
; oplot's 'linestyle' keyword, color can be added by specifying oplot's
; 'color' keyword, and the line thickness can be changed by specifying
; oplot's 'thick' keyword.
;
; This routine uses polymaskfit, a modified version of IDL's poly_fit; this
; is a lot faster than svdfit, but allegedly is less stable.  Consider using svdfit
; (with the double keyword set, or else it's not very precise!) if this turns out
; to be a persistent problem.
;
; Modified 4/12/02 to replace polyfitw (obsolete starting with IDL 5.4)
; with poly_fit, and to shift/scale the x-axis to the range [-1,1] before
; fitting so as to reduce roundoff error.  Before this I was using
; frequency bin number [0,ndata-1] as the x-coordinate; the change results
; in changed coefficients being written as extra tags.
;
; Modified 8/9/02 to add the omitf=farray keyword option for omitting
; frequency ranges from the baseline fit (in addition to the signal
; region); this was spurred by the need to omit the imperfectly
; interpolated DC region when baselining some unhopped spectra.
;
; Modified 11/20/02 to switch to using polymaskfit rather than the
; built-in poly_fit, and to use the frequency vector corresponding to
; bin *centers* as the x-variable for fits.
;
; Modified 11/17/05 to add the plotbase keyword
;
; Modified 6/25/05 to recognize that frequency *already* refers to the center
; of each bin, not to the left edge
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that the sets of frequency ranges to be omitted from the baseline fit
; have been correctly input (if input at all).  Define the number of such
; regions to be the number input plus one: one extra for the signal region
; specified via tags jsnr1 - jsnr2, which is always omitted from the fit.

if n_elements(farray) eq 0 then begin
  n_omitregions = 1
  omitOK = 1
endif else begin
  omitsize = size(farray)
  omitdims = omitsize[0]
  dim1 = omitsize[1]
  if omitdims eq 0 or omitdims gt 2 or dim1 ne 2 then begin
    omitOK = 0
  endif else begin
    omitOK = 1
    n_omitregions = (omitdims eq 1) ? 2 : (1 + omitsize[2])
    for k=0L,n_omitregions-2 do begin
      if farray[0,k] gt farray[1,k] then omitOK = 0
    endfor
  endelse
endelse

if n_params() eq 0 then degree = 0  ;  just so it's defined

if n_params() ne 1 or degree lt 0 or keyword_set(help) $
                   or (keyword_set(nosub) and keyword_set(div)) then begin
  print,' '
  print,'bline1,degree[,omitf=farray][,/stack][,/nosubtract][,/divide] $'
  print,'             [,/plotbase[,oplot keywords]][,/help]'
  print,' '
  print,'        farray is an array of frequency pairs defining regions'
  print,'             (IN ADDITION to the signal region) which should be omitted'
  print,'             from the baseline fit.  Values are in Hz and must be'
  print,'             in increasing order even if the spectrum has frequency'
  print,'             increasing leftward.  Examples of farray are [-150,-120]'
  print,'             or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  print,'        /stack operates on all stack1 spectra instead of on the loaded'
  print,'             single-channel spectrum'
  print,'        /nosubtract computes and displays the fit without subtracting'
  print,'        /divide subtracts the fit and then divides by the fit'
  print,' '
  print,'        /plotbase overplots the fitted baseline on the current plot.'
  print,'             The default is to use long dashes; this can be changed,'
  print,'             line thickness can be increased, and color can be added'
  print,"             by specifying oplot keywords 'linestyle' 'thick' and 'color'"
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Subtract a baseline from the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,"ERROR in bline1: There's no loaded single-channel spectrum"
    return
  endif
  nuse = 1L
endif else begin

  ; Subtract a baseline from each spectrum in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in bline1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
endelse

; Define the format code for displaying the fit coefficients

ndegree = round(degree)
coeff_format = strtrim(ndegree+1L, 2) + 'g'

; Process the loaded single-channel spectrum, or else loop through
; the single-channel stack

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum, including the frequency
  ; vector corresponding to bin *centers*

  if keyword_set(st) then begin
    useStruc = *stack1[n]
    startstring = 'Stack1 spectrum #' + string(n+1, format='(i0)') + ': '
  endif else begin
    useStruc = *loaded1
    startstring = ''
  endelse
  ndata = useStruc.ndata
  tags1 = useStruc.tags
  df = tags1.dfreq
  posfr = tags1.posfr
  bin0 = tags1.xjcen
  if keyword_set(st) then begin
    freq = posfr*df*(findgen(ndata) - bin0)
  endif else begin
    freq = useStruc.freq
  endelse
  f_leftmost = -posfr*df*(bin0 + 0.5)
  f_rightmost = posfr*df*(ndata - bin0 - 0.5)
  fmin = f_leftmost < f_rightmost
  fmax = f_leftmost > f_rightmost
  spec = useStruc.spec

  ; Mask out the frequency bin ranges to be omitted from the fit

  in_mask = intarr(ndata) + 1
  for k=0L,n_omitregions-2 do begin
    if posfr eq 1 then begin
      omit_left = bin0 + round(farray[0,k]/df)
      omit_right = bin0 + round(farray[1,k]/df)
    endif else begin
      omit_left = bin0 - round(farray[1,k]/df)
      omit_right = bin0 - round(farray[0,k]/df)
    endelse
    if omit_left lt 0 or omit_right ge ndata then begin 
      print,' '
      print,'WARNING for spectrum #',n+1,': Omitted frequency range (', $
            farray[0,k],' Hz to ',farray[1,k], $
            ' Hz) extends beyond the spectrum (',fmin,' Hz to ', $
            fmax,' Hz)',format='(a,i0,a,4(f9.2,a))'
    endif
    if omit_left lt ndata and omit_right ge 0 then begin
      omit_left = omit_left > 0
      omit_right = omit_right < (ndata - 1)
      in_mask[omit_left:omit_right] = 0
    endif
  endfor
  omit_left = tags1.jsnr1
  omit_right = tags1.jsnr2
  in_mask[omit_left:omit_right] = 0
  n_noise = total(in_mask)

  ; If there's any baseline to fit, go ahead and fit it, then
  ; (if specified) subtract it and (if also specified) divide by it

  if n_noise eq 0 then begin
    print,startstring,'Spectrum has no baseline points, so no baseline was computed'
  endif else begin
    fit_coeffs = polymaskfit(freq,spec,ndegree,in_mask=in_mask,yfit=baseline)
    if keyword_set(plotbase) then oplot,freq,baseline,linestyle=5,_extra=_ext
    print,startstring,'Fit coefficients (',n_noise,' points) =', $
          float(fit_coeffs),format='(2a,i0,a,'+coeff_format+')'
    if not keyword_set(nosub) then begin
      newspec = spec - baseline
      if not keyword_set(div) then begin
        print,startstring,'Subtracted baseline'

        ; If a baseline was subtracted, add an extra tag for each fit coefficient
        ;
        ; (Don't do this if normalization was performed, since this is a part of
        ; initial data reduction and no one wants a record of it.)

        fit_coeffs_str = string(float(fit_coeffs), format='(g)')
        for i=0L,ndegree do begin
          etagstem = 'bline1_c' + strtrim(i, 2)
          if keyword_set(st) then begin
            addextra1,'f',etagstem,fit_coeffs_str[i],stack=(n+1)
          endif else begin
            addextra1,'f',etagstem,fit_coeffs_str[i]
          endelse
        endfor
      endif else begin
        not_pos = where(baseline le 0.0, count)
        if count gt 0 then begin
          print,startstring,"Baseline has non-positive values, so can't divide"
        endif else begin
          newspec = newspec/baseline
          print,startstring,'Subtracted baseline, then divided by baseline'
        endelse
      endelse
      if keyword_set(st) then begin
        (*stack1[n]).spec = float(newspec)
      endif else begin
        (*loaded1).spec = float(newspec)
      endelse
    endif
  endelse
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro smoothf,efb,n=n,gauss=gauss,chan=chan,stack=st,help=help

; Smooth the loaded pair or the stack pairs to a specified
; effective frequency resolution efb.  This is done by convolving each
; spectrum with a filter whose shape is the same as the shape of a
; radar echo from a spherical target with a (cos theta)^n scattering law:
;
;    S(f) = amplitude * [ 1 - (2f/B)^2 ]^(n/2)  ;  |f| < B/2
;
; Here B is the full width of the filter; the amplitude is chosen (by function
; makefilter) so as to normalize the filter.  B is chosen such that the filter's
; equivalent bandwidth is equal to efb.
;
; Note that power-law exponent n needn't be an integer.  If not specified,
; n is set to 2 (Lambertian scattering).
;
; 2013 Jun 12: add /gauss keyword

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'smoothf,efb[,n=n OR /gauss][,chan=1 or 2][,/stack][,/help]'
  print,' '
  print,'Smooth the loaded pair to specified effective frequency resolution efb'
  print,'     by convolving with a filter whose shape is that of echoes from'
  print,'     a spherical target with a (cos theta)^n scattering law:'
  print,' '
  print,'          S(f) = amplitude * [ 1 - (2f/B)^2 ]^(n/2)  ;  |f| < B/2'
  print,' '
  print,"     Full width B is chosen such that the filter's equivalent bandwidth"
  print,'     is equal to efb.'
  print,' '
  print,"     n is the power-law exponent, which needn't be an integer:"
  print,'          n = 0.0 yields a rectangular filter;'
  print,'          large n yields a narrow central "spike" with weak "tails"'
  print,'          (default = 2.0)'
  print,' '
  print,'     /gauss uses a Gaussian filter rather than the one described above'
  print,' '
  print,'     /stack smooths each stack pair instead of the loaded pair'
  print,' '
  return
endif else if efb le 0 then begin
  print,'ERROR in smoothf: Effective resolution must be positive'
  return
endif else if n_elements(n) gt 0 and keyword_set(gauss) then begin
  print,'ERROR in smoothf: /gauss cannot be set if n is specified'
  return
endif

if n_elements(n) eq 0 then n = 2.0
if n lt 0 then begin
  print,"ERROR in smoothf: Power-law exponent n can't be negative"
  return
endif

; Do the smoothing

efb = 1.0*efb
n = 1.0*n

if not keyword_set(st) then begin

  ; Smooth the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in smoothf: No pair is loaded'
    return
  endif
  npol = n_elements((*loaded).spec[*,0])
    
  ndata = (*loaded).ndata
  dfreq = (*loaded).tags[0].dfreq
  freq = (*loaded).freq
  chanWasSmoothed = intarr(npol)

; Check which channels to smooth

  smoothpol = intarr(npol)
  if (n_elements(chan) eq 0) then begin
    smoothpol = smoothpol + 1
  endif else begin
    if chan lt 1 || chan gt npol then begin
      print, "ERROR in smoothf: request chan larger than number of channels in pair"
      return
    endif
    smoothpol[chan-1] = 1
  endelse
  
  ; Proceed according to the type of smoothing requested

  if keyword_set(gauss) then begin

    ; Gaussian smoothing

    fwhm = sqrt(2*alog(2)/!DPI)*efb
    fracNoiseReduc = sqrt(efb/dfreq)  ;  fractional noise reduction
    for ch=1,npol do begin
      if smoothpol[ch-1] then begin
        spec = reform((*loaded).spec[ch-1,*])
        sdev = (*loaded).tags[ch-1].sdev
        (*loaded).spec[ch-1,*] = fracNoiseReduc*gaussfold_cm(freq,spec,fwhm)
        (*loaded).tags[ch-1].sdev = sdev/fracNoiseReduc
        chanWasSmoothed[ch-1] = 1
      endif
    endfor

    ; Add an extra tag listing the new effective resolution

    addextra,'f','efb',efb,chanWasSmoothed

  endif else begin

    ; Matched-filter smoothing (for a spherical asteroid)

    ; Construct the smoothing filter
    ; (If there's a problem, makefilter prints a specific error message and
    ;  returns [0])

    filter = makefilter(efb,n,dfreq,ndata,efb_use)
    fracNoiseReduc = sqrt(efb_use/dfreq)  ;  fractional noise reduction

    ; If no problem, convolve the filter with the spectra and reduce sdev

    if total(filter) gt 0 then begin
      for ch=1,npol do begin
        if smoothpol[ch-1] then begin
          spec = reform((*loaded).spec[ch-1,*])
          sdev = (*loaded).tags[ch-1].sdev
          (*loaded).spec[ch-1,*] = fracNoiseReduc*convol(spec,filter,/edge_wrap)
          (*loaded).tags[ch-1].sdev = sdev/fracNoiseReduc
          chanWasSmoothed[ch-1] = 1
        endif
      endfor

      ; Add an extra tag listing the new effective resolution

      addextra,'f','efb',efb_use,chanWasSmoothed

    endif else begin
      print,'No channels of the loaded pair were smoothed'
    endelse
  endelse

endif else begin

  ; Smooth the pairs in the stack
  ;
  ; If some pairs in the stack have raw resolution coarser than the requested
  ; effective resolution, they'll be left alone and error messages will be
  ; printed, while the other pairs will be smoothed.

  if nstack eq 0 then begin
    print,'ERROR in smoothf: The pair stack is empty'
    return
  endif

  for k=0L,nstack-1 do begin
    ndata = (*stack[k]).ndata
    dfreq = (*stack[k]).tags[0].dfreq
    posfr = (*stack[k]).tags[0].posfr
    xjcen = (*stack[k]).tags[0].xjcen
    npol = n_elements((*stack[k]).spec[*,0])
; Check which channels to smooth

    smoothpol = intarr(npol)
    if (n_elements(chan) eq 0) then begin
      smoothpol = smoothpol + 1
    endif else begin
      if chan lt 1 || chan gt npol then begin
        print, "ERROR in smoothf: request chan larger than number of channels in pair"
        return
      endif
      smoothpol[chan-1] = 1
    endelse
    freq = posfr*dfreq*(findgen(ndata) - xjcen)
    chanWasSmoothed = intarr(npol)
    if keyword_set(gauss) then begin
      fwhm = sqrt(2*alog(2)/!DPI)*efb
      fracNoiseReduc = sqrt(efb/dfreq)
      for ch=1,npol do begin
        if smoothpol[ch-1] then begin
          spec = reform((*stack[k]).spec[ch-1,*])
          sdev = (*stack[k]).tags[ch-1].sdev
          (*stack[k]).spec[ch-1,*] = fracNoiseReduc*gaussfold_cm(freq,spec,fwhm)
          (*stack[k]).tags[ch-1].sdev = sdev/fracNoiseReduc
          chanWasSmoothed[ch-1] = 1
        endif
      endfor
      addextra,'f','efb',efb,chanWasSmoothed,stack=(k+1)
    endif else begin
      filter = makefilter(efb,n,dfreq,ndata,efb_use)
      fracNoiseReduc = sqrt(efb_use/dfreq)
      if total(filter) gt 0 then begin
        for ch=1,npol do begin
          if smoothpol[ch-1] then begin
            spec = reform((*stack[k]).spec[ch-1,*])
            sdev = (*stack[k]).tags[ch-1].sdev
            (*stack[k]).spec[ch-1,*] = fracNoiseReduc*convol(spec,filter,/edge_wrap)
            (*stack[k]).tags[ch-1].sdev = sdev/fracNoiseReduc
            chanWasSmoothed[ch-1] = 1
          endif
        endfor
        addextra,'f','efb',efb_use,chanWasSmoothed,stack=(k+1)
      endif else begin
        print,'No channels of stack pair #',k+1,' were smoothed',format='(a,i0,a)'
      endelse
    endelse
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro smoothf1,efb,n=n,gauss=gauss,chan=chan,stack=st,help=help

; Smooth the loaded single-channel spectrum or the single-channel stack spectra
; to a specified effective frequency resolution efb.  This is done by convolving
; each spectrum with a filter whose shape is the same as the shape of a
; radar echo from a spherical target with a (cos theta)^n scattering law:
;
;    S(f) = amplitude * [ 1 - (2f/B)^2 ]^(n/2)  ;  |f| < B/2
;
; Here B is the full width of the filter; the amplitude is chosen (by function
; makefilter) so as to normalize the filter.  B is chosen such that the filter's
; equivalent bandwidth is equal to efb.
;
; Note that power-law exponent n needn't be an integer.  If not specified,
; n is set to 2 (Lambertian scattering).
;
; 2013 Jun 12: add /gauss keyword

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'smoothf1,efb[,n=n OR /gauss][,chan=1 or 2][,/stack][,/help]'
  print,' '
  print,'Smooth the loaded single-channel spectrum to specified effective frequency'
  print,'     resolution efb by convolving with a filter whose shape is that of'
  print,'     echoes from a spherical target with a (cos theta)^n scattering law:'
  print,' '
  print,'          S(f) = amplitude * [ 1 - (2f/B)^2 ]^(n/2)  ;  |f| < B/2'
  print,' '
  print,"     Full width B is chosen such that the filter's equivalent bandwidth"
  print,'     is equal to efb.'
  print,' '
  print,"     n is the power-law exponent, which needn't be an integer:"
  print,'          n = 0.0 yields a rectangular filter;'
  print,'          large n yields a narrow central "spike" with weak "tails"'
  print,'          (default = 2.0)'
  print,' '
  print,'     /gauss uses a Gaussian filter rather than the one described above'
  print,' '
  print,'     /stack smooths each spectrum in the single-channel stack instead of'
  print,'          the loaded single-channel spectrum'
  print,' '
  return
endif else if efb le 0 then begin
  print,'ERROR in smoothf1: Effective resolution must be positive'
  return
endif else if n_elements(n) gt 0 and keyword_set(gauss) then begin
  print,'ERROR in smoothf1: /gauss cannot be set if n is specified'
  return
endif

if n_elements(n) eq 0 then n = 2.0
if n lt 0 then begin
  print,"ERROR in smoothf1: Power-law exponent n can't be negative"
  return
endif

; Do the smoothing

efb = 1.0*efb
n = 1.0*n

if not keyword_set(st) then begin

  ; Smooth the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in smoothf1: No single-channel spectrum is loaded'
    return
  endif

  ndata = (*loaded1).ndata
  dfreq = (*loaded1).tags.dfreq
  freq = (*loaded1).freq

  spec = (*loaded1).spec
  sdev = (*loaded1).tags.sdev

  ; Proceed according to the type of smoothing requested

  if keyword_set(gauss) then begin

    ; Gaussian smoothing

    fwhm = sqrt(2*alog(2)/!DPI)*efb
    fracNoiseReduc = sqrt(efb/dfreq)  ;  fractional noise reduction
    (*loaded1).spec = fracNoiseReduc*gaussfold_cm(freq,spec,fwhm)
    (*loaded1).tags.sdev = sdev/fracNoiseReduc

    ; Add an extra tag listing the new effective resolution

    addextra1,'f','efb',efb

  endif else begin

    ; Matched-filter smoothing (for a spherical asteroid)

    ; Construct the smoothing filter
    ; (If there's a problem, makefilter prints a specific error message and
    ;  returns [0])

    filter = makefilter(efb,n,dfreq,ndata,efb_use)
    fracNoiseReduc = sqrt(efb_use/dfreq)  ;  fractional noise reduction

    ; If no problem, convolve the filter with the spectrum and reduce sdev

    if total(filter) gt 0 then begin
      (*loaded1).spec = fracNoiseReduc*convol(spec,filter,/edge_wrap)
      (*loaded1).tags.sdev = sdev/fracNoiseReduc

      ; Add an extra tag listing the new effective resolution

      addextra1,'f','efb',efb_use

    endif else begin
      print,'The loaded single-channel spectrum was not smoothed'
    endelse
  endelse

endif else begin

  ; Smooth the spectra in the single-channel stack
  ;
  ; If some spectra in stack1 have raw resolution coarser than the requested
  ; effective resolution, they'll be left alone and error messages will be
  ; printed, while the other spectra will be smoothed.

  if nstack1 eq 0 then begin
    print,'ERROR in smoothf1: The single-channel stack is empty'
    return
  endif

  for k=0L,nstack1-1 do begin
    ndata = (*stack1[k]).ndata
    dfreq = (*stack1[k]).tags.dfreq
    posfr = (*stack1[k]).tags.posfr
    xjcen = (*stack1[k]).tags.xjcen
    freq = posfr*dfreq*(findgen(ndata) - xjcen)
    spec = (*stack1[k]).spec
    sdev = (*stack1[k]).tags.sdev
    if keyword_set(gauss) then begin
      fwhm = sqrt(2*alog(2)/!DPI)*efb
      fracNoiseReduc = sqrt(efb/dfreq)
      (*stack1[k]).spec = fracNoiseReduc*gaussfold_cm(freq,spec,fwhm)
      (*stack1[k]).tags.sdev = sdev/fracNoiseReduc
      addextra1,'f','efb',efb,stack=(k+1)
    endif else begin
      filter = makefilter(efb,n,dfreq,ndata,efb_use)
      fracNoiseReduc = sqrt(efb_use/dfreq)
      if total(filter) gt 0 then begin
        (*stack1[k]).spec = fracNoiseReduc*convol(spec,filter,/edge_wrap)
        (*stack1[k]).tags.sdev = sdev/fracNoiseReduc
        addextra1,'f','efb',efb_use,stack=(k+1)
      endif else begin
        print,'Single-channel stack spectrum #',k+1,' was not smoothed', $
              format='(a,i0,a)'
      endelse
    endelse
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro addextra,extraformat,extraname,extravalue,chanflags, $
             comment=extracomment,stack=n,help=help

; Add an extra tag to the loaded pair, or else to stack pair n if the 
; stack keyword is set
;
; chanflags is an optional two-element vector containing 0 and/or 1.  If only one
; channel is involved (chanflags ne [1,1]) then an appropriate suffix is added
; to the tag name.
;
; Floating-point and double-precision values are formatted with the highest
; required precision (format='(f0)').  If some other precision is desired,
; use the "string" function to create an appropriately formatted string,
; then pass that string as the extravalue argument.
;
; May 2002: Added the comment keyword for specifying the comment field

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

if n_elements(extraformat) eq 0 then extraformat = ''  ; just so it's defined
if n_elements(extraname) eq 0 then extraname = ''
if n_elements(extracomment) eq 0 then extracomment = ''
if keyword_set(help) or n_params() lt 3 or n_params() gt 4 $ 
             or size(extraformat, /type) ne 7 or size(extraname, /type) ne 7 $
             or size(extracomment, /type) ne 7 then begin
  print,' '
  print,'addextra,format,name,value[,chanflags][,comment=comment][,stack=n][,/help]'
  print,'         format, name, and comment must be quoted strings'
  print,'         chanflags = [1,0] if tag is relevant to OC, [0,1] if SC, ', $
        ' [1,1] if both (default), [0,0] if neither (no tag added)',format='(2a)'
  print,'         stack = n adds an extra tag to stack pair #n rather than ', $
        'to the loaded pair',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in addextra: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in addextra: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in addextra: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif else npol = n_elements((*stack[n-1]).tags)
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in addextra: No pair is loaded'
  return
endif else npol = n_elements((*loaded).tags)

if isnull(extraname) and notnull(extraformat) then begin
  print,'ERROR in addextra: Must specify a tag name'
  return
endif else if isnull(extraname) and isnull(extraformat) and $
              isnull(extravalue) and isnull(extracomment) then begin
  print,"ERROR in addextra: Can't specify a blank comment line"
  return
endif

if total(chanflags) eq 0 then begin
  return  ; no extra tags to set!
endif

if n_params() lt 4 then chanflags = intarr(npol) + 1

if n_elements(chanflags) lt 2 then begin
  c = chanflags
  chanflags = intarr(npol)
  chanflags[c-1] = 1
endif

if total(chanflags) eq 0 then begin
  print, 'WARNING in addextra: No channels to change'
  return
endif


; Construct a new structure containing the extra tag

chanst1 = chanstrings
chanst2 = chanstrings
for i=0, maxchan-1 do begin
  chanst1[i] = '_' + chanstrings[i]
  chanst2[i] = '# ' + chanstrings[i] + ': '
endfor
  
newcomment = strtrim(extracomment,2)
if strmid(newcomment,0,1) eq '#' then newcomment = strtrim(strmid(newcomment,1), 2)
newcomment = (isnull(newcomment)) ? '' : '# ' + newcomment
newformat = strtrim(extraformat,2)
newname = strtrim(extraname,2)

valuetype = size(extravalue, /type)
if valuetype eq 4 or valuetype eq 5 then begin
  newvalue = string(extravalue, format='(f0)')  ;  float or double
endif else begin
  newvalue = strtrim(extravalue,2)
endelse

; Add the new extra tag to the loaded pair or to the specified stack pair
;
; Have to replace the entire pair structure: IDL won't let you change the
; dimensions of the existing extratags array.
;
; IDL 5.3 is a pain about concatenating anonymous structures, so it has to be
; done in a tedious fashion.

if n_elements(n) eq 0 then oldStruc = *loaded else oldStruc = *stack[n-1]
extratags = reform(oldStruc.extratags)
npol = n_elements(oldStruc.spec[*,0])
nextra = oldStruc.nextra
;create the tags
allchan = (total(chanflags) eq npol)
whichchan = where(chanflags)
for i = 0, (allchan) ? 0 : total(chanflags)-1 do begin
  cnewname = newname
  cnewcomment = newcomment
  if ~allchan then begin
    if notnull(newname) then begin
      cnewname = newname + chanst1[whichchan[i]]
    endif else begin
      if strpos(newcomment,chanst2[whichchan[i]]) ne 0 then $
          cnewcomment = chanst2[whichchan[i]] + strmid(newcomment,2)
    endelse
  endif ; if allchan
  newextra = {format:newformat, name:cnewname, value:newvalue, comment:cnewcomment}
  extratags = [extratags, extratags[0]]
  for k=0L,3 do extratags[nextra].(k) = newextra.(k)
  nextra = nextra + 1L
  newStruc = {freq:oldStruc.freq, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname}
  if n_elements(n) eq 0 then *loaded = newStruc else *stack[n-1] = newStruc
  oldStruc = newStruc
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro addextra1,extraformat,extraname,extravalue, $
              comment=extracomment,stack=n,help=help

; Add an extra tag to the loaded single-channel spectrum, or else to
; single-channel stack spectrum n if the stack keyword is set
;
; Floating-point and double-precision values are formatted with the highest
; required precision (format='(f0)').  If some other precision is desired,
; use the "string" function to create an appropriately formatted string,
; then pass that string as the extravalue argument.
;
; May 2002: Added the comment keyword for specifying the comment field

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(extraformat) eq 0 then extraformat = ''  ; just so it's defined
if n_elements(extraname) eq 0 then extraname = ''
if n_elements(extracomment) eq 0 then extracomment = ''

if keyword_set(help) or n_params() ne 3 or size(extraformat, /type) ne 7 $
                     or size(extraname, /type) ne 7 $
                     or size(extracomment, /type) ne 7 then begin
  print,' '
  print,'addextra1,format,name,value[,comment=comment][,stack=n][,/help]'
  print,'          format, name, and comment must be quoted strings'
  print,'          stack = n adds an extra tag to stack1 spectrum #n rather ', $
        'than to the loaded single-channel spectrum',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in addextra1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in addextra1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in addextra1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in addextra1: No single-channel spectrum is loaded'
  return
endif else if isnull(extraname) and notnull(extraformat) then begin
  print,'ERROR in addextra1: Must specify a tag name'
  return
endif else if isnull(extraname) and isnull(extraformat) and $
              isnull(extravalue) and isnull(extracomment) then begin
  print,"ERROR in addextra1: Can't specify a blank comment line"
  return
endif

; Construct a new structure containing the extra tag

newcomment = strtrim(extracomment,2)
if strmid(newcomment,0,1) eq '#' then newcomment = strtrim(strmid(newcomment,1), 2)
newcomment = (isnull(newcomment)) ? '' : '# ' + newcomment
newformat = strtrim(extraformat,2)
newname = strtrim(extraname,2)
valuetype = size(extravalue, /type)
if valuetype eq 4 or valuetype eq 5 then begin
  newvalue = string(extravalue, format='(f0)')  ;  float or double
endif else begin
  newvalue = strtrim(extravalue,2)
endelse
newextra = {format:newformat, name:newname, value:newvalue, comment:newcomment}

; Add the new extra tag to the loaded single-channel spectrum
; or to the specified single-channel stack spectrum
;
; Have to replace the entire single-channel spectrum structure: IDL won't
; let you change the dimensions of the existing extratags array.
;
; IDL 5.3 is a pain about concatenating anonymous structures, so it has to be
; done in a tedious fashion.

if n_elements(n) eq 0 then begin
  oldStruc = *loaded1
  extratags = reform(oldStruc.extratags)
  nextra = oldStruc.nextra
  extratags = [extratags, extratags[0]]
  for k=0L,3 do extratags[nextra].(k) = newextra.(k)
  nextra = nextra + 1L
  newStruc = {freq:oldStruc.freq, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *loaded1 = newStruc
endif else begin
  oldStruc = *stack1[n-1]
  extratags = reform(oldStruc.extratags)
  nextra = oldStruc.nextra
  extratags = [extratags, extratags[0]]
  for k=0L,3 do extratags[nextra].(k) = newextra.(k)
  nextra = nextra + 1L
  newStruc = {group:oldStruc.group, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *stack1[n-1] = newStruc
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro addextrai,extraformat,extraname,extravalue, $
              comment=extracomment,stack=n,help=help

; Add an extra tag to the loaded image, or else to
; stacki image n if the stack keyword is set
;
; Floating-point and double-precision values are formatted with the highest
; required precision (format='(f0)').  If some other precision is desired,
; use the "string" function to create an appropriately formatted string,
; then pass that string as the extravalue argument.
;
; May 2002: Added the comment keyword for specifying the comment field

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(extraformat) eq 0 then extraformat = ''  ; just so it's defined
if n_elements(extraname) eq 0 then extraname = ''
if n_elements(extracomment) eq 0 then extracomment = ''

if keyword_set(help) or n_params() ne 3 or size(extraformat, /type) ne 7 $
                     or size(extraname, /type) ne 7 $
                     or size(extracomment, /type) ne 7 then begin
  print,' '
  print,'addextrai,format,name,value[,comment=comment][,stack=n][,/help]'
  print,'          format, name, and comment must be quoted strings'
  print,'          stack = n adds an extra tag to stacki image #n rather ', $
        'than to the loaded image',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in addextrai: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in addextrai: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in addextrai: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in addextrai: No image is loaded'
  return
endif else if isnull(extraname) and notnull(extraformat) then begin
  print,'ERROR in addextrai: Must specify a tag name'
  return
endif else if isnull(extraname) and isnull(extraformat) and $
              isnull(extravalue) and isnull(extracomment) then begin
  print,"ERROR in addextrai: Can't specify a blank comment line"
  return
endif

; Construct a new structure containing the extra tag

newcomment = strtrim(extracomment,2)
if strmid(newcomment,0,1) eq '#' then newcomment = strtrim(strmid(newcomment,1), 2)
newcomment = (isnull(newcomment)) ? '' : '# ' + newcomment
newformat = strtrim(extraformat,2)
newname = strtrim(extraname,2)
valuetype = size(extravalue, /type)
if valuetype eq 4 or valuetype eq 5 then begin
  newvalue = string(extravalue, format='(f0)')  ;  float or double
endif else begin
  newvalue = strtrim(extravalue,2)
endelse
newextra = {format:newformat, name:newname, value:newvalue, comment:newcomment}

; Add the new extra tag to the loaded image
; or to the specified stacki image
;
; Have to replace the entire image structure: IDL won't
; let you change the dimensions of the existing extratags array.
;
; IDL 5.3 is a pain about concatenating anonymous structures, so it has to be
; done in a tedious fashion.

if n_elements(n) eq 0 then begin
  oldStruc = *loadedi
  extratags = reform(oldStruc.extratags)
  nextra = oldStruc.nextra
  extratags = [extratags, extratags[0]]
  for k=0L,3 do extratags[nextra].(k) = newextra.(k)
  nextra = nextra + 1L
  newStruc = {image:oldStruc.image, extratags:extratags, $
              width:oldStruc.width, height:oldStruc.height, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *loadedi = newStruc
endif else begin
  oldStruc = *stacki[n-1]
  extratags = reform(oldStruc.extratags)
  nextra = oldStruc.nextra
  extratags = [extratags, extratags[0]]
  for k=0L,3 do extratags[nextra].(k) = newextra.(k)
  nextra = nextra + 1L
  newStruc = {group:oldStruc.group, image:oldStruc.image, extratags:extratags, $
              width:oldStruc.width, height:oldStruc.height, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *stacki[n-1] = newStruc
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro addcomment,comment,chanflags,stack=n,help=help

; Add a comment line to the extra tags of the loaded pair, or else of
; stack pair n if the stack keyword is set
;
; chanflags is an optional two-element vector containing 0 and/or 1.  If only one
; channel is involved (chanflags ne [1,1]) then an appropriate suffix is added
; to the start of the comment.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if size(comment, /type) ne 7 then comment = ''
if n_elements(chanflags) ne 2 then chanflags = [0,0]
if n_params() le 1 then chanflags = [1,1]

if keyword_set(help) or n_params() lt 1 or n_params() gt 2 or isnull(comment) $
             or (total(chanflags) lt 1) then begin
  print,' '
  print,'addcomment,comment[,chanflags][,stack=n][,/help]'
  print,'         comment must be a quoted non-null string'
  print,'         chanflags = [1,0] if tag is relevant to OC, [0,1] if SC, ', $
        ' [1,1] if both (default), [0,0] if neither (no tag added)',format='(2a)'
  print,'         stack = n adds a comment line to the extra tags of stack ', $
        'pair #n rather than of the loaded pair',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in addcomment: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in addcomment: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in addcomment: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in addcomment: No pair is loaded'
  return
endif else if total(chanflags) eq 0 then begin
  return  ; no extra tags to set!
endif

if n_elements(n) eq 0 then begin
  addextra,'','','',chanflags,comment=comment
endif else begin
  addextra,'','','',chanflags,comment=comment,stack=n
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro addcomment1,comment,stack=n,help=help

; Add a comment line to the extra tags of the loaded single-channel spectrum,
; or else of single-channel stack spectrum n if the stack keyword is set

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if size(comment, /type) ne 7 then comment = ''

if keyword_set(help) or n_params() ne 1 or isnull(comment) then begin
  print,' '
  print,'addcomment1,comment[,stack=n][,/help]'
  print,'         comment must be a quoted non-null string'
  print,'         stack = n adds an extra tag to stack1 spectrum #n rather ', $
        'than to the loaded single-channel spectrum',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in addcomment1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in addcomment1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in addcomment1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in addcomment1: No single-channel spectrum is loaded'
  return
endif

if n_elements(n) eq 0 then begin
  addextra1,'','','',comment=comment
endif else begin
  addextra1,'','','',comment=comment,stack=n
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro addcommenti,comment,stack=n,help=help

; Add a comment line to the extra tags of the loaded image,
; or else of stacki image n if the stack keyword is set

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if size(comment, /type) ne 7 then comment = ''

if keyword_set(help) or n_params() ne 1 or isnull(comment) then begin
  print,' '
  print,'addcommenti,comment[,stack=n][,/help]'
  print,'         comment must be a quoted non-null string'
  print,'         stack = n adds an extra tag to stacki image #n rather ', $
        'than to the loaded image',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in addcommenti: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in addcommenti: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in addcommenti: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in addcommenti: No image is loaded'
  return
endif

if n_elements(n) eq 0 then begin
  addextrai,'','','',comment=comment
endif else begin
  addextrai,'','','',comment=comment,stack=n
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro deleteextra,name,stack=n,silent=silent,help=help

; Delete an extra tag on the loaded pair, or else on stack pair n if the 
; stack keyword is set
;
; 2002 May: Now allow "name" to be an integer, the number of the extra tag
;           to be deleted; this will be especially useful for deleting comments
;           (for which extratags.name is the null string)

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

nametype = size(name, /type)
if nametype ne 2 and nametype ne 3 and nametype ne 7 then name = ''  ;  just so it's defined as a string

if keyword_set(help) or isnull(name) or n_params() ne 1 then begin
  print,' '
  print,'deleteextra,name[,stack=n][,/silent][,/help]'
  print,'            name is a quoted string OR the number (1-based)', $
                ' of the extra tag to be deleted',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in deleteextra: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in deleteextra: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in deleteextra: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in deleteextra: No pair is loaded'
  return
endif

; Get the array of extra tags

oldStruc = (n_elements(n) eq 0) ? *loaded : *stack[n-1]
extratags = reform(oldStruc.extratags)
nextra = oldStruc.nextra

if nextra eq 1 then begin
  print,'ERROR in deleteextra: Must leave at least one extra tag intact'
  return
endif else if nametype eq 2 or nametype eq 3 then begin
  if name lt 1 or name gt nextra then begin
    print,'ERROR in deleteextra: Extra tag number must be in the range [1, ',nextra,']', $
          format='(a,i0,a)'
    return
  endif
endif

; Search and destroy
;
; Have to replace the entire pair structure: IDL won't let you change
; the dimensions of the existing extratags array.

if nametype eq 7 then begin
  extranum = where(strlowcase(extratags.name) eq strlowcase(name), count)
  if count eq 0 then begin
    if not keyword_set(silent) then $
      print,"deleteextra: Extra tag '",name,"' is not present"
    return
  endif
  extranum = extranum[0]
endif else begin
  extranum = name - 1
  if not keyword_set(silent) then begin
    answer = ''  ;  define as a string for read
    print,'Extra tag #',extranum+1,' out of ',nextra,':',format='(a,i0,a,i0,a)'
    sformat = '(a,3x,a16,3x,a,3x,a)'
    print,strtrim(string(extratags[extranum], format=sformat), 2)
    print,'Delete this extra tag (y/n)? ',format='(a,$)'
    read,answer
    if strmid(strlowcase(answer),0,1) ne 'y' then begin
      print,'deleteextra: Operation canceled'
      return
    endif
  endif
endelse

; Tag is present, so delete it

if extranum eq 0 then begin
  extratags = extratags[1:nextra-1]
endif else if extranum eq nextra-1 then begin
  extratags = extratags[0:nextra-2]
endif else begin
  extratags = [extratags[0:extranum-1], extratags[extranum+1:nextra-1]]
endelse
nextra = nextra - 1L
if n_elements(n) eq 0 then begin
  newStruc = {freq:oldStruc.freq, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname}
  *loaded = newStruc
endif else begin
  newStruc = {group:oldStruc.group, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname}
  *stack[n-1] = newStruc
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro deleteextra1,name,stack=n,silent=silent,help=help

; Delete an extra tag on the loaded single-channel spectrum, or else on
; single-channel stack spectrum n if the stack keyword is set
;
; 2002 May: Now allow "name" to be an integer, the number of the extra tag
;           to be deleted; this will be especially useful for deleting comments
;           (for which extratags.name is the null string)

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

nametype = size(name, /type)
if nametype ne 2 and nametype ne 3 and nametype ne 7 then name = ''  ;  just so it's defined as a string

if keyword_set(help) or isnull(name) or n_params() ne 1 then begin
  print,' '
  print,'deleteextra1,name[,stack=n][,/silent][,/help]'
  print,'            name is a quoted string OR the number (1-based)', $
                ' of the extra tag to be deleted',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in deleteextra1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in deleteextra1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in deleteextra1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in deleteextra1: No single-channel spectrum is loaded'
  return
endif

; Get the array of extra tags

oldStruc = (n_elements(n) eq 0) ? *loaded1 : *stack1[n-1]
extratags = reform(oldStruc.extratags)
nextra = oldStruc.nextra

if nextra eq 1 then begin
  print,'ERROR in deleteextra1: Must leave at least one extra tag intact'
  return
endif else if nametype eq 2 or nametype eq 3 then begin
  if name lt 1 or name gt nextra then begin
    print,'ERROR in deleteextra1: Extra tag number must be in the range [1, ',nextra,']', $
          format='(a,i0,a)'
    return
  endif
endif

; Search and destroy
;
; Have to replace the entire single-channel spectrum structure: IDL won't
; let you change the dimensions of the existing extratags array.

if nametype eq 7 then begin
  extranum = where(strlowcase(extratags.name) eq strlowcase(name), count)
  if count eq 0 then begin
    if not keyword_set(silent) then $
      print,"deleteextra1: Extra tag '",name,"' is not present"
    return
  endif
  extranum = extranum[0]
endif else begin
  extranum = name - 1
  if not keyword_set(silent) then begin
    answer = ''  ;  define as a string for read
    print,'Extra tag #',extranum+1,' out of ',nextra,':',format='(a,i0,a,i0,a)'
    sformat = '(a,3x,a16,3x,a,3x,a)'
    print,strtrim(string(extratags[extranum], format=sformat), 2)
    print,'Delete this extra tag (y/n)? ',format='(a,$)'
    read,answer
    if strmid(strlowcase(answer),0,1) ne 'y' then begin
      print,'deleteextra1: Operation canceled'
      return
    endif
  endif
endelse

; Tag is present, so delete it

if extranum eq 0 then begin
  extratags = extratags[1:nextra-1]
endif else if extranum eq nextra-1 then begin
  extratags = extratags[0:nextra-2]
endif else begin
  extratags = [extratags[0:extranum-1], extratags[extranum+1:nextra-1]]
endelse
nextra = nextra - 1
if n_elements(n) eq 0 then begin
  newStruc = {freq:oldStruc.freq, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *loaded1 = newStruc
endif else begin
  newStruc = {group:oldStruc.group, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *stack1[n-1] = newStruc
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro deleteextrai,name,stack=n,silent=silent,help=help

; Delete an extra tag on the loaded image, or else on
; stacki image n if the stack keyword is set
;
; 2002 May: Now allow "name" to be an integer, the number of the extra tag
;           to be deleted; this will be especially useful for deleting comments
;           (for which extratags.name is the null string)

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

nametype = size(name, /type)
if nametype ne 2 and nametype ne 3 and nametype ne 7 then name = ''  ;  just so it's defined as a string

if keyword_set(help) or isnull(name) or n_params() ne 1 then begin
  print,' '
  print,'deleteextrai,name[,stack=n][,/silent][,/help]'
  print,'             name must be a quoted string'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in deleteextrai: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in deleteextrai: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in deleteextrai: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in deleteextrai: No image is loaded'
  return
endif

; Get the array of extra tags

oldStruc = (n_elements(n) eq 0) ? *loadedi : *stacki[n-1]
extratags = reform(oldStruc.extratags)
nextra = oldStruc.nextra

if nextra eq 1 then begin
  print,'ERROR in deleteextrai: Must leave at least one extra tag intact'
  return
endif else if nametype eq 2 or nametype eq 3 then begin
  if name lt 1 or name gt nextra then begin
    print,'ERROR in deleteextrai: Extra tag number must be in the range [1, ',nextra,']', $
          format='(a,i0,a)'
    return
  endif
endif

; Search and destroy
;
; Have to replace the entire image structure: IDL won't
; let you change the dimensions of the existing extratags array.

if nametype eq 7 then begin
  extranum = where(strlowcase(extratags.name) eq strlowcase(name), count)
  if count eq 0 then begin
    if not keyword_set(silent) then $
      print,"deleteextrai: Extra tag '",name,"' is not present"
    return
  endif
  extranum = extranum[0]
endif else begin
  extranum = name - 1
  if not keyword_set(silent) then begin
    answer = ''  ;  define as a string for read
    print,'Extra tag #',extranum+1,' out of ',nextra,':',format='(a,i0,a,i0,a)'
    sformat = '(a,3x,a16,3x,a,3x,a)'
    print,strtrim(string(extratags[extranum], format=sformat), 2)
    print,'Delete this extra tag (y/n)? ',format='(a,$)'
    read,answer
    if strmid(strlowcase(answer),0,1) ne 'y' then begin
      print,'deleteextrai: Operation canceled'
      return
    endif
  endif
endelse

; Tag is present, so delete it

if extranum eq 0 then begin
  extratags = extratags[1:nextra-1]
endif else if extranum eq nextra-1 then begin
  extratags = extratags[0:nextra-2]
endif else begin
  extratags = [extratags[0:extranum-1], extratags[extranum+1:nextra-1]]
endelse
nextra = nextra - 1
if n_elements(n) eq 0 then begin
  newStruc = {image:oldStruc.image, extratags:extratags, $
              width:oldStruc.width, height:oldStruc.height, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *loadedi = newStruc
endif else begin
  newStruc = {group:oldStruc.group, image:oldStruc.image, extratags:extratags, $
              width:oldStruc.width, height:oldStruc.height, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *stacki[n-1] = newStruc
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setextra,extraformat,extraname,extravalue,comment=extracomment,stack=n,help=help

; (Re)sets the value and/or format of an extra tag on the loaded pair,
; or else on stack pair n if the stack keyword is set.
;
; If no extra tag by the specified name is found, a new extra tag is added.
;
; Floating-point and double-precision values are formatted with the highest
; required precision (format='(f0)').  If some other precision is desired,
; use the "string" function to create an appropriately formatted string,
; then pass that string as the extravalue argument.
;
; May 2002: Added the comment keyword for specifying (or erasing) the comment field

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(extraformat) eq 0 then extraformat = ''  ; just so it's defined
if n_elements(extraname) eq 0 then extraname = ''
extracomment_specified = n_elements(extracomment) gt 0
if not extracomment_specified then extracomment = ''

if keyword_set(help) or n_params() ne 3 or size(extraformat, /type) ne 7  $
                     or size(extraname, /type) ne 7  $
                     or size(extracomment, /type) ne 7 then begin
  print,' '
  print,'setextra,format,name,value[,comment=comment][,stack=n][,/help]'
  print,'         format, name, and comment must be quoted strings'
  print,'         stack=n (re)sets an extra tag to stack pair #n rather than ', $
        'to the loaded pair',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in setextra: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in setextra: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in setextra: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in setextra: No pair is loaded'
  return
endif

; Construct a new structure containing the extra tag

newformat = strtrim(extraformat,2)
newname = strtrim(extraname,2)
valuetype = size(extravalue, /type)
if valuetype eq 4 or valuetype eq 5 then begin
  newvalue = string(extravalue, format='(f0)')  ;  float or double
endif else begin
  newvalue = strtrim(extravalue,2)
endelse
newcomment = strtrim(extracomment,2)
if strmid(newcomment,0,1) eq '#' then newcomment = strtrim(strmid(newcomment,1), 2)
newcomment = (isnull(newcomment)) ? '' : '# ' + newcomment
newextra = {format:newformat, name:newname, value:newvalue, comment:newcomment}

; Get the array of extra tags

oldStruc = (n_elements(n) gt 0) ? *stack[n-1] : *loaded
extratags = reform(oldStruc.extratags)  ;  get rid of useless extra dimension
nextra = oldStruc.nextra

; Search for an existing tag with that name

extranum = where(strlowcase(extratags.name) eq strlowcase(newname), count)

if count gt 0 then begin

  ; Tag present: change format and/or value
  ;              change comment IF a new one (even '') was explicitly given

  extranum = extranum[0]
  for k=0L,2 do extratags[extranum].(k) = newextra.(k)
  if extracomment_specified then extratags[extranum].(3) = newextra.(3)
  
endif else begin

  ; Tag not present: add a new extra tag

  extratags = [extratags, extratags[0]]
  for k=0L,3 do extratags[nextra].(k) = newextra.(k)
  nextra = nextra + 1L

endelse

; Create a new pair structure and load it or put it back in the stack
;
; Have to replace the entire pair structure: If you've added a tag, IDL
; won't let you change the dimensions of the existing extratags array

if n_elements(n) eq 0 then begin
  newStruc = {freq:oldStruc.freq, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname}
  *loaded = newStruc
endif else begin
  newStruc = {group:oldStruc.group, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname}
  *stack[n-1] = newStruc
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setextra1,extraformat,extraname,extravalue,comment=extracomment,stack=n,help=help

; (Re)sets the value and/or format of an extra tag on the loaded
; single-channel spectrum, or else on stack1 spectrum n if the stack
; keyword is set.
;
; If no extra tag by the specified name is found, a new extra tag is added.
;
; Floating-point and double-precision values are formatted with the highest
; required precision (format='(f0)').  If some other precision is desired,
; use the "string" function to create an appropriately formatted string,
; then pass that string as the extravalue argument.
;
; May 2002: Added the comment keyword for specifying (or erasing) the comment field

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(extraformat) eq 0 then extraformat = ''  ; just so it's defined
if n_elements(extraname) eq 0 then extraname = ''
extracomment_specified = n_elements(extracomment) gt 0
if not extracomment_specified then extracomment = ''

if keyword_set(help) or n_params() ne 3 or size(extraformat, /type) ne 7  $
                     or size(extraname, /type) ne 7  $
                     or size(extracomment, /type) ne 7 then begin
  print,' '
  print,'setextra1,format,name,value[,comment=comment][,stack=n][,/help]'
  print,'          format, name, and comment must be quoted strings'
  print,'          stack=n adds an extra tag to stack1 spectrum #n rather than ', $
        'to the loaded single-channel spectrum',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in setextra1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in setextra1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in setextra1: There are only ',nstack1, $
          ' single-channel spectra in stack1',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in setextra1: No single-channel spectrum is loaded'
  return
endif

; Construct a new structure containing the extra tag

newformat = strtrim(extraformat,2)
newname = strtrim(extraname,2)
valuetype = size(extravalue, /type)
if valuetype eq 4 or valuetype eq 5 then begin
  newvalue = string(extravalue, format='(f0)')  ;  float or double
endif else begin
  newvalue = strtrim(extravalue,2)
endelse
newcomment = strtrim(extracomment,2)
if strmid(newcomment,0,1) eq '#' then newcomment = strtrim(strmid(newcomment,1), 2)
newcomment = (isnull(newcomment)) ? '' : '# ' + newcomment
newextra = {format:newformat, name:newname, value:newvalue, comment:newcomment}

; Get the array of extra tags

oldStruc = (n_elements(n) gt 0) ? *stack1[n-1] : *loaded1
extratags = reform(oldStruc.extratags)  ;  get rid of useless extra dimension
nextra = oldStruc.nextra

; Search for an existing tag with that name

extranum = where(strlowcase(extratags.name) eq strlowcase(newname), count)

if count gt 0 then begin

  ; Tag present: change format and/or value
  ;              change comment IF a new one (even '') was explicitly given

  extranum = extranum[0]
  for k=0L,2 do extratags[extranum].(k) = newextra.(k)
  if extracomment_specified then extratags[extranum].(3) = newextra.(3)
  
endif else begin

  ; Tag not present: add a new extra tag

  extratags = [extratags, extratags[0]]
  for k=0L,3 do extratags[nextra].(k) = newextra.(k)
  nextra = nextra + 1L

endelse

; Create a new single-channel structure and load it or put it back in stack1
;
; Have to replace the entire single-channel structure: If you've added a tag,
; IDL won't let you change the dimensions of the existing extratags array

if n_elements(n) eq 0 then begin
  newStruc = {freq:oldStruc.freq, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *loaded1 = newStruc
endif else begin
  newStruc = {group:oldStruc.group, spec:oldStruc.spec, tags:oldStruc.tags, $
              extratags:extratags, ndata:oldStruc.ndata, ntags:oldStruc.ntags, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *stack1[n-1] = newStruc
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setextrai,extraformat,extraname,extravalue,comment=extracomment,stack=n,help=help

; (Re)sets the value and/or format of an extra tag on the loaded
; image, or else on stacki image n if the stack keyword is set.
;
; If no extra tag by the specified name is found, a new extra tag is added.
;
; Floating-point and double-precision values are formatted with the highest
; required precision (format='(f0)').  If some other precision is desired,
; use the "string" function to create an appropriately formatted string,
; then pass that string as the extravalue argument.
;
; May 2002: Added the comment keyword for specifying (or erasing) the comment field

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(extraformat) eq 0 then extraformat = ''  ; just so it's defined
if n_elements(extraname) eq 0 then extraname = ''
extracomment_specified =  n_elements(extracomment) gt 0
if not extracomment_specified then extracomment = ''

if keyword_set(help) or n_params() ne 3 or size(extraformat, /type) ne 7  $
                     or size(extraname, /type) ne 7  $
                     or size(extracomment, /type) ne 7 then begin
  print,' '
  print,'setextrai,format,name,value[,comment=comment][,stack=n][,/help]'
  print,'          format, name, and comment must be quoted strings'
  print,'          stack=n adds an extra tag to stacki image #n rather than ', $
        'to the loaded image',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in setextrai: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in setextrai: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in setextrai: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in setextrai: No image is loaded'
  return
endif

; Construct a new structure containing the extra tag

newformat = strtrim(extraformat,2)
newname = strtrim(extraname,2)
valuetype = size(extravalue, /type)
if valuetype eq 4 or valuetype eq 5 then begin
  newvalue = string(extravalue, format='(f0)')  ;  float or double
endif else begin
  newvalue = strtrim(extravalue,2)
endelse
newcomment = strtrim(extracomment,2)
if strmid(newcomment,0,1) eq '#' then newcomment = strtrim(strmid(newcomment,1), 2)
newcomment = (isnull(newcomment)) ? '' : '# ' + newcomment
newextra = {format:newformat, name:newname, value:newvalue, comment:newcomment}

; Get the array of extra tags

oldStruc = (n_elements(n) gt 0) ? *stacki[n-1] : *loadedi
extratags = reform(oldStruc.extratags)  ;  get rid of useless extra dimension
nextra = oldStruc.nextra

; Search for an existing tag with that name

extranum = where(strlowcase(extratags.name) eq strlowcase(newname), count)

if count gt 0 then begin

  ; Tag present: change format and/or value
  ;              change comment IF a new one (even '') was explicitly given

  extranum = extranum[0]
  for k=0L,2 do extratags[extranum].(k) = newextra.(k)
  if extracomment_specified then extratags[extranum].(3) = newextra.(3)
  
endif else begin

  ; Tag not present: add a new extra tag

  extratags = [extratags, extratags[0]]
  for k=0L,3 do extratags[nextra].(k) = newextra.(k)
  nextra = nextra + 1L

endelse

; Create a new image structure and load it or put it back in stacki
;
; Have to replace the entire image structure: If you've added a tag,
; IDL won't let you change the dimensions of the existing extratags array

if n_elements(n) eq 0 then begin
  newStruc = {image:oldStruc.image, extratags:extratags, $
              width:oldStruc.width, height:oldStruc.height, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *loadedi = newStruc
endif else begin
  newStruc = {group:oldStruc.group, image:oldStruc.image, extratags:extratags, $
              width:oldStruc.width, height:oldStruc.height, $
              nextra:nextra, tname:oldStruc.tname, pol:oldStruc.pol}
  *stacki[n-1] = newStruc
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setextraname,oldname,newname,stack=n,silent=silent,help=help

; Resets the name of an extra tag on the loaded pair,
; or else on stack pair n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(oldname) eq 0 then oldname = ''  ; just so it's defined
if n_elements(newname) eq 0 then newname = ''

if keyword_set(help) or n_params() ne 2 or size(oldname, /type) ne 7  $
                     or size(newname, /type) ne 7 then begin
  print,' '
  print,'setextraname,oldname,newname[,stack=n][,/silent][,/help]'
  print,'         oldname and newname must be quoted strings'
  print,'         stack=n changes the name of an extra tag to stack pair #n ', $
        'rather than to the loaded pair',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in setextraname: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in setextraname: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in setextraname: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in setextraname: No pair is loaded'
  return
endif

; Get the array of extra tags and search for the one whose name is oldname

extratags = (n_elements(n) gt 0) ? (*stack[n-1]).extratags : (*loaded).extratags
extranum = where(strlowcase(extratags.name) eq strlowcase(oldname), count)
if count eq 0 then begin
  if not keyword_set(silent) then $
    print,'setextraname: Extra tag ',oldname,' is not present'
  return
endif

; Extra tag is present: change its name to newname

extranum = extranum[0]
if n_elements(n) eq 0 then $
  (*loaded).extratags[extranum].name = strtrim(newname,2) $
else $
  (*stack[n-1]).extratags[extranum].name = strtrim(newname,2)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setextraname1,oldname,newname,stack=n,silent=silent,help=help

; Resets the name of an extra tag on the loaded single-channel spectrum,
; or else on stack1 spectrum n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(oldname) eq 0 then oldname = ''  ; just so it's defined
if n_elements(newname) eq 0 then newname = ''

if keyword_set(help) or n_params() ne 2 or size(oldname, /type) ne 7  $
                     or size(newname, /type) ne 7 then begin
  print,' '
  print,'setextraname1,oldname,newname[,stack=n][,/silent][,/help]'
  print,'         oldname and newname must be quoted strings'
  print,'         stack=n changes the name of an extra tag to stack1 spectrum #n ', $
        'rather than to the loaded single-channel spectrum',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in setextraname1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in setextraname1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in setextraname1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in setextraname1: No single-channel spectrum is loaded'
  return
endif

; Get the array of extra tags and search for the one whose name is oldname

extratags = (n_elements(n) gt 0) ? (*stack1[n-1]).extratags : (*loaded1).extratags
extranum = where(strlowcase(extratags.name) eq strlowcase(oldname), count)
if count eq 0 then begin
  if not keyword_set(silent) then $
    print,'setextraname1: Extra tag ',oldname,' is not present'
  return
endif

; Extra tag is present: change its name to newname

extranum = extranum[0]
if n_elements(n) eq 0 then $
  (*loaded1).extratags[extranum].name = strtrim(newname,2) $
else $
  (*stack1[n-1]).extratags[extranum].name = strtrim(newname,2)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setextranamei,oldname,newname,stack=n,silent=silent,help=help

; Resets the name of an extra tag on the loaded image,
; or else on stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(oldname) eq 0 then oldname = ''  ; just so it's defined
if n_elements(newname) eq 0 then newname = ''

if keyword_set(help) or n_params() ne 2 or size(oldname, /type) ne 7  $
                     or size(newname, /type) ne 7 then begin
  print,' '
  print,'setextranamei,oldname,newname[,stack=n][,/silent][,/help]'
  print,'         oldname and newname must be quoted strings'
  print,'         stack=n changes the name of an extra tag to stacki image #n ', $
        'rather than to the loaded image',format='(2a)'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in setextranamei: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in setextranamei: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in setextranamei: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in setextranamei: No image is loaded'
  return
endif

; Get the array of extra tags and search for the one whose name is oldname

extratags = (n_elements(n) gt 0) ? (*stacki[n-1]).extratags : (*loadedi).extratags
extranum = where(strlowcase(extratags.name) eq strlowcase(oldname), count)
if count eq 0 then begin
  if not keyword_set(silent) then $
    print,'setextranamei: Extra tag ',oldname,' is not present'
  return
endif

; Extra tag is present: change its name to newname

extranum = extranum[0]
if n_elements(n) eq 0 then $
  (*loadedi).extratags[extranum].name = strtrim(newname,2) $
else $
  (*stacki[n-1]).extratags[extranum].name = strtrim(newname,2)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showextra,extraname,stack=n,help=help

; Displays an extra tag on the loaded pair, or else on
; stack pair n if the stack keyword is set.
;
; If extraname is omitted (or set to 'all') then all extra tags are displayed.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then extraname = 'all'

if keyword_set(help) or n_params() gt 1 $
                     or size(extraname, /type) ne 7 then begin
  print,' '
  print,'showextra[,name][,stack=n][,/help])'
  print,'          name must be a quoted string'
  print,' '
  print,'If no extra tag name is given, all extra (image) tags are displayed'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in showextra: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showextra: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in showextra: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in showextra: No pair is loaded'
  return
endif

; Get the extra tags to be displayed

if n_elements(n) eq 0 then begin
  extratags = (*loaded).extratags
  nextra = (*loaded).nextra
endif else begin
  extratags = (*stack[n-1]).extratags
  nextra = (*stack[n-1]).nextra
endelse
name = strlowcase(strtrim(extraname,2))
sformat = '(a,3x,a16,3x,a,3x,a)'
if name eq 'all' then begin
  displayextra = strtrim(string(extratags, format=sformat), 2)
  ndisplay = nextra
endif else begin

  ; Search for the desired extra tag

  extranum = where(strlowcase(extratags.name) eq name, count)

  if count gt 0 then begin

    ; Tag present

    extranum = extranum[0]
    displayextra = [strtrim(string(extratags[extranum], format=sformat), 2)]
    ndisplay = 1L
  
  endif else begin

    ; Tag not present

    print,"showextra: Extra tag '",name,"' isn't present"
    return

  endelse
endelse

; Do the display, three extra tags to a line
; (format each column separately)

imax = 2 < (ndisplay-1)
for i=0L,imax-1 do begin
  indmin = i*(ndisplay/3) + (i < (ndisplay mod 3))
  indmax = (i+1)*(ndisplay/3) + ((i+1) < (ndisplay mod 3)) - 1
  maxlen = max(strlen(displayextra[indmin:indmax]))
  newlen = 33 > (maxlen + 6)
  for k=indmin,indmax do begin
    nblanks = newlen - strlen(displayextra[k])
    for j=0L,nblanks-1 do displayextra[k] = displayextra[k] + ' '
  endfor
endfor
nlines = ceil(ndisplay/3.0)
print,' '
for k=0L,nlines-1 do begin
  printstring = displayextra[k]
  if k+nlines lt ndisplay then printstring = printstring + displayextra[k+nlines]
  if k+2*nlines lt ndisplay then printstring = printstring + displayextra[k+2*nlines]
  print,printstring
endfor
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showextra1,extraname,stack=n,help=help

; Displays an extra tag on the loaded single-channel spectrum, or else on
; stack1 spectrum n if the stack keyword is set.
;
; If extraname is omitted (or set to 'all') then all extra tags are displayed.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then extraname = 'all'

if keyword_set(help) or n_params() gt 1 $
                     or size(extraname, /type) ne 7 then begin
  print,' '
  print,'showextra1[,name][,stack=n][,/help])'
  print,'           name must be a quoted string'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in showextra1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showextra1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in showextra1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in showextra1: No single-channel spectrum is loaded'
  return
endif

; Get the extra tags to be displayed

if n_elements(n) eq 0 then begin
  extratags = (*loaded1).extratags
  nextra = (*loaded1).nextra
endif else begin
  extratags = (*stack1[n-1]).extratags
  nextra = (*stack1[n-1]).nextra
endelse
name = strlowcase(strtrim(extraname,2))
sformat = '(a,3x,a16,3x,a,3x,a)'
if name eq 'all' then begin
  displayextra = strtrim(string(extratags, format=sformat), 2)
  ndisplay = nextra
endif else begin

  ; Search for the desired extra tag

  extranum = where(strlowcase(extratags.name) eq name, count)

  if count gt 0 then begin

    ; Tag present

    extranum = extranum[0]
    displayextra = [strtrim(string(extratags[extranum], format=sformat), 2)]
    ndisplay = 1L
  
  endif else begin

    ; Tag not present

    print,"showextra1: Extra tag '",name,"' isn't present"
    return

  endelse
endelse

; Do the display, three extra tags to a line
; (format each column separately)

imax = 2 < (ndisplay-1)
for i=0L,imax-1 do begin
  indmin = i*(ndisplay/3) + (i < (ndisplay mod 3))
  indmax = (i+1)*(ndisplay/3) + ((i+1) < (ndisplay mod 3)) - 1
  maxlen = max(strlen(displayextra[indmin:indmax]))
  newlen = 33 > (maxlen + 6)
  for k=indmin,indmax do begin
    nblanks = newlen - strlen(displayextra[k])
    for j=0L,nblanks-1 do displayextra[k] = displayextra[k] + ' '
  endfor
endfor
nlines = ceil(ndisplay/3.0)
print,' '
for k=0L,nlines-1 do begin
  printstring = displayextra[k]
  if k+nlines lt ndisplay then printstring = printstring + displayextra[k+nlines]
  if k+2*nlines lt ndisplay then printstring = printstring + displayextra[k+2*nlines]
  print,printstring
endfor
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showextrai,extraname,stack=n,help=help

; Displays an extra tag on the loaded image, or else on
; stacki image n if the stack keyword is set.
;
; If extraname is omitted (or set to 'all') then all extra tags are displayed.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then extraname = 'all'

if keyword_set(help) or n_params() gt 1 $
                     or size(extraname, /type) ne 7 then begin
  print,' '
  print,'showextrai[,name][,stack=n][,/help])'
  print,'           name must be a quoted string'
  print,' '
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in showextrai: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showextrai: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in showextrai: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in showextrai: No image is loaded'
  return
endif

; Get the extra tags to be displayed

if n_elements(n) eq 0 then begin
  extratags = (*loadedi).extratags
  nextra = (*loadedi).nextra
endif else begin
  extratags = (*stacki[n-1]).extratags
  nextra = (*stacki[n-1]).nextra
endelse
name = strlowcase(strtrim(extraname,2))
sformat = '(a,3x,a16,3x,a,3x,a)'
if name eq 'all' then begin
  displayextra = strtrim(string(extratags, format=sformat), 2)
  ndisplay = nextra
endif else begin

  ; Search for the desired extra tag

  extranum = where(strlowcase(extratags.name) eq name, count)

  if count gt 0 then begin

    ; Tag present

    extranum = extranum[0]
    displayextra = [strtrim(string(extratags[extranum], format=sformat), 2)]
    ndisplay = 1L
  
  endif else begin

    ; Tag not present

    print,"showextrai: Extra tag '",name,"' isn't present"
    return

  endelse
endelse

; Do the display, one extra tag to a line

print,' '
for k=0L,ndisplay-1 do print,displayextra[k]
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showtarget,stack=n,help=help

; Displays the target name for the loaded pair, or else for stack pair n
; if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showtarget[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in showtarget: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showtarget: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in showtarget: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in showtarget: No pair is loaded'
  return
endif

; Get the target name to be displayed

tname = (n_elements(n) gt 0) ? (*stack[n-1]).tname : (*loaded).tname
print,' '
print,'Target: ',tname
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showtarget1,stack=n,help=help

; Displays the target name for the loaded single-channel spectrum,
; or else for stack1 spectrum n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showtarget1[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in showtarget1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showtarget1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in showtarget1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in showtarget1: No single-channel spectrum is loaded'
  return
endif

; Get the target name to be displayed

tname = (n_elements(n) gt 0) ? (*stack1[n-1]).tname : (*loaded1).tname
print,' '
print,'Target: ',tname
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showtargeti,stack=n,help=help

; Displays the target name for the loaded image,
; or else for stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showtargeti[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in showtargeti: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showtargeti: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in showtargeti: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in showtargeti: No image is loaded'
  return
endif

; Get the target name to be displayed

tname = (n_elements(n) gt 0) ? (*stacki[n-1]).tname : (*loadedi).tname
print,' '
print,'Target: ',tname
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro settarget,tname,stack=n,help=help

; Changes the target name for the loaded pair, or else for stack pair n
; if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(tname) eq 0 then tname = ''  ;  just so it's defined

if keyword_set(help) or n_params() ne 1 $
                     or size(tname, /type) ne 7 then begin
  print,'settarget,name[,stack=n][,/help]'
  print,'          name must be a quoted string'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in settarget: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in settarget: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in settarget: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in settarget: No pair is loaded'
  return
endif

; Change the target name

if n_elements(n) eq 0 then $
  (*loaded).tname = tname $
else $
  (*stack[n-1]).tname = tname

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro settarget1,tname,stack=n,help=help

; Changes the target name for the loaded single-channel spectrum,
; or else for stack1 spectrum n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(tname) eq 0 then tname = ''  ;  just so it's defined

if keyword_set(help) or n_params() ne 1 $
                     or size(tname, /type) ne 7 then begin
  print,'settarget1,name[,stack=n][,/help]'
  print,'           name must be a quoted string'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in settarget1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in settarget1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in settarget1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in settarget1: No single-channel spectrum is loaded'
  return
endif

; Change the target name

if n_elements(n) eq 0 then $
  (*loaded1).tname = tname $
else $
  (*stack1[n-1]).tname = tname

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro settargeti,tname,stack=n,help=help

; Changes the target name for the loaded image,
; or else for stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_elements(tname) eq 0 then tname = ''  ;  just so it's defined

if keyword_set(help) or n_params() ne 1 $
                     or size(tname, /type) ne 7 then begin
  print,'settargeti,name[,stack=n][,/help])'
  print,'           name must be a quoted string'
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in settargeti: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in settargeti: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in settargeti: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in settargeti: No image is loaded'
  return
endif

; Change the target name

if n_elements(n) eq 0 then $
  (*loadedi).tname = tname $
else $
  (*stacki[n-1]).tname = tname

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showpol1,stack=n,help=help

; Displays the polarization for the loaded single-channel spectrum,
; or else for stack1 spectrum n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showpol1[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in showpol1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showpol1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in showpol1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in showpol1: No single-channel spectrum is loaded'
  return
endif

; Get the polarization string to be displayed

pol = (n_elements(n) gt 0) ? (*stack1[n-1]).pol : (*loaded1).pol
print,' '
print,'Polarization: ',pol
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showpoli,stack=n,help=help

; Displays the polarization for the loaded image,
; or else for stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showpoli[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in showpoli: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showpoli: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in showpoli: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in showpoli: No image is loaded'
  return
endif

; Get the polarization string to be displayed

pol = (n_elements(n) gt 0) ? (*stacki[n-1]).pol : (*loadedi).pol
print,' '
print,'Polarization: ',pol
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showstruc,stack=n,help=help

; Displays information about the loaded pair structure, or else about
; the structure of stack pair n if the stack keyword is set.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showstruc[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in showstruc: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showstruc: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in showstruc: There are only ',nstack,' pairs in the stack', $
          format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in showstruc: No pair is loaded'
  return
endif

; Get the structure whose parameters are to be displayed

dispStruc = (n_elements(n) gt 0) ? *stack[n-1] : *loaded
tags = dispStruc.tags
dfreq = tags[0].dfreq
posfr = tags[0].posfr
xjcen = tags[0].xjcen
ndata = dispStruc.ndata
npol = n_elements(tags)
f = (n_elements(n) gt 0) ? posfr*dfreq*(findgen(ndata) - xjcen) : dispStruc.freq

; Display the parameters

print,' '
print,'Target:    ',dispStruc.tname,format='(2a)'
if n_elements(n) gt 0 then begin
  print,'Group:     ',dispStruc.group,format='(a,i0)'
endif else begin
  print,'Frequency: ',ndata,' points at ',dfreq,' Hz resolution  (', $
        f[0]-0.5*posfr*dfreq,' Hz to ',f[ndata-1]+0.5*posfr*dfreq,' Hz)', $
        format='(a,i0,a,f7.3,a,f9.2,a,f9.2,a)'
endelse
print,'Spectra:   ',dispStruc.ndata,' points for each of ',strtrim(string(npol),2),' channels', $
      format='(a,i0,a,a,a)'
print,'Tags:      ',dispStruc.ntags,' cw tags and ',dispStruc.nextra, $
      ' extra tags',format='(a,i0,a,i0,a)'
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showstruc1,stack=n,help=help

; Displays information about the loaded single-channel structure,
; or else about the structure of stack1 spectrum n if the stack keyword is set.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showstruc1[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in showstruc1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showstruc1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in showstruc1: There are only ',nstack1, $
          ' spectra in the single-channel stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in showstruc1: No single-channel spectrum is loaded'
  return
endif

; Get the structure whose parameters are to be displayed

dispStruc = (n_elements(n) gt 0) ? *stack1[n-1] : *loaded1
tags1 = dispStruc.tags
dfreq = tags1.dfreq
posfr = tags1.posfr
xjcen = tags1.xjcen
ndata = dispStruc.ndata
f = (n_elements(n) gt 0) ? posfr*dfreq*(findgen(ndata) - xjcen) : dispStruc.freq

; Display the parameters

print,' '
print,'Target:    ',dispStruc.tname,format='(2a)'
if n_elements(n) gt 0 then begin
  print,'Group:     ',dispStruc.group,format='(a,i0)'
endif else begin
  print,'Frequency: ',ndata,' points at ',dfreq,' Hz resolution  (', $
        f[0]-0.5*posfr*dfreq,' Hz to ',f[ndata-1]+0.5*posfr*dfreq,' Hz)', $
        format='(a,i0,a,f7.3,a,f9.2,a,f9.2,a)'
endelse
print,'Spectrum:  ',ndata,' points',format='(a,i0,a)'
print,'Tags:      ',dispStruc.ntags,' cw tags and ',dispStruc.nextra, $
      ' extra tags',format='(a,i0,a,i0,a)'
print,'Channel:   ',dispStruc.pol,format='(2a)'
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showstruci,stack=n,help=help

; Displays information about the loaded image stucture,
; or else about the structure of stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showstruci[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in showstruci: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showstruci: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in showstruci: There are only ',nstacki, $
          ' images in the image stack',format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in showstruci: No image is loaded'
  return
endif

; Get the parameters to be displayed

if n_elements(n) gt 0 then begin
  group = (*stacki[n-1]).group
  width = (*stacki[n-1]).width
  height = (*stacki[n-1]).height
  nextra = (*stacki[n-1]).nextra
  tname = (*stacki[n-1]).tname
  pol = (*stacki[n-1]).pol
  dfreq = getextrai('fres',stack=n)
  eph_col = getextrai('eph_col',stack=n)
  eph_row = getextrai('eph_row',stack=n)
  delayunit = getextrai('delayunit',stack=n)
  spb = getextrai('samples_per_baud',stack=n)
  rows_per_baud = getextrai('rows_per_baud',stack=n)
  baud = getextrai('baudlen',stack=n)
endif else begin
  width = (*loadedi).width
  height = (*loadedi).height
  nextra = (*loadedi).nextra
  tname = (*loadedi).tname
  pol = (*loadedi).pol
  dfreq = getextrai('fres')
  eph_col = getextrai('eph_col')
  eph_row = getextrai('eph_row')
  delayunit = getextrai('delayunit')
  spb = getextrai('samples_per_baud')
  rows_per_baud = getextrai('rows_per_baud')
  baud = getextrai('baudlen')
endelse
if isnull(delayunit) and notnull(baud) then begin
  if isnull(spb) then spb = 1L
  if isnull(rows_per_baud) then rows_per_baud = spb
  delayunit = baud/rows_per_baud
endif
if notnull(eph_col) and notnull(eph_row) and notnull(dfreq) and notnull(delayunit) then begin
  f = dfreq*(findgen(width) - eph_col)       ; Doppler (Hz)
  d = delayunit*(findgen(height) - eph_row)  ; delay (usec)
endif else begin
  print,"ERROR in showstruci: Need the 'eph_col' and 'eph_row' and 'fres' tags,"
  print,"                     plus either the 'delayunit' or 'baudlen' tag"
  return
endelse

; Construct the output lines for the image's delay-Doppler parameters

nchars = strlen(string((width > height), format='(i0)'))
formatstring = '(i' + string(nchars,format='(i0)') + ')'
printstring1 = 'Image:     '                                              $
               + string(width,format=formatstring) + ' Doppler cols at '  $
               + string(dfreq,format='(f7.3)') + ' Hz resolution  ('      $
               + string(f[0]-0.5*dfreq,format='(f8.2)') + ' Hz to '       $
               + string(f[width-1]+0.5*dfreq,format='(f8.2)') + ' Hz)'
printstring2 = '           '                                                $
               + string(height,format=formatstring) + ' delay   rows at '   $
               + string(delayunit,format='(f7.3)') + ' us resolution  ('    $
               + string(d[0]-0.5*delayunit,format='(f8.2)') + ' us to '     $
               + string(d[height-1]+0.5*delayunit,format='(f8.2)')          $
               + ' us)'

; Display the structure's parameters

print,' '
print,'Target:    ',tname
if n_elements(n) gt 0 then print,'Group:     ',group,format='(a,i0)'
print,printstring1
print,printstring2
print,'Tags:      ',nextra,' extra tags',format='(a,i0,a)'
print,'Channel:   ',pol
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setgroups,group=g,n=n,first=first,last=last,arr=arr,phasewidth=pw,phase0=p0, $
              jdwidth=jdw,jd0=jd0,silent=silent,help=help

; Assign group numbers to some or all pairs in the stack, either by directly
; specifying the numbers or else by grouping by phase or by time

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that a valid combination of keywords has been used

if n_elements(g) gt 0 then begin
  keywordsOK = $
    (n_elements(pw) eq 0 and n_elements(p0) eq 0  $
                         and n_elements(jdw) eq 0 and n_elements(jd0) eq 0) $
    and $
    ( (n_elements(n) gt 0 and $
       (n_elements(first) eq 0 and n_elements(last) eq 0 and n_elements(arr) eq 0)) $
      or $
      (n_elements(first) gt 0 and n_elements(last) gt 0 and $
       n_elements(n) eq 0 and n_elements(arr) eq 0) $
      or $
      (n_elements(arr) gt 0 and $
       (n_elements(n) eq 0 and n_elements(first) eq 0 and n_elements(last) eq 0)) )
endif else if n_elements(n) gt 0 or n_elements(first) gt 0 or $
              n_elements(last) gt 0 or n_elements(arr) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(pw) gt 0 then begin
  keywordsOK = n_elements(jdw) eq 0 and n_elements(jd0) eq 0
endif else if n_elements(p0) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(jdw) gt 0 then begin
  keywordsOK = 1
endif else begin
  keywordsOK = 0
endelse

if n_params() ne 0 or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'setgroups,group=g,n=n,first=first,last=last,arr=arr,phasewidth=pw,', $
        'phase0=p0,jdwidth=jdw,jd0=jd0[,/silent][,/help]',format='(2a)'
  print,'   Can assign a particular group number by specifying g plus ', $
        'one and only one of the following:',format='(2a)'
  print,'   -  n is the number (counting from 1) of a single stack pair'
  print,'   -  first and last together specify a range of pairs to group'
  print,'      (last will be truncated if it points beyond the end of the stack)'
  print,'   -  arr is an integer vector containing the numbers of the pairs to group'
  print,'   Can instead assign group numbers based on rotation phase (deg) by ', $
        'specifying pw',format='(2a)'
  print,'   -  Optionally p0 gives the starting phase for group 1 ', $
        '(default = 0)',format='(2a)'
  print,'   Can instead assign group numbers based on mean Julian date (days) by ', $
        'specifying jdw',format='(2a)'
  print,'   -  Optionally jd0 gives the starting Julian date for group 1 ', $
        '(default = earliest mean JD in stack)',format='(2a)'
  print,' '
  return
endif

if nstack eq 0 then begin
  print,'ERROR in setgroups: There are no pairs in the stack,', $
        ' so none can be grouped',format='(2a)'
  return
endif

; Create a sorted vector containing the numbers of the stack pairs which will be
; given new group numbers, and another vector containing those group numbers
;
; ngroup      = number of pairs to be regrouped
; groupvec    = ordered list of those pairs' numbers within the stack
; groupnumber = list of those pairs' new group numbers
;
; If grouping the entire single-channel stack by rotation phase or Julian date,
; the group numbers are 1, 2, 3, ..., unless some phases or dates are "missing"
; from the stack, in which case the corresponding group numbers will be skipped.

if n_elements(n) gt 0 then begin
  ngroup = 1L
  groupvec = [long(n)]
  groupnumber = [long(g)]
endif else if n_elements(first) gt 0 then begin
  if first le last then begin
    uselast = long(last) < nstack
    ngroup = uselast - long(first) + 1L
    groupvec = lindgen(ngroup) + long(first)
    groupnumber = replicate(long(g), ngroup)
  endif else begin
    print,'ERROR in setgroups: Must specify range with last >= first'
    return
  endelse
endif else if n_elements(arr) gt 0 then begin
  ngroup = n_elements(arr)
  longarray = long(arr)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroup then begin
    print,"ERROR in setgroups: Don't include duplicate pair numbers in array arr"
    return
  endif
  groupnumber = replicate(long(g), ngroup)
endif else if n_elements(pw) gt 0 then begin
  if n_elements(p0) eq 0 then p0 = 0
  if pw lt 0 or pw ge 360 or p0 lt 0 or p0 ge 360 then begin
    print,'ERROR in setgroups: Both pw and p0 must be in the range [0,360)'
    return
  endif
  ngroup = nstack
  groupvec = lindgen(ngroup) + 1L
  groupnumber = lonarr(ngroup)
  for k=0L,nstack-1 do begin
    phase_diff = (*stack[k]).tags[0].phase - p0
    phase_diff = phase_diff - 360*floor(phase_diff/360.0)  ;  put in range [0,360)
    groupnumber[k] = floor(phase_diff/pw) + 1L
  endfor
endif else begin
  ngroup = nstack
  groupvec = lindgen(nstack) + 1L
  jd = dblarr(ngroup)
  for k=0L,nstack-1 do begin
    jd[k] = getextra('jdmean', stack=(k+1))
    if isnull(jd[k]) then begin
      tags1 = (*stack[k]).tags[0]
      jd[k] = julday(tags1.imm, tags1.idd, tags1.iyy, tags1.rchour, $
                     tags1.rcmin, (tags1.rcsec + tags1.rcnsec/1.0e9))
    endif
  endfor
  if n_elements(jd0) eq 0 then jd0 = min(jd)
  groupnumber = floor( (jd - jd0)/jdw ) - floor( (min(jd) - jd0)/jdw ) + 1L
endelse

; Check that all stack pairs involved actually exist

nmin = min(groupvec)
nmax = max(groupvec)
if nmin lt 1 then begin
  print,'ERROR in setgroups: Pair #',nmin,' is out of the valid stack range (1 - ', $
        nstack,')',format='(a,i0,a,i0,a)'
  return
endif else if nmax gt nstack then begin
  print,'ERROR in setgroups: Pair #',nmax,' is out of the valid stack range (1 - ', $
        nstack,')',format='(a,i0,a,i0,a)'
  return
endif

; Make the changes

for k=0L,ngroup-1 do (*stack[groupvec[k]-1]).group = groupnumber[k]
if not keyword_set(silent) then begin
  print,'Regrouped ',ngroup,' spectra in the pair stack',format='(a,i0,a)'
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setgroups1,group=g,n=n,first=first,last=last,arr=arr,phasewidth=pw,phase0=p0, $
               jdwidth=jdw,jd0=jd0,chan=chan,silent=silent,help=help

; Assign group numbers to some or all spectra in the single-channel stack, either
; by directly specifying the numbers or else by grouping by phase or by time.
; The procedure creates separate groups for the OC vs. SC polarization channel,
; over and above the explicit grouping criterion used.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that a valid combination of keywords has been used

if n_elements(g) gt 0 then begin
  keywordsOK = $
    (n_elements(pw) eq 0 and n_elements(p0) eq 0  $
                         and n_elements(jdw) eq 0 and n_elements(jd0) eq 0) $
    and $
    ( (n_elements(n) gt 0 and $
       (n_elements(first) eq 0 and n_elements(last) eq 0 and n_elements(arr) eq 0)) $
      or $
      (n_elements(first) gt 0 and n_elements(last) gt 0 and $
       n_elements(n) eq 0 and n_elements(arr) eq 0) $
      or $
      (n_elements(arr) gt 0 and $
       (n_elements(n) eq 0 and n_elements(first) eq 0 and n_elements(last) eq 0)) )
endif else if n_elements(n) gt 0 or n_elements(first) gt 0 or $
              n_elements(last) gt 0 or n_elements(arr) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(pw) gt 0 then begin
  keywordsOK = n_elements(jdw) eq 0 and keyword_set(jd0) eq 0
endif else if n_elements(p0) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(jdw) gt 0 then begin
  keywordsOK = 1
endif else begin
  keywordsOK = 0
endelse

if n_params() ne 0 or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'setgroups1,group=g,n=n,first=first,last=last,arr=arr,phasewidth=pw,', $
        'phase0=p0,jdwidth=jdw,jd0=jd0[,chan=1 or 2][,/silent][,/help]',format='(2a)'
  print,'   Can assign a particular group number by specifying g plus ', $
        'one and only one of the following:',format='(2a)'
  print,'   -  n is the number (counting from 1) of a single stack1 spectrum'
  print,'   -  first and last together specify a range of single-channel ', $
        'spectra to group',format='(2a)'
  print,'      (last will be truncated if it points beyond the end of the stack)'
  print,'   -  arr is an integer vector containing the numbers of the spectra ', $
        'to group',format='(2a)'
  print,'   Can instead assign group numbers based on rotation phase (deg) by ', $
        'specifying pw',format='(2a)'
  print,'   -  Optionally p0 gives the starting phase for group 1 ', $
        '(default = 0)',format='(2a)'
  print,'   Can instead assign group numbers based on mean Julian date (days) by ', $
        'specifying jdw',format='(2a)'
  print,'   -  Optionally jd0 gives the starting Julian date for group 1 ', $
        '(default = earliest mean JD in stack)',format='(2a)'
  print,'   Over and above these grouping criteria, separate groups are created ', $
        'for OC vs. SC spectra',format='(2a)'
  print,' '
  return
endif

if nstack1 eq 0 then begin
  print,'ERROR in setgroups1: There are no spectra in the single-channel stack,', $
        ' so none can be grouped',format='(2a)'
  return
endif

; Check which channels to group

if n_elements(chan) eq 0 then begin
  chanstring = ''
endif else if chan eq 1 or chan eq 2 then begin
  chanstring = (chan eq 1) ? 'OC' : 'SC'
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse

; Create a sorted vector containing the numbers of the stack1 spectra which *might*
; be given new group numbers; don't yet take the chan keyword into account.
;
; ngroup      = number of spectra to be regrouped IF all meet the chan criterion
; groupvec    = ordered list of those spectra's numbers within stack1
;
; Revised values and new group numbers will be assigned later after dealing with chan.

if n_elements(n) gt 0 then begin
  ngroup = 1L
  groupvec = [long(n)]
endif else if n_elements(first) gt 0 then begin
  if first gt last then begin
    print,'ERROR in setgroups1: Must specify range with last >= first'
    return
  endif
  uselast = long(last) < nstack1
  ngroup = uselast - long(first) + 1L
  groupvec = lindgen(ngroup) + long(first)
endif else if n_elements(arr) gt 0 then begin
  ngroup = n_elements(arr)
  longarray = long(arr)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroup then begin
    print,"ERROR in setgroups1: Don't include duplicate spectrum numbers ", $
          "in array arr",format='(2a)'
    return
  endif
endif else if n_elements(pw) gt 0 then begin
  if n_elements(p0) eq 0 then p0 = 0
  if pw lt 0 or pw ge 360 or p0 lt 0 or p0 ge 360 then begin
    print,'ERROR in setgroups1: Both pw and p0 must be in the range [0,360)'
    return
  endif
  ngroup = nstack1
  groupvec = lindgen(ngroup) + 1L
endif else begin
  ngroup = nstack1
  groupvec = lindgen(ngroup) + 1L
endelse

; Check that all spectra to be grouped actually exist in stack1

nmin = min(groupvec)
nmax = max(groupvec)
if nmin lt 1 then begin
  print,'ERROR in setgroups1: Spectrum #',nmin, $
        ' is out of the valid stack1 range (1 - ',nstack1,')',format='(a,i0,a,i0,a)'
  return
endif else if nmax gt nstack1 then begin
  print,'ERROR in setgroups1: Spectrum #',nmax, $
        ' is out of the valid stack1 range (1 - ',nstack1,')',format='(a,i0,a,i0,a)'
  return
endif

; Check whether or not all spectra in to be grouped have the same polarization
; and whether or not all spectra in stack1 have the same polarization

firstpol = (*stack1[groupvec[0]-1]).pol
onepol = 1
k = 1L
while onepol and (k lt ngroup) do begin
  onepol = (*stack1[groupvec[k]-1]).pol eq firstpol
  k = k + 1L
endwhile

firstpol = (*stack1[0]).pol
onepolInStack1 = 1
k = 1L
while onepolInStack1 and (k lt nstack1) do begin
  onepolInStack1 = (*stack1[k]).pol eq firstpol
  k = k + 1L
endwhile

; Give up if only one group number has been specified for two polarizations

if (not onepol) and (n_elements(chan) eq 0) $
                and (n_elements(first) gt 0 or n_elements(arr) gt 0) then begin
  print,"ERROR in setgroups1: Can't assign OC and SC spectra to the same group"
  print,'                     Rerun setgroups1 with chan=1 or chan=2 set'
  return
endif

; Create a vector containing the new group numbers, finally taking
; the chan keyword into account
;
; newngroup      = number of spectra to be regrouped
; newgroupvec    = ordered list of those spectra's numbers within stack1
; groupnumber    = list of those spectra's new group numbers
;
; If grouping the entire single-channel stack by rotation phase or Julian date,
; the group numbers are 1, 2, 3, ... if there's only OC or only SC present,
; else (1, 3, 5, ...) for OC and (2, 4, 6, ...) for SC.
;
; (If some phases or dates are "missing" from stack1 then the corresponding
; group numbers will be skipped.)

newgroupvec = lonarr(ngroup)
newngroup = 0L

if n_elements(n) gt 0 or n_elements(first) gt 0 or n_elements(arr) gt 0 then begin
  groupnumber = replicate(long(g), ngroup)
  for k=0L,ngroup-1 do begin
    if (n_elements(chan) eq 0) or ((*stack1[groupvec[k]-1]).pol eq chanstring) then begin
      newgroupvec[newngroup] = groupvec[k]
      newngroup = newngroup + 1L
    endif
  endfor
endif else if n_elements(pw) gt 0 then begin
  groupnumber = lonarr(ngroup)
  for k=0L,nstack1-1 do begin
    if (n_elements(chan) eq 0) or ((*stack1[k]).pol eq chanstring) then begin
      phase_diff = (*stack1[k]).tags[0].phase - p0
      phase_diff = phase_diff - 360*floor(phase_diff/360.0)  ;  put in range [0,360)
      if onepolInStack1 then begin
        groupnumber[newngroup] = floor(phase_diff/pw) + 1L
      endif else begin
        channel = ((*stack1[k]).pol eq 'OC') ? 1L : 2L
        groupnumber[newngroup] = 2*floor(phase_diff/pw) + channel
      endelse
      newgroupvec[newngroup] = k + 1L
      newngroup = newngroup + 1L
    endif
  endfor
endif else begin
  jd = dblarr(ngroup)
  channel = lonarr(ngroup)
  for k=0L,nstack1-1 do begin
    if (n_elements(chan) eq 0) or ((*stack1[k]).pol eq chanstring) then begin
      jd[newngroup] = getextra1('jdmean', stack=(k+1))
      if isnull(jd[newngroup]) then begin
        tags1 = (*stack1[k]).tags
        jd[newngroup] = julday(tags1.imm, tags1.idd, tags1.iyy, tags1.rchour, $
                               tags1.rcmin, (tags1.rcsec + tags1.rcnsec/1.0e9))
      endif
      channel[newngroup] = ((*stack1[k]).pol eq 'OC') ? 1L : 2L
      newgroupvec[newngroup] = k + 1L
      newngroup = newngroup + 1L
    endif
  endfor
  jd = (newngroup gt 0) ? jd[0:newngroup-1] : [0.0D]
  if n_elements(jd0) eq 0 then jd0 = min(jd)
  if onepolInStack1 then $
    groupnumber = floor( (jd - jd0)/jdw ) - floor( (min(jd) - jd0)/jdw ) + 1L  $
  else $
    groupnumber = 2*(floor( (jd - jd0)/jdw ) - floor( (min(jd) - jd0)/jdw )) + channel
endelse

; Make sure that applying the chan keyword didn't eliminate all spectra

if newngroup gt 0 then begin
  newgroupvec = newgroupvec[0:newngroup-1]
  groupnumber = groupnumber[0:newngroup-1]
endif else begin
  print,'setgroups1: No spectra in the single-channel stack meet your criteria'
  return
endelse

; Make the changes

for k=0L,newngroup-1 do (*stack1[newgroupvec[k]-1]).group = groupnumber[k]
if not keyword_set(silent) then begin
 print,'Regrouped ',newngroup,' spectra in the single-channel stack', $
       format='(a,i0,a)'
endif

end
               
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setgroupsi,group=g,n=n,first=first,last=last,arr=arr,phasewidth=pw,phase0=p0, $
               jdwidth=jdw,jd0=jd0,chan=chan,silent=silent,help=help

; Assign group numbers to some or all images in the image stack, either
; by directly specifying the numbers or else by grouping by phase or by time.
; The procedure creates separate groups for the OC vs. SC polarization channel,
; over and above the explicit grouping criterion used.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that a valid combination of keywords has been used

if n_elements(g) gt 0 then begin
  keywordsOK = $
    (n_elements(pw) eq 0 and n_elements(p0) eq 0  $
                         and n_elements(jdw) eq 0 and n_elements(jd0) eq 0) $
    and $
    ( (n_elements(n) gt 0 and $
       (n_elements(first) eq 0 and n_elements(last) eq 0 and n_elements(arr) eq 0)) $
      or $
      (n_elements(first) gt 0 and n_elements(last) gt 0 and $
       n_elements(n) eq 0 and n_elements(arr) eq 0) $
      or $
      (n_elements(arr) gt 0 and $
       (n_elements(n) eq 0 and n_elements(first) eq 0 and n_elements(last) eq 0)) )
endif else if n_elements(n) gt 0 or n_elements(first) gt 0 or $
              n_elements(last) gt 0 or n_elements(arr) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(pw) gt 0 then begin
  keywordsOK = n_elements(jdw) eq 0 and keyword_set(jd0) eq 0
endif else if n_elements(p0) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(jdw) gt 0 then begin
  keywordsOK = 1
endif else begin
  keywordsOK = 0
endelse

if n_params() ne 0 or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'setgroupsi,group=g,n=n,first=first,last=last,arr=arr,phasewidth=pw,', $
        'phase0=p0,jdwidth=jdw,jd0=jd0[,chan=1 or 2][,/silent][,/help]',format='(2a)'
  print,'   Can assign a particular group number by specifying g plus ', $
        'one and only one of the following:',format='(2a)'
  print,'   -  n is the number (counting from 1) of a single stacki image'
  print,'   -  first and last together specify a range of images to group'
  print,'      (last will be truncated if it points beyond the end of the stack)'
  print,'   -  arr is an integer vector containing the numbers of the images ', $
        'to group',format='(2a)'
  print,'   Can instead assign group numbers based on rotation phase (deg) by ', $
        'specifying pw',format='(2a)'
  print,'   -  Optionally p0 gives the starting phase for group 1 ', $
        '(default = 0)',format='(2a)'
  print,'   Can instead assign group numbers based on mean Julian date (days) by ', $
        'specifying jdw',format='(2a)'
  print,'   -  Optionally jd0 gives the starting Julian date for group 1 ', $
        '(default = earliest mean JD in stack)',format='(2a)'
  print,'   Over and above these grouping criteria, separate groups are created ', $
        'for OC vs. SC images',format='(2a)'
  print,' '
  return
endif

if nstacki eq 0 then begin
  print,'ERROR in setgroupsi: There are no images in the image stack,', $
        ' so none can be grouped',format='(2a)'
  return
endif

; Check which channels to group

if n_elements(chan) eq 0 then begin
  chanstring = ''
endif else if chan eq 1 or chan eq 2 then begin
  chanstring = (chan eq 1) ? 'OC' : 'SC'
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse

; Create a sorted vector containing the numbers of the stacki images which *might*
; be given new group numbers; don't yet take the chan keyword into account.
;
; ngroup      = number of images to be regrouped IF all meet the chan criterion
; groupvec    = ordered list of those images's numbers within stacki
;
; Revised values and new group numbers will be assigned later after dealing with chan.

if n_elements(n) gt 0 then begin
  ngroup = 1L
  groupvec = [long(n)]
endif else if n_elements(first) gt 0 then begin
  if first gt last then begin
    print,'ERROR in setgroupsi: Must specify range with last >= first'
    return
  endif
  uselast = long(last) < nstacki
  ngroup = uselast - long(first) + 1L
  groupvec = lindgen(ngroup) + long(first)
endif else if n_elements(arr) gt 0 then begin
  ngroup = n_elements(arr)
  longarray = long(arr)
  groupvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupvec) lt ngroup then begin
    print,"ERROR in setgroupsi: Don't include duplicate image numbers ", $
          "in array arr",format='(2a)'
    return
  endif
endif else if n_elements(pw) gt 0 then begin
  if n_elements(p0) eq 0 then p0 = 0
  if pw lt 0 or pw ge 360 or p0 lt 0 or p0 ge 360 then begin
    print,'ERROR in setgroupsi: Both pw and p0 must be in the range [0,360)'
    return
  endif
  ngroup = nstacki
  groupvec = lindgen(ngroup) + 1L
endif else begin
  ngroup = nstacki
  groupvec = lindgen(ngroup) + 1L
endelse

; Check that all images to be grouped actually exist in stacki

nmin = min(groupvec)
nmax = max(groupvec)
if nmin lt 1 then begin
  print,'ERROR in setgroupsi: Image #',nmin, $
        ' is out of the valid stacki range (1 - ',nstacki,')',format='(a,i0,a,i0,a)'
  return
endif else if nmax gt nstacki then begin
  print,'ERROR in setgroupsi: Image #',nmax, $
        ' is out of the valid stacki range (1 - ',nstacki,')',format='(a,i0,a,i0,a)'
  return
endif

; Check whether or not all images in to be grouped have the same polarization
; and whether or not all images in stacki have the same polarization

firstpol = (*stacki[groupvec[0]-1]).pol
onepol = 1
k = 1L
while onepol and (k lt ngroup) do begin
  onepol = (*stacki[groupvec[k]-1]).pol eq firstpol
  k = k + 1L
endwhile

firstpol = (*stacki[0]).pol
onepolInStacki = 1
k = 1L
while onepolInStacki and (k lt nstacki) do begin
  onepolInStacki = (*stacki[k]).pol eq firstpol
  k = k + 1L
endwhile

; Give up if only one group number has been specified for two polarizations

if (not onepol) and (n_elements(chan) eq 0) $
                and (n_elements(first) gt 0 or n_elements(arr) gt 0) then begin
  print,"ERROR in setgroupsi: Can't assign OC and SC images to the same group"
  print,'                     Rerun setgroupsi with chan=1 or chan=2 set'
  return
endif

; Create a vector containing the new group numbers, finally taking
; the chan keyword into account
;
; newngroup      = number of images to be regrouped
; newgroupvec    = ordered list of those images's numbers within stacki
; groupnumber    = list of those images's new group numbers
;
; If grouping the entire image stack by rotation phase or Julian date,
; the group numbers are 1, 2, 3, ... if there's only OC or only SC present,
; else (1, 3, 5, ...) for OC and (2, 4, 6, ...) for SC.
;
; (If some phases or dates are "missing" from stacki then the corresponding
; group numbers will be skipped.)

newgroupvec = lonarr(ngroup)
newngroup = 0L

if n_elements(n) gt 0 or n_elements(first) gt 0 or n_elements(arr) gt 0 then begin
  groupnumber = replicate(long(g), ngroup)
  for k=0L,ngroup-1 do begin
    if (n_elements(chan) eq 0) or ((*stacki[groupvec[k]-1]).pol eq chanstring) then begin
      newgroupvec[newngroup] = groupvec[k]
      newngroup = newngroup + 1L
    endif
  endfor
endif else if n_elements(pw) gt 0 then begin
  groupnumber = lonarr(ngroup)
  for k=0L,nstacki-1 do begin
    if (n_elements(chan) eq 0) or ((*stacki[k]).pol eq chanstring) then begin
      phase_diff = (*stacki[k]).tags[0].phase - p0
      phase_diff = phase_diff - 360*floor(phase_diff/360.0)  ;  put in range [0,360)
      if onepolInStacki then begin
        groupnumber[newngroup] = floor(phase_diff/pw) + 1L
      endif else begin
        channel = ((*stacki[k]).pol eq 'OC') ? 1L : 2L
        groupnumber[newngroup] = 2*floor(phase_diff/pw) + channel
      endelse
      newgroupvec[newngroup] = k + 1L
      newngroup = newngroup + 1L
    endif
  endfor
endif else begin
  jd = dblarr(ngroup)
  channel = lonarr(ngroup)
  for k=0L,nstacki-1 do begin
    if (n_elements(chan) eq 0) or ((*stacki[k]).pol eq chanstring) then begin
      jd[newngroup] = getextrai('jdmean', stack=(k+1))
      if isnull(jd[newngroup]) then begin
        print,'setgroupsi: stacki image #',k+1," lacks the 'jdmean' extra tag", $
              format='(a,i0,a)'
        return
      endif
      channel[newngroup] = ((*stacki[k]).pol eq 'OC') ? 1L : 2L
      newgroupvec[newngroup] = k + 1L
      newngroup = newngroup + 1L
    endif
  endfor
  jd = (newngroup gt 0) ? jd[0:newngroup-1] : [0.0D]
  if n_elements(jd0) eq 0 then jd0 = min(jd)
  if onepolInStacki then $
    groupnumber = floor( (jd - jd0)/jdw ) - floor( (min(jd) - jd0)/jdw ) + 1L  $
  else $
    groupnumber = 2*(floor( (jd - jd0)/jdw ) - floor( (min(jd) - jd0)/jdw )) + channel
endelse

; Make sure that applying the chan keyword didn't eliminate all images

if newngroup gt 0 then begin
  newgroupvec = newgroupvec[0:newngroup-1]
  groupnumber = groupnumber[0:newngroup-1]
endif else begin
  print,'setgroupsi: No images in the image stack meet your criteria'
  return
endelse

; Make the changes

for k=0L,newngroup-1 do (*stacki[newgroupvec[k]-1]).group = groupnumber[k]
if not keyword_set(silent) then begin
  print,'Regrouped ',newngroup,' images in the image stack', $
        format='(a,i0,a)'
endif

end
               
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sumspec,group=g,ming=ming,maxg=maxg,arrayg=ag,all=all, $
            noweight=noweight,chan=chan,help=help

; Create one or more weighted spectral sums of stack pairs in one or more groups,
; then add the sum(s) to the stack with new group number(s) attached.
;
; Calling sumspec with no arguments creates a separate sum for each group which 
; contains two or more pairs.  Calling sumspec with keywords g, ming and maxg,
; or ag set creates a single sum of all pairs in the specified group(s).  Calling 
; "sumspec,/all" creates a single sum of all pairs in the stack.
;
; /noweight creates unweighted sums
;
; If keyword chan is set then the the single-channel stack gets the resulting
; single-channel sum(s).
;
; 05 Jan 2002: Fixed bug in computing phase2_min, phase2_max, az2_min, az2_max
; 13 Jul 2002: Fixed bug in computing new sdev for sum of unhopped spectra
; 29 Nov 2002: Fixed bug in computing rmsm
; 17 Nov 2005: Fixed bug in assigning doppl
; 02 May 2005: Added noweight keyword
; 25 Jun 2006: Frequency refers to center of bin, not left edge
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

maxextratags = 1000
too_low_sdev = 1.0d-10

monthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

; Check if a valid set of keywords is set

if n_elements(g) gt 0 then begin
  keywordsOK = n_elements(ming) eq 0 and n_elements(maxg) eq 0 and n_elements(ag) eq 0  $
                                                               and not keyword_set(all)
endif else if n_elements(ming) gt 0 then begin
  keywordsOK = n_elements(maxg) gt 0 and n_elements(ag) eq 0 and not keyword_set(all)
endif else if n_elements(maxg) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(ag) gt 0 then begin
  keywordsOK = not keyword_set(all)
endif else begin
  keywordsOK = 1
endelse

if n_params() ne 0 or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'sumspec,group=g,ming=ming,maxg=maxg,arrayg=ag,all=all[,/noweight]', $
        '[,chan=1 or 2][,/help]',format='(2a)'
  print,' '
  print,'    Calling sumspec,/all creates a single sum of all pairs in the stack'
  print,'    Setting keywords g, ming and maxg, or ag yields a single sum of ', $
        'all pairs in the specified group(s)',format='(2a)'
  print,'    Keywords all, g, ming and maxg, and ag are four mutually ', $
        'exclusive choices',format='(2a)'
  print,'    Calling sumspec with none of those keywords set creates ', $
        'a sum for each group in the stack',format='(2a)'
  print,' '
  print,'    /noweight produces unweighted sums'
  print,' '
  return
endif

if nstack eq 0 then begin
  print,'ERROR in sumspec: The pair stack is empty'
  return
endif

; Create an array which contains, for each output sum,
; a sorted vector of the group numbers which will contribute to that sum
;
; nsums       = number of sums to do
; ngroups     = number of groups to include in each sum (scalar, not vector)
; groupnumber = array of those group numbers:
;                    groupnumber[i,*] = ordered list of group numbers which
;                                       contribute to the ith sum

allgroups = stackgroups()
n_all = n_elements(allgroups)

if n_elements(g) gt 0 then begin
  nsums = 1L
  ngroups = 1L
  groupnumber = lonarr(nsums,ngroups)
  groupnumber[0,0] = long(g)
endif else if n_elements(ming) gt 0 then begin
  if ming gt maxg then begin
    print,'ERROR in sumspec: Must specify range with ming <= maxg'
    return
  endif
  nsums = 1L
  groupsvec = lonarr(n_all)
  ngroups = 0L
  for k=0L,n_all-1 do begin
    if allgroups[k] ge ming and allgroups[k] le maxg then begin
      groupsvec[ngroups] = allgroups[k]
      ngroups = ngroups + 1L
    endif
  endfor
  if ngroups eq 0 then begin
    print,'sumspec: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
  groupnumber = lonarr(nsums,ngroups)
  groupnumber[0,*] = groupsvec[0:ngroups-1]
endif else if n_elements(ag) gt 0 then begin
  nsums = 1L
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupsvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupsvec) lt ngroups then begin
    print,"ERROR in sumspec: Don't include duplicate group numbers ", $
          "in array ag",format='(2a)'
    return
  endif
  groupnumber = lonarr(nsums,ngroups)
  groupnumber[0,*] = groupsvec
endif else if keyword_set(all) then begin
  nsums = 1L
  ngroups = n_all
  groupnumber = lonarr(nsums,ngroups)
  groupnumber[0,*] = allgroups
endif else begin
  nsums = n_all
  ngroups = 1L
  groupnumber = lonarr(nsums,ngroups)
  groupnumber[*,0] = allgroups
endelse

; Check that all specified groups actually exist

for i=0L,nsums-1 do begin
  for j=0L,ngroups-1 do begin
    groupToCheck = groupnumber[i,j]
    groupnum = where(allgroups eq groupToCheck, count)
    if count eq 0 then begin
      print,'ERROR in sumspec: Group #',groupToCheck, $
            ' is not present in the pair stack',format='(a,i0,a)'
      return
    endif
  endfor
endfor

; Initialize some arrays which will be used for constructing the sums

n_in = lonarr(nsums)
newndata = lonarr(nsums)
newnpol = lonarr(nsums)
newtags = replicate(blanktags1(), maxchan, nsums)
newtags.jsnr1 = 999999999L
newtags.jsnr2 = -999999999L
newntags = lonarr(nsums)
extratag = {format:'', name:'', value:'', comment:''}
newextratags = replicate(extratag, nsums, maxextratags)
newnextra = lonarr(nsums)
newtname = strarr(nsums)
timezone = strarr(nsums)
timezoneOK = intarr(nsums) + 1L
sumweights = dblarr(maxchan,nsums)
sumvariance = dblarr(maxchan,nsums)
rmsc_sum = dblarr(maxchan,nsums)
phase1_mean = fltarr(maxchan,nsums)
phase2_mean = fltarr(maxchan,nsums)
phase1_min = fltarr(maxchan,nsums) + 999.9
phase1_max = fltarr(maxchan,nsums) - 999.9
phase2_min = fltarr(maxchan,nsums) + 999.9
phase2_max = fltarr(maxchan,nsums) - 999.9
az1_mean = fltarr(maxchan,nsums)
az2_mean = fltarr(maxchan,nsums)
az1_min = fltarr(maxchan,nsums) + 999.9
az1_max = fltarr(maxchan,nsums) - 999.9
az2_min = fltarr(maxchan,nsums) + 999.9
az2_max = fltarr(maxchan,nsums) - 999.9
jdstart = dblarr(nsums) + 999999999L
jdmean = dblarr(nsums)
calmean = strarr(nsums)
jdend = dblarr(nsums) - 999999999L
jdmeanOK = intarr(nsums) + 1L
dist_min = fltarr(nsums) + 999999.9
dist_max = fltarr(nsums) - 999999.9
dist_mean = fltarr(nsums)
dist_meanOK = intarr(nsums) + 1L
ra1_min = fltarr(nsums) + 999.9
ra1_max = fltarr(nsums) - 999.9
ra2_min = fltarr(nsums) + 999.9
ra2_max = fltarr(nsums) - 999.9
ra1_mean = fltarr(nsums)
ra2_mean = fltarr(nsums)
ra_meanOK = intarr(nsums) + 1L
dec_min = fltarr(nsums) + 999.9
dec_max = fltarr(nsums) - 999.9
dec_mean = fltarr(nsums)
dec_meanOK = intarr(nsums) + 1L

; Initialize the weighted sums; use pointers in case different sums have
; different numbers of spectral values

wsum = ptrarr(nsums, /allocate_heap)

; Go through the stack, pair by pair, and construct the sum(s)

for n=0L,nstack-1 do begin

  ; See if this pair is included in one of the groups which will be summed

  gnum1D = where(groupnumber eq (*stack[n]).group, count)

  if count gt 0 then begin

    ; This pair contributes to a sum: figure out WHICH sum

    gnum1D = gnum1D[0]
    nsum = gnum1D mod nsums
    ngroup = gnum1D/nsums

    ; Get the various elements of this pair

    pair = (*stack[n]).spec
    ndata = (*stack[n]).ndata
    tags = (*stack[n]).tags
    ntags = (*stack[n]).ntags
    extratags = (*stack[n]).extratags
    nextra = (*stack[n]).nextra
    tname = (*stack[n]).tname
    npol = n_elements(tags)
    np = npol-1 ; for indexing

    ; If this is the first pair contributing to this sum, initialize some values

    if n_in[nsum] eq 0 then begin
; Check which channels to sum

      sumpol = intarr(npol)
      if n_elements(chan) eq 0 then begin
        sumpol = sumpol + 1
      endif else if chan ge 1 and chan le npol then begin
        sumpol[chan-1] = 1
        polstring = chanstrings[chan-1]
      endif else begin
        print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
        return
      endelse

; We can use the sumpol vector as a mask to ensure that we only flag parameter
; differences if they occur in the polarization channel(s) of interest

      pmask = transpose(sumpol)


      ; Create an array which will hold the weighted spectra for this sum
      ;
      ; Since different sums may involve spectra of different lengths, and
      ; since we don't want to presuppose a maximum spectrum length, use a
      ; pointer to reference the array.  But since there doesn't seem to be
      ; any way to do that directly, fool IDL by making a structure whose
      ; only field is the spectral array, then point to that structure.
      ; (Can't yet create the complete output structure: We don't know
      ; how long the extratags array will be.)

      dummyStruc = {spec:fltarr(npol,ndata)}
      *wsum[nsum] = dummyStruc

      ; Use the first pair contributing to this sum to get values for various
      ; pair parameters, so that we can check whether or not later pairs which
      ; contribute to the same sum have the same values

      newndata[nsum] = ndata
      newntags[nsum] = ntags
      newnpol[nsum]  = npol
      newtags[0:np,nsum].xjcen = tags.xjcen
      newtags[0:np,nsum].posfr = tags.posfr
      newtags[0:np,nsum].dfreq = tags.dfreq
      newtags[0:np,nsum].doppl = tags.doppl
      newtags[0:np,nsum].zepch = tags.zepch
      newtags[0:np,nsum].obs = tags.obs
      newtags[0:np,nsum].itar = tags.itar
      newtags[0:np,nsum].lfft = tags.lfft
      newtags[0:np,nsum].igw = tags.igw
      newtags[0:np,nsum].kpts = tags.kpts
      newtags[0:np,nsum].nfreq = tags.nfreq
      newtags[0:np,nsum].frstep = tags.frstep
      newtags[0:np,nsum].color = tags.color
      newtags[0:np,nsum].freq1 = tags.freq1
      for k=0L,nextra-1 do begin
        for j=0L,3 do newextratags[nsum,k].(j) = extratags[k].(j)
      endfor
      newnextra[nsum] = nextra
      newtname[nsum] = tname
      timezone[nsum] = getextra('timezone',stack=(n+1))

    endif

    ; Since this is a weighted sum, and weight = 1/sdev^2,
    ; make sure that sdev isn't zero (or really tiny)

    if not keyword_set(noweight) then begin
      if total(pmask # (abs(1.0D*tags.sdev) lt too_low_sdev)) gt 0 then begin
        print,'ERROR in sumspec: Stack pair #',n+1,' has channel(s) with sdev < ', $
              too_low_sdev,' km^2',format='(a,i0,a,e7.1,a)'
        return
      endif
    endif

    ; Insist that essential parameters have the same values
    ; for all pairs contributing to this sum

    if tname ne newtname[nsum] then begin
      print,'ERROR in sumspec: Pairs with different target names contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if ndata ne newndata[nsum] then begin
      print,'ERROR in sumspec: Pairs with different numbers of spectral points ', $
            'contribute to sum #',nsum+1,format='(2a,i0)'
      return
    endif else if ntags ne newntags[nsum] then begin
      print,'ERROR in sumspec: Pairs with different numbers of tags contribute to sum #', $
            nsum+1,format='(a,i0)'
      return
    endif else if npol ne newnpol[nsum] then begin
      print,'ERROR in sumspec: Pairs with different numbers of pols contribute to sum #', $
            nsum+1,format='(a,i0)'
      return
    endif else if total(pmask # (tags.xjcen ne newtags[0:np,nsum].xjcen)) gt 0 then begin
      print,'ERROR in sumspec: Pairs with different xjcen tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if total(pmask # (tags.posfr ne newtags[0:np,nsum].posfr)) gt 0 then begin
      print,'ERROR in sumspec: Pairs with different posfr tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if total(pmask # (tags.dfreq ne newtags[0:np,nsum].dfreq)) gt 0 then begin
      print,'ERROR in sumspec: Pairs with different dfreq tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif

    ; Construct the weighted mean spectra
    ;
    ; Note that multiplying the unnormalized spectra (sdev * pair)
    ; by the weight (1/sdev^2) is equivalent to dividing pair/sdev

    if keyword_set(noweight) then begin
      (*wsum[nsum]).spec = (*wsum[nsum]).spec $
                           + float( pair*(tags.sdev # replicate(1.0D, ndata)) )
    endif else begin
      (*wsum[nsum]).spec = (*wsum[nsum]).spec $
                           + float( pair/(tags.sdev # replicate(1.0D, ndata)) )
    endelse

    ; Accumulate some parameters

    weight = (keyword_set(noweight)) ? dblarr(npol) + 1.d0 : 1/(1.0D*tags.sdev)^2
    sumweights[0:np,nsum] = sumweights[0:np,nsum] + weight
    rmsc_sum[0:np,nsum] = rmsc_sum[0:np,nsum] + 1/(1.0D*tags.rmsc)^2
    newtags[0:np,nsum].nffts = newtags[0:np,nsum].nffts + tags.nffts
    newtags[0:np,nsum].tau = newtags[0:np,nsum].tau + tags.tau

    ; Work with two versions of the phase and the azimuth,
    ; using ranges [0,360) vs. [180,540)

    phase1 = tags.phase - 360*floor(tags.phase/360.0)
    phase2 = phase1 - 360*floor( (phase1 - 180)/360.0 )
    az1 = tags.azim - 360*floor(tags.azim/360.0)
    az2 = az1 - 360*floor( (az1 - 180)/360.0 )

    ; Look for minimum or maximum values for some parameters

    newtags[0:np,nsum].jsnr1 = newtags[0:np,nsum].jsnr1 < tags.jsnr1
    newtags[0:np,nsum].jsnr2 = newtags[0:np,nsum].jsnr2 > tags.jsnr2
    jd_in = getextra('jdstart', stack=(n+1))
    jdstart[nsum] = (notnull(jd_in)) ? (jdstart[nsum] < jd_in) : -999999999L
    jd_in = getextra('jdend', stack=(n+1))
    jdend[nsum] = (notnull(jd_in)) ? (jdend[nsum] > jd_in) : 999999999L
    phase1_min[0:np,nsum] = phase1_min[0:np,nsum] < phase1
    phase1_max[0:np,nsum] = phase1_max[0:np,nsum] > phase1
    phase2_min[0:np,nsum] = phase2_min[0:np,nsum] < phase2
    phase2_max[0:np,nsum] = phase2_max[0:np,nsum] > phase2
    az1_min[0:np,nsum] = az1_min[0:np,nsum] < az1
    az1_max[0:np,nsum] = az1_max[0:np,nsum] > az1
    az2_min[0:np,nsum] = az2_min[0:np,nsum] < az2
    az2_max[0:np,nsum] = az2_max[0:np,nsum] > az2
    dist_in = getextra('distmin', stack=(n+1))
    dist_min[nsum] = (notnull(dist_in)) ? (dist_min[nsum] < dist_in) : -999999.9
    dist_in = getextra('distmax', stack=(n+1))
    dist_max[nsum] = (notnull(dist_in)) ? (dist_max[nsum] > dist_in) : 999999.9
    ra1_in = getextra('ramin', stack=(n+1))
    if notnull(ra1_in) then begin
      ra1_min[nsum] = ra1_min[nsum] < ra1_in
      ra1_max[nsum] = ra1_max[nsum] > ra1_in
      ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
      ra2_min[nsum] = ra2_min[nsum] < ra2_in
      ra2_max[nsum] = ra2_max[nsum] > ra2_in
    endif else begin
      ra1_min[nsum] = -999.9
      ra1_max[nsum] = 999.9
      ra2_min[nsum] = -999.9
      ra2_max[nsum] = 999.9
    endelse
    ra1_in = getextra('ramax', stack=(n+1))
    if notnull(ra1_in) then begin
      ra1_min[nsum] = ra1_min[nsum] < ra1_in
      ra1_max[nsum] = ra1_max[nsum] > ra1_in
      ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
      ra2_min[nsum] = ra2_min[nsum] < ra2_in
      ra2_max[nsum] = ra2_max[nsum] > ra2_in
    endif else begin
      ra1_min[nsum] = -999.9
      ra1_max[nsum] = 999.9
      ra2_min[nsum] = -999.9
      ra2_max[nsum] = 999.9
    endelse
    dec_in = getextra('decmin', stack=(n+1))
    dec_min[nsum] = (notnull(dec_in)) ? (dec_min[nsum] < dec_in) : -999.9
    dec_in = getextra('decmax', stack=(n+1))
    dec_max[nsum] = (notnull(dec_in)) ? (dec_max[nsum] > dec_in) : 999.9

    ; Construct weighted means of some parameters

    newtags[0:np,nsum].elev = newtags[0:np,nsum].elev + weight*tags.elev
    newtags[0:np,nsum].rttim = newtags[0:np,nsum].rttim + weight*tags.rttim
    newtags[0:np,nsum].trpwr = newtags[0:np,nsum].trpwr + weight*tags.trpwr
    newtags[0:np,nsum].tsys = newtags[0:np,nsum].tsys + weight*tags.tsys
    newtags[0:np,nsum].gain = newtags[0:np,nsum].gain + weight*tags.gain
    sumvariance[0:np,nsum] = sumvariance[0:np,nsum] + (1.0D*tags.sdev)^2
    phase1_mean[0:np,nsum] = phase1_mean[0:np,nsum] + weight*phase1
    phase2_mean[0:np,nsum] = phase2_mean[0:np,nsum] + weight*phase2
    az1_mean[0:np,nsum] = az1_mean[0:np,nsum] + weight*az1
    az2_mean[0:np,nsum] = az2_mean[0:np,nsum] + weight*az2
    if timezoneOK[nsum] then begin
      tempzone = getextra('timezone',stack=(n+1))
      timezoneOK[nsum] = (tempzone eq timezone[nsum])
    endif
    if jdmeanOK[nsum] then begin
      jd_in = getextra('jdmean', stack=(n+1))
      if notnull(jd_in) then begin
        jdmean[nsum] = jdmean[nsum] + weight[0]*jd_in
      endif else begin
        jdmeanOK[nsum] = 0
      endelse
    endif
    if dist_meanOK[nsum] then begin
      dist_in = getextra('distmean', stack=(n+1))
      if notnull(dist_in) then begin
        dist_mean[nsum] = dist_mean[nsum] + weight[0]*dist_in
      endif else begin
        dist_meanOK[nsum] = 0
      endelse
    endif
    if ra_meanOK[nsum] then begin
      ra1_in = getextra('ramean', stack=(n+1))
      if notnull(ra1_in) then begin
        ra1_mean[nsum] = ra1_mean[nsum] + weight[0]*ra1_in
        ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
        ra2_mean[nsum] = ra2_mean[nsum] + weight[0]*ra2_in
      endif else begin
        ra_meanOK[nsum] = 0
      endelse
    endif
    if dec_meanOK[nsum] then begin
      dec_in = getextra('decmean', stack=(n+1))
      if notnull(dec_in) then begin
        dec_mean[nsum] = dec_mean[nsum] + weight[0]*dec_in
      endif else begin
        dec_meanOK[nsum] = 0
      endelse
    endif

    ; Check whether nonessential tags, such as doppl (tx offset), have the same values
    ; for all pairs contributing to this sum.  If not, reset that tag to a dummy value.

    if total(pmask # (tags.doppl ne newtags[0:np,nsum].doppl)) gt 0 then $
           newtags[0:np,nsum].doppl = 0.0
    if total(pmask # (tags.zepch ne newtags[0:np,nsum].zepch)) gt 0 then $
           newtags[0:np,nsum].zepch = 0.0
    if total(pmask # (tags.obs ne newtags[0:np,nsum].obs)) gt 0 then $
           newtags[0:np,nsum].obs = 0L
    if total(pmask # (tags.itar ne newtags[0:np,nsum].itar)) gt 0 then $
           newtags[0:np,nsum].itar = 0L
    if total(pmask # (tags.lfft ne newtags[0:np,nsum].lfft)) gt 0 then $
           newtags[0:np,nsum].lfft = 0L
    if total(pmask # (tags.igw ne newtags[0:np,nsum].igw)) gt 0 then $
           newtags[0:np,nsum].igw = 0L
    if total(pmask # (tags.kpts ne newtags[0:np,nsum].kpts)) gt 0 then $
           newtags[0:np,nsum].kpts = 0L
    if total(pmask # (tags.nfreq ne newtags[0:np,nsum].nfreq)) gt 0 then $
           newtags[0:np,nsum].nfreq = 0L
    if total(pmask # (tags.frstep ne newtags[0:np,nsum].frstep)) gt 0 then $
           newtags[0:np,nsum].frstep = 0.0
    if total(pmask # (tags.color ne newtags[0:np,nsum].color)) gt 0 then $
           newtags[0:np,nsum].color = -1L
    if total(pmask # (tags.freq1 ne newtags[0:np,nsum].freq1)) gt 0 then $
           newtags[0:np,nsum].freq1 = 0.0

    ; Go through the array of extra tags and find the ones which are present and have
    ; the same format and value for all pairs contributing to this sum: those are the
    ; ones we'll keep.  Similarly, keep a comment line only if every pair has the same
    ; line.  Finally, if formats and values match for a given extra tag but trailing
    ; comments don't match, keep the extra tag but erase the trailing comment.

    for k=newnextra[nsum]-1L,0,-1L do begin

      newextraformat = strlowcase(newextratags[nsum,k].format)
      newextraname = strlowcase(newextratags[nsum,k].name)
      newextravalue = strlowcase(newextratags[nsum,k].value)
      newextracomment = strlowcase(newextratags[nsum,k].comment)
      keepextra = 0
      keepcomment = 0
      extranum = where(extratags.name eq newextraname, extracount)
      if extracount gt 0 then begin
        for j=0L,extracount-1 do begin
          extraformat = strlowcase(extratags[extranum[j]].format)
          extraname = newextraname
          extravalue = strlowcase(extratags[extranum[j]].value)
          extracomment = strlowcase(extratags[extranum[j]].comment)
          match = ((notnull(extraname)) and $
                   (extraformat eq newextraformat) and (extravalue eq newextravalue)) or $
                  ((isnull(extraname)) and (extracomment eq newextracomment))
          keepextra = keepextra or match
          keepcomment = keepcomment or (match and (extracomment eq newextracomment))
        endfor
      endif

      if not keepextra then begin

        ; Don't create an error by trying to remove the last array element

        if newnextra[nsum] gt 1 then begin
          if k eq 0 then begin
            newextratags[nsum,0:newnextra[nsum]-2] = newextratags[nsum,1:newnextra[nsum]-1]
          endif else if k lt newnextra[nsum]-1 then begin
            newextratags[nsum,k:newnextra[nsum]-2] = newextratags[nsum,k+1:newnextra[nsum]-1]
          endif
        endif

        newnextra[nsum] = newnextra[nsum] - 1

      endif else if not keepcomment then begin
        newextratags[nsum,k].comment = ''
      endif

    endfor

    ; Increment the number of pairs contributing to this sum

    n_in[nsum] = n_in[nsum] + 1L

  endif  ;  count gt 0?

; Look at the next pair in the stack

endfor

; All stack pairs have been inspected

for nsum=0L,nsums-1 do begin

  if n_in[nsum] ge 2 then begin

    ; Normalize the weighted mean spectra and get the rms noise (sdev)
    ;
    ; Note that dividing by the sum of all weights, then normalizing
    ; by dividing by sdev = 1/sqrt(sum of all weights), is equivalent
    ; to multiplying by sdev

    np = newnpol[nsum]-1
    if keyword_set(noweight) then begin
      sdev = sqrt(sumvariance[0:np,nsum])
      (*wsum[nsum]).spec = float( (*wsum[nsum]).spec / (sdev # replicate(1.0D,newndata[nsum])) )
      newtags[0:np,nsum].sdev = sdev/n_in[nsum]
    endif else begin
      sdev = 1.0/sqrt(sumweights[0:np,nsum])
      (*wsum[nsum]).spec = (sdev # replicate(1.0,newndata[nsum])) * (*wsum[nsum]).spec
      newtags[0:np,nsum].sdev = sdev
    endelse

    ; Complete the weighted means for various tags

    newtags[0:np,nsum].elev = newtags[0:np,nsum].elev / sumweights[0:np,nsum]
    newtags[0:np,nsum].rttim = newtags[0:np,nsum].rttim / sumweights[0:np,nsum]
    newtags[0:np,nsum].trpwr = newtags[0:np,nsum].trpwr / sumweights[0:np,nsum]
    newtags[0:np,nsum].tsys = newtags[0:np,nsum].tsys / sumweights[0:np,nsum]
    newtags[0:np,nsum].gain = newtags[0:np,nsum].gain / sumweights[0:np,nsum]
    if jdmeanOK[nsum] then jdmean[nsum] = jdmean[nsum] / sumweights[0,nsum]
    if dist_meanOK[nsum] then dist_mean[nsum] = dist_mean[nsum] / sumweights[0,nsum]
    if ra_meanOK[nsum] then begin
      ra1_mean[nsum] = ra1_mean[nsum] / sumweights[0,nsum]
      ra2_mean[nsum] = ra2_mean[nsum] / sumweights[0,nsum]
    endif
    if dec_meanOK[nsum] then dec_mean[nsum] = dec_mean[nsum] / sumweights[0,nsum]

    ; Deal with rotation phase and azimuth:
    ;
    ; If the pairs contributing to this sum have phases which are more
    ; "concentrated" when considered over the range [180,540) rather than [0,360),
    ; then compute the mean phase using the former range (and then take mod 360).
    ; For example, a bunch of phases between 330-360 and 0-40 should yield a
    ; mean phase close to 360/0, not 180.  Ditto azimuth.

    phase1_range = phase1_max[0:np,nsum] - phase1_min[0:np,nsum]
    phase2_range = phase2_max[0:np,nsum] - phase2_min[0:np,nsum]
    az1_range = az1_max[0:np,nsum] - az1_min[0:np,nsum]
    az2_range = az2_max[0:np,nsum] - az2_min[0:np,nsum]
    for ch=1,npol do begin
      if sumpol[ch-1] then begin
        if phase1_range[ch-1] gt phase2_range[ch-1] then begin
          newtags[ch-1,nsum].phase = $
                        (phase2_mean[ch-1,nsum]/sumweights[ch-1,nsum]) mod 360.0
        endif else begin
          newtags[ch-1,nsum].phase = phase1_mean[ch-1,nsum]/sumweights[ch-1,nsum]
        endelse
        if az1_range[ch-1] gt az2_range[ch-1] then begin
          newtags[ch-1,nsum].azim = (az2_mean[ch-1,nsum]/sumweights[ch-1,nsum]) mod 360.0
        endif else begin
          newtags[ch-1,nsum].azim = az1_mean[ch-1,nsum]/sumweights[ch-1,nsum]
        endelse
      endif
    endfor

    ; Set the jcp (polarization) tags

    newtags[0:np,nsum].jcp = indgen(npol)+1 ; makes assumptions about what's there.

    ; Get the calculated fractional rms noise deviation

    newtags[0:np,nsum].rmsc = 1/sqrt(rmsc_sum[0:np,nsum])

    ; Get the measured rms noise

    freq = (newtags[0,nsum].posfr) * (newtags[0,nsum].dfreq)  $
                                   * (findgen(newndata[nsum]) - newtags[0,nsum].xjcen)

    for ch=1,npol do begin
      if sumpol[ch-1] then begin
        spec = reform((*wsum[nsum]).spec[ch-1,*])
        in_mask = intarr(newndata[nsum]) + 1
        lim1 = newtags[ch-1,nsum].jsnr1
        lim2 = newtags[ch-1,nsum].jsnr2
        in_mask[lim1:lim2] = 0
        if total(in_mask) ge 2 then begin
          specmean = polymaskfit(freq,spec,0,in_mask=in_mask,yerror=specrms)
          newtags[ch-1,nsum].rmsm = specrms
        endif
      endif
    endfor

    ; If all pairs contributing to this sum had mean Julian date extra tags,
    ; use the new mean Julian date to set some time/date tags

    if jdmeanOK[nsum] and timezoneOK[nsum] then begin
      jd_midnight = floor(jdmean[nsum] - 0.5) + 0.5D
      rctime = round( 86400*((jdmean[nsum] - 0.5D) mod 1) )  ;  nearest sec
      jduse = jd_midnight + rctime/86400.0
      caldat_roundsec,jdmean[nsum],mon,day,year,hr,min,sec
      newtags[0:np,nsum].iyy = year
      newtags[0:np,nsum].imm = mon
      newtags[0:np,nsum].idd = day
      newtags[0:np,nsum].rchour = hr
      newtags[0:np,nsum].rcmin = min
      newtags[0:np,nsum].rcsec = sec
      newtags[0:np,nsum].rcnsec = 0L
      calmean[nsum] = string(year,format='(i4)') + ' '     $
                      + monthnames[mon-1] + ' '            $
                      + string(day,format='(i2.2)') + ' '  $
                      + string(hr,format='(i2.2)') + ':'   $
                      + string(min,format='(i2.2)') + ':'  $
                      + string(sec,format='(i2.2)') + ' '  $
                      + timezone[nsum]
    endif else begin
      newtags[0:np,nsum].iyy = 1L
      newtags[0:np,nsum].imm = 1L  ;  so showstack can call it 'Jan'
      newtags[0:np,nsum].idd = 1L
      if jdmeanOK[nsum] then begin
        print,'WARNING in sumspec: No time information output for sum #', $
              nsum+1,format='(a,i0)'
        print,"                    due to discrepant 'timezone' extra tags"
      endif
    endelse

    ; Make sure there's at least one extra tag (so that the extratags
    ; array is defined): Add the cw tag names if necessary.

    if newnextra[nsum] eq 0 then begin
      newnextra[nsum] = newntags[nsum]
      newextratags[nsum,0:newntags[nsum]-1] = replicate(extratag, ntags[nsum])
      newextratags[nsum,0:newntags[nsum]-1].format = replicate('t', ntags[nsum])
      newextratags[nsum,0:newntags[nsum]-1].name = strlowcase(tag_names(newtags))
      newextratags[nsum,0:newntags[nsum]-1].value = $
                   string(indgen(newntags[nsum]), format='(i0)')
      newextratags[nsum,0:newntags[nsum]-1].comment = ''
    endif

    ; Put it all together and add it to the stack

    if total(sumpol) ge 2 then begin

      ; Pairs

      ; Figure out what the sums' group numbers will be

      maxgroup = max(allgroups)
      newstart = 100*( 1L + maxgroup/100 ) + 1L
      newgroupnumber = lindgen(nsums) + newstart

      stackStruc = {group:newgroupnumber[nsum], $
                    spec:(*wsum[nsum]).spec, $
                    tags:newtags[0:np,nsum], $
                    extratags:newextratags[nsum,0:newnextra[nsum]-1], $ 
                    ndata:newndata[nsum], $
                    ntags:newntags[nsum], $
                    nextra:newnextra[nsum], $
                    tname:newtname[nsum]}
      stack = [stack, ptr_new(stackStruc)]
      nstack = nstack + 1L
      print,'Sum #',nsum+1,' (',n_in[nsum],' pairs) is stack pair #',nstack, $
            ' (group #',newgroupnumber[nsum],')',format='(4(a,i0),a)'

      ; If all pairs contributing to this sum had extra tags for Julian dates,
      ; get rid of any old Julian date extra tags and create new ones
      ;
      ; Ditto for distance, RA, and dec

      if timezoneOK[nsum] then begin
        if jdstart[nsum] gt 0 then $
           setextra,'d','jdstart',string(jdstart[nsum], format='(d13.5)'),stack=nstack
        if jdmeanOK[nsum] then begin
           setextra,'d','jdmean',string(jdmean[nsum], format='(d13.5)'),stack=nstack
           setextra,'s','calmean',calmean[nsum],stack=nstack
        endif
        if jdend[nsum] lt 10000000L then $
           setextra,'d','jdend',string(jdend[nsum], format='(d13.5)'),stack=nstack
      endif
      if dist_min[nsum] gt 0 then $
         setextra,'f','distmin',dist_min[nsum],stack=nstack
      if dist_meanOK[nsum] then $
         setextra,'f','distmean',dist_mean[nsum],stack=nstack
      if dist_max[nsum] lt 100000.0 then $
         setextra,'f','distmax',dist_max[nsum],stack=nstack
      if ra1_min[nsum] ge  0 and ra1_max[nsum] lt 24 and $
         ra2_min[nsum] ge 12 and ra2_max[nsum] lt 36 and ra_meanOK[nsum] then begin
        ra1_range = ra1_max[nsum] - ra1_min[nsum]
        ra2_range = ra2_max[nsum] - ra2_min[nsum]
        if ra2_range lt ra1_range then begin
          setextra,'f','ramin',(ra2_min[nsum] mod 24),stack=nstack
          setextra,'f','ramean',(ra2_mean[nsum] mod 24),stack=nstack
          setextra,'f','ramax',(ra2_max[nsum] mod 24),stack=nstack
        endif else begin
          setextra,'f','ramin',ra1_min[nsum],stack=nstack
          setextra,'f','ramean',ra1_mean[nsum],stack=nstack
          setextra,'f','ramax',ra1_max[nsum],stack=nstack
        endelse
      endif
      if dec_min[nsum] ge -90.0 then $
         setextra,'f','decmin',dec_min[nsum],stack=nstack
      if dec_meanOK[nsum] then $
         setextra,'f','decmean',dec_mean[nsum],stack=nstack
      if dec_max[nsum] le 90.0 then $
         setextra,'f','decmax',dec_max[nsum],stack=nstack
      if keyword_set(noweight) then setextra,'','','',comment='unweighted sum',stack=nstack

    endif else begin

      ; Single channels

      maxgroup = max(stackgroups1())
      newstart = 100*( 1L + maxgroup/100 ) + 1L
      newgroupnumber = lindgen(nsums) + newstart

      stack1Struc = {group:newgroupnumber[nsum], $
                     spec:(*wsum[nsum]).spec[chan-1,*], $
                     tags:newtags[chan-1,nsum], $
                     extratags:newextratags[nsum,0:newnextra[nsum]-1], $ 
                     ndata:newndata[nsum], $
                     ntags:newntags[nsum], $
                     nextra:newnextra[nsum], $
                     tname:newtname[nsum], $
                     pol:polstring}
      if nstack1 eq 0 then begin
        stack1 = ptrarr(1, /allocate_heap)
        *stack1[0] = stack1Struc
      endif else begin
        stack1 = [stack1, ptr_new(stack1Struc)]
      endelse
      nstack1 = nstack1 + 1L
      print,'Sum #',nsum+1,' (',n_in[nsum],' spectra) channel ',chan, $
            ' is stack1 spectrum #',nstack1, $
            ' (group #',newgroupnumber[nsum],')',format='(5(a,i0),a)'
      if timezoneOK[nsum] then begin
        if jdstart[nsum] gt 0 then $
           setextra1,'d','jdstart',string(jdstart[nsum], format='(d13.5)'),stack=nstack1
        if jdmeanOK[nsum] then begin
           setextra1,'d','jdmean',string(jdmean[nsum], format='(d13.5)'),stack=nstack1
           setextra1,'s','calmean',calmean[nsum],stack=nstack1
        endif
        if jdend[nsum] lt 10000000L then $
           setextra1,'d','jdend',string(jdend[nsum], format='(d13.5)'),stack=nstack1
      endif
      if dist_min[nsum] gt 0 then $
         setextra1,'f','distmin',dist_min[nsum],stack=nstack1
      if dist_meanOK[nsum] then $
         setextra1,'f','distmean',dist_mean[nsum],stack=nstack1
      if dist_max[nsum] lt 100000.0 then $
         setextra1,'f','distmax',dist_max[nsum],stack=nstack1
      if ra1_min[nsum] ge  0 and ra1_max[nsum] lt 24 and $
         ra2_min[nsum] ge 12 and ra2_max[nsum] lt 36 and ra_meanOK[nsum] then begin
        ra1_range = ra1_max[nsum] - ra1_min[nsum]
        ra2_range = ra2_max[nsum] - ra2_min[nsum]
        if ra2_range lt ra1_range then begin
          setextra1,'f','ramin',(ra2_min[nsum] mod 24),stack=nstack1
          setextra1,'f','ramean',(ra2_mean[nsum] mod 24),stack=nstack1
          setextra1,'f','ramax',(ra2_max[nsum] mod 24),stack=nstack1
        endif else begin
          setextra1,'f','ramin',ra1_min[nsum],stack=nstack1
          setextra1,'f','ramean',ra1_mean[nsum],stack=nstack1
          setextra1,'f','ramax',ra1_max[nsum],stack=nstack1
        endelse
      endif
      if dec_min[nsum] ge -90.0 then $
         setextra1,'f','decmin',dec_min[nsum],stack=nstack1
      if dec_meanOK[nsum] then $
         setextra1,'f','decmean',dec_mean[nsum],stack=nstack1
      if dec_max[nsum] le 90.0 then $
         setextra1,'f','decmax',dec_max[nsum],stack=nstack1
      if keyword_set(noweight) then setextra1,'','','',comment='unweighted sum',stack=nstack1
    endelse

  endif else begin

    ; Sums with zero or one pair contributing: Just print a message

      print,'Sum #',nsum+1,' (',n_in[nsum],' pairs) will not be output', $
            format='(2(a,i0),a)'

  endelse

endfor

; Clean up pointers

ptr_free,wsum

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sumspec1,group=g,ming=ming,maxg=maxg,arrayg=ag,all=all, $
             noweight=noweight,chan=chan,help=help

; Create one or more weighted spectral sums of stack1 spectra in one or more groups,
; then add the sum(s) to the single-channel stack with new group number(s) attached.
;
; Calling sumspec1 with no arguments creates a separate sum for each group which 
; contains two or more spectra.  Calling sumspec1 with keywords g, ming and maxg,
; or ag set creates a single sum of all spectra in the specified group(s).  Calling 
; "sumspec1,/all" creates a single sum of all spectra in the single-channel stack.
;
; /noweight creates unweighted sums
;
; If the chan keyword is not used, sumspec1 creates separate sums for the OC vs. SC
; polarization channel, over and above the explicit grouping criterion used.  If
; chan is set to 1 or 2, only the OC or SC spectra, respectively, are considered.
;
; 25 Jun 2006: Frequency refers to center of bin, not left edge
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

maxextratags = 1000
too_low_sdev = 1.0d-10

monthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

; Check if a valid set of keywords is set

if n_elements(g) gt 0 then begin
  keywordsOK = n_elements(ming) eq 0 and n_elements(maxg) eq 0 and n_elements(ag) eq 0 $
                                                               and not keyword_set(all)
endif else if n_elements(ming) gt 0 then begin
  keywordsOK = n_elements(maxg) gt 0 and keyword_set(ag) eq 0 and not keyword_set(all)
endif else if n_elements(maxg) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(ag) gt 0 then begin
  keywordsOK = not keyword_set(all)
endif else begin
  keywordsOK = 1
endelse

if n_params() ne 0 or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'sumspec1,group=g,ming=ming,maxg=maxg,arrayg=ag,all=all[,/noweight]', $
        '[,chan=1 or 2][,/help]',format='(2a)'
  print,' '
  print,'    Calling sumspec1,/all creates a single sum of all spectra in the ', $
        'single-channel stack',format='(2a)'
  print,'    Setting keywords g, ming and maxg, or ag yields a single sum of ', $
        'all spectra in the specified group(s)',format='(2a)'
  print,'    Keywords all, g, ming and maxg, and ag are four mutually ', $
        'exclusive choices',format='(2a)'
  print,'    Calling sumspec1 with none of those keywords set creates ', $
        'a sum for each group in the single-channel stack',format='(2a)'
  print,' '
  print,'    /noweight produces unweighted sums'
  print,' '
  print,'    If chan is not specified, separate sums are created for ', $
        'the OC vs. SC',format='(2a)'
  print,'         polarization channel, over and above the explicit grouping ', $
        'criteria used',format='(2a)'
  print,' '
  return
endif

if nstack1 eq 0 then begin
  print,'ERROR in sumspec1: The single-channel stack is empty'
  return
endif

; Check which channels to sum

if n_elements(chan) eq 0 then begin
  sumpol = [1,1]
endif else if chan eq 1 or chan eq 2 then begin
  sumpol = (chan eq 1) ? [1,0] : [0,1]
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse
npols = round(total(sumpol))

; Create an array which contains, for each output sum,
; a sorted vector of the group numbers which will contribute to that sum;
; each sum gets contributions from just one polarization channel.
;
; nsums       = number of sums to do (with OC and SC sums counted separately)
; ngroups     = number of groups to include in each sum (scalar, not vector)
; groupnumber = array of those group numbers:
;                    groupnumber[i,*] = ordered list of group numbers which
;                                       contribute to the ith sum
;                    If both channels are used, these lists are the same for
;                         OC sums (i = 0, 2, 4, ...) and SC sums (1, 3, 5, ...)

allgroups = stackgroups1()
n_all = n_elements(allgroups)

if n_elements(g) gt 0 then begin
  nsums = 1L*npols
  ngroups = 1L
  groupnumber = lonarr(nsums,ngroups)
  for ch=1L,npols do groupnumber[ch-1,0] = long(g)
endif else if n_elements(ming) gt 0 then begin
  if ming gt maxg then begin
    print,'ERROR in sumspec1: Must specify range with ming <= maxg'
    return
  endif
  nsums = 1L*npols
  groupsvec = lonarr(n_all)
  ngroups = 0L
  for k=0L,n_all-1 do begin
    if allgroups[k] ge ming and allgroups[k] le maxg then begin
      groupsvec[ngroups] = allgroups[k]
      ngroups = ngroups + 1L
    endif
  endfor
  if ngroups eq 0 then begin
    print,'sumspec1: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
  groupnumber = lonarr(nsums,ngroups)
  for ch=1L,npols do groupnumber[ch-1,*] = groupsvec[0:ngroups-1]
endif else if n_elements(ag) gt 0 then begin
  nsums = 1L*npols
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupsvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupsvec) lt ngroups then begin
    print,"ERROR in sumspec1: Don't include duplicate group numbers ", $
          "in array ag",format='(2a)'
    return
  endif
  groupnumber = lonarr(nsums,ngroups)
  for ch=1L,npols do groupnumber[ch-1,*] = groupsvec
endif else if keyword_set(all) then begin
  nsums = 1L*npols
  ngroups = n_all
  groupnumber = lonarr(nsums,ngroups)
  for ch=1L,npols do groupnumber[ch-1,*] = allgroups
endif else begin
  nsums = n_all*npols
  ngroups = 1L
  groupnumber = lonarr(nsums,ngroups)
  tempindex = lindgen(n_all)*npols
  for ch=1L,npols do groupnumber[tempindex+ch-1,0] = allgroups
endelse

; Check that all specified groups actually exist

for i=0L,nsums-1,npols do begin
  for j=0L,ngroups-1 do begin
    groupToCheck = groupnumber[i,j]
    groupnum = where(allgroups eq groupToCheck, count)
    if count eq 0 then begin
      print,'ERROR in sumspec1: Group #',groupToCheck, $
            ' is not present in the single-channel stack',format='(a,i0,a)'
      return
    endif
  endfor
endfor

; Initialize some arrays which will be used for constructing the sums

n_in = lonarr(nsums)
newndata = lonarr(nsums)
newtags1 = replicate(blanktags1(), nsums)
newtags1.jsnr1 = 999999999L
newtags1.jsnr2 = -999999999L
newntags = lonarr(nsums)
extratag = {format:'', name:'', value:'', comment:''}
newextratags = replicate(extratag, nsums, maxextratags)
newnextra = lonarr(nsums)
newtname = strarr(nsums)
timezone = strarr(nsums)
timezoneOK = intarr(nsums) + 1L
sumweights = dblarr(nsums)
sumvariance = dblarr(nsums)
rmsc_sum = dblarr(nsums)
phase1_mean = fltarr(nsums)
phase2_mean = fltarr(nsums)
phase1_min = fltarr(nsums) + 999.9
phase1_max = fltarr(nsums) - 999.9
phase2_min = fltarr(nsums) + 999.9
phase2_max = fltarr(nsums) - 999.9
az1_mean = fltarr(nsums)
az2_mean = fltarr(nsums)
az1_min = fltarr(nsums) + 999.9
az1_max = fltarr(nsums) - 999.9
az2_min = fltarr(nsums) + 999.9
az2_max = fltarr(nsums) - 999.9
jdstart = dblarr(nsums) + 999999999L
jdmean = dblarr(nsums)
calmean = strarr(nsums)
jdend = dblarr(nsums) - 999999999L
jdmeanOK = intarr(nsums) + 1L
dist_min = fltarr(nsums) + 999999.9
dist_max = fltarr(nsums) - 999999.9
dist_mean = fltarr(nsums)
dist_meanOK = intarr(nsums) + 1L
ra1_min = fltarr(nsums) + 999.9
ra1_max = fltarr(nsums) - 999.9
ra2_min = fltarr(nsums) + 999.9
ra2_max = fltarr(nsums) - 999.9
ra1_mean = fltarr(nsums)
ra2_mean = fltarr(nsums)
ra_meanOK = intarr(nsums) + 1L
dec_min = fltarr(nsums) + 999.9
dec_max = fltarr(nsums) - 999.9
dec_mean = fltarr(nsums)
dec_meanOK = intarr(nsums) + 1L

; Initialize the polarization tags

newpol = strarr(nsums)
if npols eq 2 then begin
  newpol[2*lindgen(nsums/2)] = 'OC'
  newpol[2*lindgen(nsums/2) + 1] = 'SC'
endif else if chan eq 1 then begin
  newpol[lindgen(nsums)] = 'OC'
endif else begin
  newpol[lindgen(nsums)] = 'SC'
endelse

; Initialize the weighted sums; use pointers in case different sums have
; different numbers of spectral values

wsum = ptrarr(nsums, /allocate_heap)

; Go through the single-channel stack, spectrum by spectrum, and construct the sum(s)

for n=0L,nstack1-1 do begin

  ; See if this spectrum is included in one of the groups which will be summed

  gnum1D = where(groupnumber eq (*stack1[n]).group, count)
  pol = (*stack1[n]).pol
  ch = (pol eq 'OC') ? 1L : 2L

  if count gt 0 and sumpol[ch-1] then begin

    ; This spectrum contributes to a sum: figure out WHICH sum

    gnum1D = gnum1D[(ch-1)*(npols/2)]
    nsum = gnum1D mod nsums
    ngroup = gnum1D/nsums

    ; Get the various elements of this spectrum

    spec = (*stack1[n]).spec
    ndata = (*stack1[n]).ndata
    tags1 = (*stack1[n]).tags
    ntags = (*stack1[n]).ntags
    extratags = (*stack1[n]).extratags
    nextra = (*stack1[n]).nextra
    tname = (*stack1[n]).tname

    ; If this is the first spectrum contributing to this sum, initialize some values

    if n_in[nsum] eq 0 then begin

      ; Create an array which will hold the weighted spectra for this sum
      ;
      ; Since different sums may involve spectra of different lengths, and
      ; since we don't want to presuppose a maximum spectrum length, use a
      ; pointer to reference the array.  But since there doesn't seem to be
      ; any way to do that directly, fool IDL by making a structure whose
      ; only field is the spectral array, then point to that structure.
      ; (Can't yet create the complete output structure: We don't know
      ; how long the extratags array will be.)

      dummyStruc = {spec:fltarr(ndata)}
      *wsum[nsum] = dummyStruc

      ; Use the first spectrum contributing to this sum to get values for various
      ; spectral parameters, so that we can check whether or not later spectra which
      ; contribute to the same sum have the same values

      newndata[nsum] = ndata
      newntags[nsum] = ntags
      newtags1[nsum].xjcen = tags1.xjcen
      newtags1[nsum].posfr = tags1.posfr
      newtags1[nsum].dfreq = tags1.dfreq
      newtags1[nsum].doppl = tags1.doppl
      newtags1[nsum].zepch = tags1.zepch
      newtags1[nsum].obs = tags1.obs
      newtags1[nsum].itar = tags1.itar
      newtags1[nsum].lfft = tags1.lfft
      newtags1[nsum].igw = tags1.igw
      newtags1[nsum].kpts = tags1.kpts
      newtags1[nsum].nfreq = tags1.nfreq
      newtags1[nsum].frstep = tags1.frstep
      newtags1[nsum].color = tags1.color
      newtags1[nsum].freq1 = tags1.freq1
      for k=0L,nextra-1 do begin
        for j=0L,3 do newextratags[nsum,k].(j) = extratags[k].(j)
      endfor
      newnextra[nsum] = nextra
      newtname[nsum] = tname
      timezone[nsum] = getextra1('timezone',stack=(n+1))

    endif

    ; Since this is a weighted sum, and weight = 1/sdev^2,
    ; make sure that sdev isn't zero (or really tiny)

    if not keyword_set(noweight) then begin
      if abs(tags1.sdev) lt too_low_sdev then begin
        print,'ERROR in sumspec1: stack1 spectrum #',n+1,' has sdev < ', $
              too_low_sdev,' km^2',format='(a,i0,a,e7.1,a)'
        return
      endif
    endif

    ; Insist that essential parameters have the same values
    ; for all spectra contributing to this sum

    if tname ne newtname[nsum] then begin
      print,'ERROR in sumspec1: Spectra with different target names contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if ndata ne newndata[nsum] then begin
      print,'ERROR in sumspec1: Spectra with different numbers of spectral points ', $
            'contribute to sum #',nsum+1,format='(2a,i0)'
      return
    endif else if ntags ne newntags[nsum] then begin
      print,'ERROR in sumspec1: Spectra with different numbers of tags contribute to sum #', $
            nsum+1,format='(a,i0)'
      return
    endif else if tags1.xjcen ne newtags1[nsum].xjcen then begin
      print,'ERROR in sumspec1: Spectra with different xjcen tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if tags1.posfr ne newtags1[nsum].posfr then begin
      print,'ERROR in sumspec1: Spectra with different posfr tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if tags1.dfreq ne newtags1[nsum].dfreq then begin
      print,'ERROR in sumspec1: Spectra with different dfreq tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif

    ; Construct the weighted mean spectra
    ;
    ; Note that multiplying the unnormalized spectrum (sdev * spec)
    ; by the weight (1/sdev^2) is equivalent to dividing spec/sdev

    if keyword_set(noweight) then begin
      (*wsum[nsum]).spec = (*wsum[nsum]).spec + float( spec*(1.0D*tags1.sdev) )
    endif else begin
      (*wsum[nsum]).spec = (*wsum[nsum]).spec + float( spec/(1.0D*tags1.sdev) )
    endelse

    ; Accumulate some parameters

    weight = (keyword_set(noweight)) ? 1.0D : 1/(1.0D*tags1.sdev)^2
    sumweights[nsum] = sumweights[nsum] + weight
    rmsc_sum[nsum] = rmsc_sum[nsum] + 1/(1.0D*tags1.rmsc)^2
    newtags1[nsum].nffts = newtags1[nsum].nffts + tags1.nffts
    newtags1[nsum].tau = newtags1[nsum].tau + tags1.tau

    ; Work with two versions of the phase and the azimuth,
    ; using ranges [0,360) vs. [180,540)

    phase1 = tags1.phase - 360*floor(tags1.phase/360.0)
    phase2 = phase1 - 360*floor( (phase1 - 180)/360.0 )
    az1 = tags1.azim - 360*floor(tags1.azim/360.0)
    az2 = az1 - 360*floor( (az1 - 180)/360.0 )

    ; Look for minimum or maximum values for some parameters

    newtags1[nsum].jsnr1 = newtags1[nsum].jsnr1 < tags1.jsnr1
    newtags1[nsum].jsnr2 = newtags1[nsum].jsnr2 > tags1.jsnr2
    jd_in = getextra1('jdstart', stack=(n+1))
    jdstart[nsum] = (notnull(jd_in)) ? (jdstart[nsum] < jd_in) : -999999999L
    jd_in = getextra1('jdend', stack=(n+1))
    jdend[nsum] = (notnull(jd_in)) ? (jdend[nsum] > jd_in) : 999999999L
    phase1_min[nsum] = phase1_min[nsum] < phase1
    phase1_max[nsum] = phase1_max[nsum] > phase1
    phase2_min[nsum] = phase2_min[nsum] < phase2
    phase2_max[nsum] = phase2_max[nsum] > phase2
    az1_min[nsum] = az1_min[nsum] < az1
    az1_max[nsum] = az1_max[nsum] > az1
    az2_min[nsum] = az2_min[nsum] < az2
    az2_max[nsum] = az2_max[nsum] > az2
    dist_in = getextra1('distmin', stack=(n+1))
    dist_min[nsum] = (notnull(dist_in)) ? (dist_min[nsum] < dist_in) : -999999.9
    dist_in = getextra1('distmax', stack=(n+1))
    dist_max[nsum] = (notnull(dist_in)) ? (dist_max[nsum] > dist_in) : 999999.9
    ra1_in = getextra1('ramin', stack=(n+1))
    if notnull(ra1_in) then begin
      ra1_min[nsum] = ra1_min[nsum] < ra1_in
      ra1_max[nsum] = ra1_max[nsum] > ra1_in
      ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
      ra2_min[nsum] = ra2_min[nsum] < ra2_in
      ra2_max[nsum] = ra2_max[nsum] > ra2_in
    endif else begin
      ra1_min[nsum] = -999.9
      ra1_max[nsum] = 999.9
      ra2_min[nsum] = -999.9
      ra2_max[nsum] = 999.9
    endelse
    ra1_in = getextra1('ramax', stack=(n+1))
    if notnull(ra1_in) then begin
      ra1_min[nsum] = ra1_min[nsum] < ra1_in
      ra1_max[nsum] = ra1_max[nsum] > ra1_in
      ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
      ra2_min[nsum] = ra2_min[nsum] < ra2_in
      ra2_max[nsum] = ra2_max[nsum] > ra2_in
    endif else begin
      ra1_min[nsum] = -999.9
      ra1_max[nsum] = 999.9
      ra2_min[nsum] = -999.9
      ra2_max[nsum] = 999.9
    endelse
    dec_in = getextra1('decmin', stack=(n+1))
    dec_min[nsum] = (notnull(dec_in)) ? (dec_min[nsum] < dec_in) : -999.9
    dec_in = getextra1('decmax', stack=(n+1))
    dec_max[nsum] = (notnull(dec_in)) ? (dec_max[nsum] > dec_in) : 999.9

    ; Construct weighted means of some parameters

    newtags1[nsum].elev = newtags1[nsum].elev + weight*tags1.elev
    newtags1[nsum].rttim = newtags1[nsum].rttim + weight*tags1.rttim
    newtags1[nsum].trpwr = newtags1[nsum].trpwr + weight*tags1.trpwr
    newtags1[nsum].tsys = newtags1[nsum].tsys + weight*tags1.tsys
    newtags1[nsum].gain = newtags1[nsum].gain + weight*tags1.gain
    sumvariance[nsum] = sumvariance[nsum] + (1.0D*tags1.sdev)^2
    phase1_mean[nsum] = phase1_mean[nsum] + weight*phase1
    phase2_mean[nsum] = phase2_mean[nsum] + weight*phase2
    az1_mean[nsum] = az1_mean[nsum] + weight*az1
    az2_mean[nsum] = az2_mean[nsum] + weight*az2
    if timezoneOK[nsum] then begin
      tempzone = getextra1('timezone',stack=(n+1))
      timezoneOK[nsum] = (tempzone eq timezone[nsum])
    endif
    if jdmeanOK[nsum] then begin
      jd_in = getextra1('jdmean', stack=(n+1))
      if notnull(jd_in) then begin
        jdmean[nsum] = jdmean[nsum] + weight*jd_in
      endif else begin
        jdmeanOK[nsum] = 0
      endelse
    endif
    if dist_meanOK[nsum] then begin
      dist_in = getextra1('distmean', stack=(n+1))
      if notnull(dist_in) then begin
        dist_mean[nsum] = dist_mean[nsum] + weight*dist_in
      endif else begin
        dist_meanOK[nsum] = 0
      endelse
    endif
    if ra_meanOK[nsum] then begin
      ra1_in = getextra1('ramean', stack=(n+1))
      if notnull(ra1_in) then begin
        ra1_mean[nsum] = ra1_mean[nsum] + weight*ra1_in
        ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
        ra2_mean[nsum] = ra2_mean[nsum] + weight[0]*ra2_in
      endif else begin
        ra_meanOK[nsum] = 0
      endelse
    endif
    if dec_meanOK[nsum] then begin
      dec_in = getextra1('decmean', stack=(n+1))
      if notnull(dec_in) then begin
        dec_mean[nsum] = dec_mean[nsum] + weight*dec_in
      endif else begin
        dec_meanOK[nsum] = 0
      endelse
    endif

    ; Check whether nonessential tags, such as doppl (tx offset), have the same values
    ; for all spectra contributing to this sum.  If not, reset that tag to a dummy value.

    if tags1.doppl ne newtags1[nsum].doppl then newtags1[nsum].doppl = 0.0
    if tags1.zepch ne newtags1[nsum].zepch then newtags1[nsum].zepch = 0.0
    if tags1.obs ne newtags1[nsum].obs then newtags1[nsum].obs = 0L
    if tags1.itar ne newtags1[nsum].itar then newtags1[nsum].itar = 0L
    if tags1.lfft ne newtags1[nsum].lfft then newtags1[nsum].lfft = 0L
    if tags1.igw ne newtags1[nsum].igw then newtags1[nsum].igw = 0L
    if tags1.kpts ne newtags1[nsum].kpts then newtags1[nsum].kpts = 0L
    if tags1.nfreq ne newtags1[nsum].nfreq then newtags1[nsum].nfreq = 0L
    if tags1.frstep ne newtags1[nsum].frstep then newtags1[nsum].frstep = 0.0
    if tags1.color ne newtags1[nsum].color then newtags1[nsum].color = -1L
    if tags1.freq1 ne newtags1[nsum].freq1 then newtags1[nsum].freq1 = 0.0

    ; Go through the array of extra tags and find the ones which are present and have
    ; the same format and value for all spectra contributing to this sum: those are the
    ; ones we'll keep.  Similarly, keep a comment line only if every spectrum has the same
    ; line.  Finally, if formats and values match for a given extra tag but trailing
    ; comments don't match, keep the extra tag but erase the trailing comment.

    for k=newnextra[nsum]-1L,0,-1L do begin

      newextraformat = strlowcase(newextratags[nsum,k].format)
      newextraname = strlowcase(newextratags[nsum,k].name)
      newextravalue = strlowcase(newextratags[nsum,k].value)
      newextracomment = strlowcase(newextratags[nsum,k].comment)
      keepextra = 0
      keepcomment = 0
      extranum = where(extratags.name eq newextraname, extracount)
      if extracount gt 0 then begin
        for j=0L,extracount-1 do begin
          extraformat = strlowcase(extratags[extranum[j]].format)
          extraname = newextraname
          extravalue = strlowcase(extratags[extranum[j]].value)
          extracomment = strlowcase(extratags[extranum[j]].comment)
          match = ((notnull(extraname)) and $
                   (extraformat eq newextraformat) and (extravalue eq newextravalue)) or $
                  ((isnull(extraname)) and (extracomment eq newextracomment))
          keepextra = keepextra or match
          keepcomment = keepcomment or (match and (extracomment eq newextracomment))
        endfor
      endif

      if not keepextra then begin

        ; Don't create an error by trying to remove the last array element

        if newnextra[nsum] gt 1 then begin
          if k eq 0 then begin
            newextratags[nsum,0:newnextra[nsum]-2] = newextratags[nsum,1:newnextra[nsum]-1]
          endif else if k lt newnextra[nsum]-1 then begin
            newextratags[nsum,k:newnextra[nsum]-2] = newextratags[nsum,k+1:newnextra[nsum]-1]
          endif
        endif

        newnextra[nsum] = newnextra[nsum] - 1

      endif else if not keepcomment then begin
        newextratags[nsum,k].comment = ''
      endif

    endfor

    ; Increment the number of spectra contributing to this sum

    n_in[nsum] = n_in[nsum] + 1L

  endif  ;  count gt 0 and sumpol[ch-1] ?

; Look at the next spectrum in the single-channel stack

endfor

; All stack1 spectra have been inspected

for nsum=0L,nsums-1 do begin

  if n_in[nsum] ge 2 then begin

    ; Normalize the weighted mean spectra and get the rms noise (sdev)
    ;
    ; Note that dividing by the sum of all weights, then normalizing
    ; by dividing by sdev = 1/sqrt(sum of all weights), is equivalent
    ; to multiplying by sdev

    if keyword_set(noweight) then begin
      sdev = sqrt(sumvariance[nsum])
      (*wsum[nsum]).spec = float( (*wsum[nsum]).spec / sdev )
      newtags1[nsum].sdev = sdev/n_in[nsum]
    endif else begin
      sdev = 1.0/sqrt(sumweights[nsum])
      (*wsum[nsum]).spec = sdev * (*wsum[nsum]).spec
      newtags1[nsum].sdev = sdev
    endelse

    ; Complete the weighted means for various tags

    newtags1[nsum].elev = newtags1[nsum].elev / sumweights[nsum]
    newtags1[nsum].rttim = newtags1[nsum].rttim / sumweights[nsum]
    newtags1[nsum].trpwr = newtags1[nsum].trpwr / sumweights[nsum]
    newtags1[nsum].tsys = newtags1[nsum].tsys / sumweights[nsum]
    newtags1[nsum].gain = newtags1[nsum].gain / sumweights[nsum]
    if jdmeanOK[nsum] then jdmean[nsum] = jdmean[nsum] / sumweights[nsum]
    if dist_meanOK[nsum] then dist_mean[nsum] = dist_mean[nsum] / sumweights[nsum]
    if ra_meanOK[nsum] then begin
      ra1_mean[nsum] = ra1_mean[nsum] / sumweights[nsum]
      ra2_mean[nsum] = ra2_mean[nsum] / sumweights[nsum]
    endif
    if dec_meanOK[nsum] then dec_mean[nsum] = dec_mean[nsum] / sumweights[nsum]

    ; Deal with rotation phase and azimuth:
    ;
    ; If the spectra contributing to this sum have phases which are more
    ; "concentrated" when considered over the range [180,540) rather than [0,360),
    ; then compute the mean phase using the former range (and then take mod 360).
    ; For example, a bunch of phases between 330-360 and 0-40 should yield a
    ; mean phase close to 360/0, not 180.  Ditto azimuth.

    phase1_range = phase1_max[nsum] - phase1_min[nsum]
    phase2_range = phase2_max[nsum] - phase2_min[nsum]
    az1_range = az1_max[nsum] - az1_min[nsum]
    az2_range = az2_max[nsum] - az2_min[nsum]
    if phase1_range gt phase2_range then begin
      newtags1[nsum].phase = (phase2_mean[nsum]/sumweights[nsum]) mod 360.0
    endif else begin
      newtags1[nsum].phase = phase1_mean[nsum]/sumweights[nsum]
    endelse
    if az1_range gt az2_range then begin
      newtags1[nsum].azim = (az2_mean[nsum]/sumweights[nsum]) mod 360.0
    endif else begin
      newtags1[nsum].azim = az1_mean[nsum]/sumweights[nsum]
    endelse

    ; Set the jcp (polarization) tags

    ch = (newpol[nsum] eq 'OC') ? 1L : 2L
    newtags1[nsum].jcp = ch

    ; Get the calculated fractional rms noise deviation

    newtags1[nsum].rmsc = 1/sqrt(rmsc_sum[nsum])

    ; Get the measured rms noise

    freq = (newtags1[nsum].posfr) * (newtags1[nsum].dfreq)  $
                                  * (findgen(newndata[nsum]) - newtags1[nsum].xjcen)
    spec = (*wsum[nsum]).spec
    in_mask = intarr(newndata[nsum]) + 1
    lim1 = newtags1[nsum].jsnr1
    lim2 = newtags1[nsum].jsnr2
    in_mask[lim1:lim2] = 0
    if total(in_mask) ge 2 then begin
      specmean = polymaskfit(freq,spec,0,in_mask=in_mask,yerror=specrms)
      newtags1[nsum].rmsm = specrms
    endif

    ; If all spectra contributing to this sum had mean Julian date extra tags,
    ; use the new mean Julian date to set some time/date tags

    if jdmeanOK[nsum] and timezoneOK[nsum] then begin
      jd_midnight = floor(jdmean[nsum] - 0.5) + 0.5D
      rctime = round( 86400*((jdmean[nsum] - 0.5D) mod 1) )  ;  nearest sec
      jduse = jd_midnight + rctime/86400.0
      caldat_roundsec,jdmean[nsum],mon,day,year,hr,min,sec
      newtags1[nsum].iyy = year
      newtags1[nsum].imm = mon
      newtags1[nsum].idd = day
      newtags1[nsum].rchour = hr
      newtags1[nsum].rcmin = min
      newtags1[nsum].rcsec = sec
      newtags1[nsum].rcnsec = 0L
      calmean[nsum] = string(year,format='(i4)') + ' '     $
                      + monthnames[mon-1] + ' '            $
                      + string(day,format='(i2.2)') + ' '  $
                      + string(hr,format='(i2.2)') + ':'   $
                      + string(min,format='(i2.2)') + ':'  $
                      + string(sec,format='(i2.2)') + ' '  $
                      + timezone[nsum]
    endif else begin
      newtags1[nsum].iyy = 1L
      newtags1[nsum].imm = 1L  ;  so showstack can call it 'Jan'
      newtags1[nsum].idd = 1L
      if jdmeanOK[nsum] then begin
        print,'WARNING in sumspec1: No time information output for sum #', $
              nsum+1,format='(a,i0)'
        print,"                     due to discrepant 'timezone' extra tags"
      endif
    endelse

    ; Make sure there's at least one extra tag (so that the extratags
    ; array is defined): Add the cw tag names if necessary.

    if newnextra[nsum] eq 0 then begin
      newnextra[nsum] = newntags[nsum]
      newextratags[nsum,0:newntags[nsum]-1] = replicate(extratag, ntags[nsum])
      newextratags[nsum,0:newntags[nsum]-1].format = replicate('t', ntags[nsum])
      newextratags[nsum,0:newntags[nsum]-1].name = strlowcase(tag_names(newtags))
      newextratags[nsum,0:newntags[nsum]-1].value = $
                   string(indgen(newntags[nsum]), format='(i0)')
      newextratags[nsum,0:newntags[nsum]-1].comment = ''
    endif

    ; Put it all together and add it to the single-channel stack

    ; Figure out what the sums' group numbers will be

    maxgroup = max(allgroups)
    newstart = 100*( 1L + maxgroup/100 ) + 1L
    newgroupnumber = lindgen(nsums) + newstart

    stack1Struc = {group:newgroupnumber[nsum], $
                   spec:(*wsum[nsum]).spec, $
                   tags:newtags1[nsum], $
                   extratags:newextratags[nsum,0:newnextra[nsum]-1], $ 
                   ndata:newndata[nsum], $
                   ntags:newntags[nsum], $
                   nextra:newnextra[nsum], $
                   tname:newtname[nsum], $
                   pol:newpol[nsum]}
    stack1 = [stack1, ptr_new(stack1Struc)]
    nstack1 = nstack1 + 1L
    print,'Sum #',nsum+1,' (',n_in[nsum],' ',newpol[nsum], $
          ' spectra) is stack1 spectrum #',nstack1,' (group #', $
          newgroupnumber[nsum],')',format='(a,i0,a,i0,3a,i0,a,i0,a)'

    ; If all spectra contributing to this sum had extra tags for Julian dates,
    ; get rid of any old Julian date extra tags and create new ones
    ;
    ; Ditto for distance, RA, and dec

    if timezoneOK[nsum] then begin
      if jdstart[nsum] gt 0 then $
         setextra1,'d','jdstart',string(jdstart[nsum], format='(d13.5)'),stack=nstack1
      if jdmeanOK[nsum] then begin
         setextra1,'d','jdmean',string(jdmean[nsum], format='(d13.5)'),stack=nstack1
         setextra1,'s','calmean',calmean[nsum],stack=nstack1
      endif
      if jdend[nsum] lt 10000000L then $
         setextra1,'d','jdend',string(jdend[nsum], format='(d13.5)'),stack=nstack1
    endif
    if dist_min[nsum] gt 0 then $
       setextra1,'f','distmin',dist_min[nsum],stack=nstack1
    if dist_meanOK[nsum] then $
       setextra1,'f','distmean',dist_mean[nsum],stack=nstack1
    if dist_max[nsum] lt 100000.0 then $
       setextra1,'f','distmax',dist_max[nsum],stack=nstack1
    if ra1_min[nsum] ge  0 and ra1_max[nsum] lt 24 and $
       ra2_min[nsum] ge 12 and ra2_max[nsum] lt 36 and ra_meanOK[nsum] then begin
      ra1_range = ra1_max[nsum] - ra1_min[nsum]
      ra2_range = ra2_max[nsum] - ra2_min[nsum]
      if ra2_range lt ra1_range then begin
        setextra1,'f','ramin',(ra2_min[nsum] mod 24),stack=nstack1
        setextra1,'f','ramean',(ra2_mean[nsum] mod 24),stack=nstack1
        setextra1,'f','ramax',(ra2_max[nsum] mod 24),stack=nstack1
      endif else begin
        setextra1,'f','ramin',ra1_min[nsum],stack=nstack1
        setextra1,'f','ramean',ra1_mean[nsum],stack=nstack1
        setextra1,'f','ramax',ra1_max[nsum],stack=nstack1
      endelse
    endif
    if dec_min[nsum] ge -90.0 then $
       setextra1,'f','decmin',dec_min[nsum],stack=nstack1
    if dec_meanOK[nsum] then $
       setextra1,'f','decmean',dec_mean[nsum],stack=nstack1
    if dec_max[nsum] le 90.0 then $
       setextra1,'f','decmax',dec_max[nsum],stack=nstack1
    if keyword_set(noweight) then setextra1,'','','',comment='unweighted sum',stack=nstack1

  endif else begin

    ; Sums with zero or one spectrum contributing: Just print a message

    print,'Sum #',nsum+1,' (',n_in[nsum],' ',newpol[nsum], $
          ' spectra) will not be output',format='(a,i0,a,i0,3a)'

  endelse

endfor

; Clean up pointers

ptr_free,wsum

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sumim,group=g,ming=ming,maxg=maxg,arrayg=ag,all=all, $
    eph_row_new=eph_row_new,interp=interp,cubic=cubic,   $
    noweight=noweight,calibrated=calibrated,chan=chan,help=help

; Create one or more weighted sums of stacki images in one or more groups,
; then add the sum(s) to the image stack with new group number(s) attached.
;
; Calling sumim with no arguments creates a separate sum for each group which 
; contains two or more images.  Calling sumim with keywords g, ming and maxg,
; or ag set creates a single sum of all images in the specified group(s).  Calling 
; "sumim,/all" creates a single sum of all images in the image stack.
;
; If the chan keyword is not used, sumim creates separate sums for the OC vs. SC
; polarization channel, over and above the explicit grouping criterion used.  If
; chan is set to 1 or 2, only the OC or SC images, respectively, are considered.
;
; Modified 2009 Sep 7 by CM:
;     When /interp is set, correct for reduced r.m.s. noise that results from
;         interpolation, so that the output sums are properly normalized if the
;         input images are properly normalized
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

maxextratags = 1000
too_low_sdev = 1.0d-10

monthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

; Check if a valid set of keywords is set

if n_elements(g) gt 0 then begin
  keywordsOK = n_elements(ming) eq 0 and n_elements(maxg) eq 0 and n_elements(ag) eq 0 $
                                                               and not keyword_set(all)
endif else if n_elements(ming) gt 0 then begin
  keywordsOK = n_elements(maxg) gt 0 and keyword_set(ag) eq 0 and not keyword_set(all)
endif else if n_elements(maxg) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(ag) gt 0 then begin
  keywordsOK = not keyword_set(all)
endif else begin
  keywordsOK = 1
endelse
if not keyword_set(interp) and n_elements(cubic) gt 0 then keywordsOK = 0
if not keyword_set(noweight) and keyword_set(calibrated) then keywordsOK = 0

if n_params() ne 0 or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'sumim,group=g,ming=ming,maxg=maxg,arrayg=ag,all=all     $'
  print,'     [,eph_row_new=eph_row_new][,/interp[,cubic=cubic]] $'
  print,'     [,/noweight[,/calibrated]][,chan=1 or 2][,/help]'
  print,' '
  print,'    Calling sumim,/all creates a single sum of all images in the ', $
        'image stack',format='(2a)'
  print,'    Setting keywords g, ming and maxg, or ag yields a single sum of ', $
        'all images in the specified group(s)',format='(2a)'
  print,'    Keywords all, g, ming and maxg, and ag are four mutually ', $
        'exclusive choices',format='(2a)'
  print,'    Calling sumim with none of those keywords set creates ', $
        'a sum for each group in the image stack',format='(2a)'
  print,' '
  print,'    eph_row_new is the eph_row value for each sum; if omitted,'
  print,'        the eph_row value for the first image in each sum is used'
  print,'        as the eph_row value for the sum'
  print,' '
  print,'    If neither /interp nor cubic is used, each image is shifted by an'
  print,'        integer number of delay rows so that its eph_row value is within'
  print,'        half a row of the eph_row value for that sum'
  print,"    If /interp is set, linear interpolation (IDL's 'interpolate' routine)"
  print,'        is then applied to the shifted image to match eph_row exactly'
  print,'    If /interp is set and the cubic keyword is used, cubic convolution'
  print,"        rather than linear interpolation is applied; the value of 'cubic'"
  print,'        can be between -1.0 and 0.0, with any positive value reset to -1.0'
  print,"        (see documentation for IDL's 'interpolate' routine)"
  print,' '
  print,'    /noweight produces an unweighted sum.  If /calibrated is NOT set,'
  print,"        input images needn't have sdev tags, and no sdev tag is output"
  print,' '
  print,'    /calibrated assumes that all input images are calibrated and uses their'
  print,'         sdev values to calculate and output sdev for the unweighted sum'
  print,' '
  print,'    If chan is not specified, separate sums are created for ', $
        'the OC vs. SC',format='(2a)'
  print,'         polarization channel, over and above the explicit grouping ', $
        'criteria used',format='(2a)'
  print,' '
  return
endif

if nstacki eq 0 then begin
  print,'ERROR in sumim: The image stack is empty'
  return
endif

; Check which channels to sum

if n_elements(chan) eq 0 then begin
  sumpol = [1,1]
endif else if chan eq 1 or chan eq 2 then begin
  sumpol = (chan eq 1) ? [1,0] : [0,1]
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse
npols = round(total(sumpol))

; Create an array which contains, for each output sum,
; a sorted vector of the group numbers which will contribute to that sum;
; each sum gets contributions from just one polarization channel.
;
; nsums       = number of sums to do (with OC and SC sums counted separately)
; ngroups     = number of groups to include in each sum (scalar, not vector)
; groupnumber = array of those group numbers:
;                    groupnumber[i,*] = ordered list of group numbers which
;                                       contribute to the ith sum
;                    If both channels are used, these lists are the same for
;                         OC sums (i = 0, 2, 4, ...) and SC sums (1, 3, 5, ...)

allgroups = stackgroupsi()
n_all = n_elements(allgroups)

if n_elements(g) gt 0 then begin
  nsums = 1L*npols
  ngroups = 1L
  groupnumber = lonarr(nsums,ngroups)
  for ch=1L,npols do groupnumber[ch-1,0] = long(g)
endif else if n_elements(ming) gt 0 then begin
  if ming gt maxg then begin
    print,'ERROR in sumim: Must specify range with ming <= maxg'
    return
  endif
  nsums = 1L*npols
  groupsvec = lonarr(n_all)
  ngroups = 0L
  for k=0L,n_all-1 do begin
    if allgroups[k] ge ming and allgroups[k] le maxg then begin
      groupsvec[ngroups] = allgroups[k]
      ngroups = ngroups + 1L
    endif
  endfor
  if ngroups eq 0 then begin
    print,'sumim: There are no groups in the specified range (', $
          ming,'-',maxg,')',format='(a,i0,a,i0,a)'
    return
  endif
  groupnumber = lonarr(nsums,ngroups)
  for ch=1L,npols do groupnumber[ch-1,*] = groupsvec[0:ngroups-1]
endif else if n_elements(ag) gt 0 then begin
  nsums = 1L*npols
  ngroups = n_elements(ag)
  longarray = long(ag)
  groupsvec = longarray[ uniq(longarray, sort(longarray)) ]
  if n_elements(groupsvec) lt ngroups then begin
    print,"ERROR in sumim: Don't include duplicate group numbers ", $
          "in array ag",format='(2a)'
    return
  endif
  groupnumber = lonarr(nsums,ngroups)
  for ch=1L,npols do groupnumber[ch-1,*] = groupsvec
endif else if keyword_set(all) then begin
  nsums = 1L*npols
  ngroups = n_all
  groupnumber = lonarr(nsums,ngroups)
  for ch=1L,npols do groupnumber[ch-1,*] = allgroups
endif else begin
  nsums = n_all*npols
  ngroups = 1L
  groupnumber = lonarr(nsums,ngroups)
  tempindex = lindgen(n_all)*npols
  for ch=1L,npols do groupnumber[tempindex+ch-1,0] = allgroups
endelse

; Check that all specified groups actually exist

for i=0L,nsums-1,npols do begin
  for j=0L,ngroups-1 do begin
    groupToCheck = groupnumber[i,j]
    groupnum = where(allgroups eq groupToCheck, count)
    if count eq 0 then begin
      print,'ERROR in sumim: Group #',groupToCheck, $
            ' is not present in the image stack',format='(a,i0,a)'
      return
    endif
  endfor
endfor

; Initialize some arrays which will be used for constructing the sums

n_in = lonarr(nsums)
newwidth = lonarr(nsums)
newheight = lonarr(nsums)
newdfreq = fltarr(nsums)
neweph_row = fltarr(nsums)
neweph_col = fltarr(nsums)
newdelayunit = fltarr(nsums)
newlambda = fltarr(nsums)
newelev = fltarr(nsums)
newrtt = fltarr(nsums)
newtxpower = fltarr(nsums)
newtsys = fltarr(nsums)
newgain = fltarr(nsums)
newnlooks = lonarr(nsums)
newtau = fltarr(nsums)
newdelbin1 = lonarr(nsums) + 999999999L
newdelbin2 = lonarr(nsums) - 999999999L
newdopbin1 = lonarr(nsums) + 999999999L
newdopbin2 = lonarr(nsums) - 999999999L
elevOK = intarr(nsums) + 1
azimOK = intarr(nsums) + 1
rttOK = intarr(nsums) + 1
txpowerOK = intarr(nsums) + 1
tsysOK = intarr(nsums) + 1
gainOK = intarr(nsums) + 1
nlooksOK = intarr(nsums) + 1
tauOK = intarr(nsums) + 1
delbinOK = intarr(nsums) + 1
dopbinOK = intarr(nsums) + 1
phaseOK = intarr(nsums) + 1
extratag = {format:'', name:'', value:'', comment:''}
newextratags = replicate(extratag, nsums, maxextratags)
newnextra = lonarr(nsums)
newtname = strarr(nsums)
timezone = strarr(nsums)
timezoneOK = intarr(nsums) + 1
sumweights = dblarr(nsums)
sumvariance = dblarr(nsums)
phase1_mean = fltarr(nsums)
phase2_mean = fltarr(nsums)
phase1_min = fltarr(nsums) + 999.9
phase1_max = fltarr(nsums) - 999.9
phase2_min = fltarr(nsums) + 999.9
phase2_max = fltarr(nsums) - 999.9
az1_mean = fltarr(nsums)
az2_mean = fltarr(nsums)
az1_min = fltarr(nsums) + 999.9
az1_max = fltarr(nsums) - 999.9
az2_min = fltarr(nsums) + 999.9
az2_max = fltarr(nsums) - 999.9
jdstart = dblarr(nsums) + 999999999L
jdmean = dblarr(nsums)
calmean = strarr(nsums)
jdend = dblarr(nsums) - 999999999L
jdmeanOK = intarr(nsums) + 1
dist_min = fltarr(nsums) + 999999.9
dist_max = fltarr(nsums) - 999999.9
dist_mean = fltarr(nsums)
dist_meanOK = intarr(nsums) + 1
ra1_min = fltarr(nsums) + 999.9
ra1_max = fltarr(nsums) - 999.9
ra2_min = fltarr(nsums) + 999.9
ra2_max = fltarr(nsums) - 999.9
ra1_mean = fltarr(nsums)
ra2_mean = fltarr(nsums)
ra_meanOK = intarr(nsums) + 1
dec_min = fltarr(nsums) + 999.9
dec_max = fltarr(nsums) - 999.9
dec_mean = fltarr(nsums)
dec_meanOK = intarr(nsums) + 1

; Initialize the polarization tags

newpol = strarr(nsums)
if npols eq 2 then begin
  newpol[2*lindgen(nsums/2)] = 'OC'
  newpol[2*lindgen(nsums/2) + 1] = 'SC'
endif else if chan eq 1 then begin
  newpol[lindgen(nsums)] = 'OC'
endif else begin
  newpol[lindgen(nsums)] = 'SC'
endelse

; Initialize the weighted sums; use pointers in case different sums have
; different image dimensions

wsum = ptrarr(nsums, /allocate_heap)

; Go through the image stack, image by image, and construct the sum(s)

for n=0L,nstacki-1 do begin

  ; See if this image is included in one of the groups which will be summed

  gnum1D = where(groupnumber eq (*stacki[n]).group, count)
  pol = (*stacki[n]).pol
  ch = (pol eq 'OC') ? 1L : 2L

  if count gt 0 and sumpol[ch-1] then begin

    ; This image contributes to a sum: figure out WHICH sum

    gnum1D = gnum1D[(ch-1)*(npols/2)]
    nsum = gnum1D mod nsums
    ngroup = gnum1D/nsums

    ; Get the various elements of this image

    image = (*stacki[n]).image
    width = (*stacki[n]).width
    height = (*stacki[n]).height
    extratags = (*stacki[n]).extratags
    nextra = (*stacki[n]).nextra
    tname = (*stacki[n]).tname
    eph_row = getextrai('eph_row',stack=(n+1))
    eph_col = getextrai('eph_col',stack=(n+1))
    dfreq = getextrai('fres',stack=(n+1))
    delayunit = getextrai('delayunit',stack=(n+1))
    spb = getextrai('samples_per_baud',stack=(n+1))
    rows_per_baud = getextrai('rows_per_baud',stack=(n+1))
    baud = getextrai('baudlen',stack=(n+1))
    codemethod = getextrai('codemethod',stack=(n+1))
    if isnull(spb) then spb = 1L
    if isnull(rows_per_baud) then rows_per_baud = spb
    stride = spb/rows_per_baud
    if isnull(delayunit) and notnull(baud) then delayunit = baud/rows_per_baud
    if isnull(codemethod) then codemethod = 'short'
    if isnull(eph_col) or isnull(eph_row) or isnull(dfreq) $
                       or isnull(delayunit) then begin
      print,"ERROR in sumim: stacki image #",n+1, $
            " needs the 'eph_col' and 'eph_row' and 'fres' tags,",format='(2a)'
      print,"                    plus either the 'delayunit' or 'baudlen' tag"
      return
    endif
    delbin1 = getextrai('delbin1',stack=(n+1))
    delbin2 = getextrai('delbin2',stack=(n+1))
    dopbin1 = getextrai('dopbin1',stack=(n+1))
    dopbin2 = getextrai('dopbin2',stack=(n+1))
    azim = getextrai('azimuth',stack=(n+1))
    elev = getextrai('elevation',stack=(n+1))
    rtt = getextrai('rtt',stack=(n+1))
    txpower = getextrai('tx_power',stack=(n+1))
    tsys = getextrai('tsys',stack=(n+1))
    gain = getextrai('gain_rxtx',stack=(n+1))
    phase = getextrai('phase',stack=(n+1))
    nlooks = getextrai('nlooks',stack=(n+1))
    tau = getextrai('int_time',stack=(n+1))
    lambda = getextrai('lambda',stack=(n+1))
    sdev = getextrai('sdev',stack=(n+1))
    if isnull(sdev) and $
           not (keyword_set(noweight) and not keyword_set(calibrated)) then begin
      print,'ERROR in sumim: stacki image #',n+1," needs the 'sdev' tag", $
            format='(a,i0,a)'
      print,'                -- must set /noweight (and NOT set /calibrated)'
      return
    endif

    ; Check whether or not various nonessential extra tags are present

    if isnull(delbin1) or isnull(delbin2) then delbinOK[nsum] = 0
    if isnull(dopbin1) or isnull(dopbin2) then dopbinOK[nsum] = 0
    if isnull(azim) then azimOK[nsum] = 0
    if isnull(elev) then elevOK[nsum] = 0
    if isnull(rtt) then rttOK[nsum] = 0
    if isnull(txpower) then txpowerOK[nsum] = 0
    if isnull(tsys) then tsysOK[nsum] = 0
    if isnull(gain) then gainOK[nsum] = 0
    if isnull(phase) then phaseOK[nsum] = 0
    if isnull(nlooks) then nlooksOK[nsum] = 0
    if isnull(tau) then tauOK[nsum] = 0

    ; If this is the first image contributing to this sum, initialize some values

    if n_in[nsum] eq 0 then begin

      ; Create an array which will hold the weighted images for this sum
      ;
      ; Since different sums may involve images of different dimensions, and
      ; since we don't want to presuppose a maximum image size, use a
      ; pointer to reference the array.  But since there doesn't seem to be
      ; any way to do that directly, fool IDL by making a structure whose
      ; only field is the image array, then point to that structure.
      ; (Can't yet create the complete output structure: We don't know
      ; how long the extratags array will be.)

      dummyStruc = {image:fltarr(width, height)}
      *wsum[nsum] = dummyStruc

      ; Use the first image contributing to this sum to get the value
      ; of eph_row to which all images in this sum will be interpolated

      neweph_row[nsum] = (keyword_set(eph_row_new)) ? eph_row_new : eph_row

      ; Also use this first image in the sum to get values for various
      ; image parameters, so that we can check whether or not later images which
      ; contribute to the same sum have the same values

      newwidth[nsum] = width
      newheight[nsum] = height
      newdfreq[nsum] = dfreq
      neweph_col[nsum] = eph_col
      newdelayunit[nsum] = delayunit
      newlambda[nsum] = lambda
      for k=0L,nextra-1 do begin
        for j=0L,3 do newextratags[nsum,k].(j) = extratags[k].(j)
        if strlowcase(extratags[k].(1)) eq 'eph_row' then begin
          newextratags[nsum,k].(2) = string(neweph_row[nsum],format='(f0)')
        endif
      endfor
      newnextra[nsum] = nextra
      newtname[nsum] = tname
      timezone[nsum] = getextrai('timezone',stack=(n+1))

    endif

    ; Since this is a weighted sum, and weight = 1/sdev^2,
    ; make sure that sdev isn't zero (or really tiny)

    if not keyword_set(noweight) then begin
      if abs(sdev) lt too_low_sdev then begin
        print,'ERROR in sumim: stacki image #',n+1,' has sdev < ', $
              too_low_sdev,' km^2',format='(a,i0,a,e7.1,a)'
        return
      endif
    endif

    ; Insist that essential parameters (EXCEPT eph_row) have the same values
    ; for all images contributing to this sum

    if tname ne newtname[nsum] then begin
      print,'ERROR in sumim: Images with different target names contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if width ne newwidth[nsum] then begin
      print,'ERROR in sumim: Images with different widths ', $
            'contribute to sum #',nsum+1,format='(2a,i0)'
      return
    endif else if height ne newheight[nsum] then begin
      print,'ERROR in sumim: Images with different heights ', $
            'contribute to sum #',nsum+1,format='(2a,i0)'
      return
    endif else if eph_col ne neweph_col[nsum] then begin
      print,'ERROR in sumim: Images with different eph_col tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if delayunit ne newdelayunit[nsum] then begin
      print,'ERROR in sumim: Images with different delayunit tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if dfreq ne newdfreq[nsum] then begin
      print,'ERROR in sumim: Images with different fres tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif else if lambda ne newlambda[nsum] then begin
      print,'ERROR in sumim: Images with different lambda tags contribute to sum #',nsum+1, $
            format='(a,i0)'
      return
    endif

    ; Shift the image by an integer number of delay rows so that
    ; eph_row is within half a row of the desired value.
    ;
    ; Note that the IDL "shift" function does a circular shift,
    ; so we'll change that to replace pixel values with zero.

    nshift = round(neweph_row[nsum] - eph_row)
    if nshift ne 0 then begin
      image = shift(image,0,nshift)
      if nshift gt 0 then begin
        image[*,0:nshift-1] = 0.0
      endif else begin
        image[*,height+nshift:height-1] = 0.0
      endelse
      eph_row = eph_row + nshift
      if delbinOK[nsum] then begin
        delbin1 = (delbin1 + nshift) > 0L
        delbin2 = (delbin2 + nshift) < (height - 1)
      endif
    endif

    ; If desired, interpolate the image in delay in order to get all
    ; eph_row values to match exactly for each sum, and then correct
    ; for the r.m.s. noise reduction produced by this interpolation
    ; ("sdevfactor" = noise reduction factor: 0.0 < sdevfactor <= 1.0)

    if keyword_set(interp) then begin
      fshift = neweph_row[nsum] - eph_row
      if fshift ne 0.0 then begin
        if n_elements(cubic) gt 0 then begin
          image = interpolate(image, findgen(width), findgen(height)-fshift, $
                              /grid, missing=0.0, cubic=cubic)
          sdevfactor = sdevfactor_cubicconv(fshift,0.0,codemethod,spb,stride, $
                                            cubic=cubic)
        endif else begin
          image = interpolate(image, findgen(width), findgen(height)-fshift, $
                              /grid, missing=0.0)
          sdevfactor = sdevfactor_bilinear(fshift,0.0,codemethod,spb,stride)
        endelse
        image = image/sdevfactor
        if not keyword_set(noweight) then sdev = sdev*sdevfactor
      endif
      eph_row = eph_row + fshift
      if delbinOK[nsum] then begin
        delbin1 = floor(delbin1 + fshift) > 0L
        delbin2 = ceil(delbin2 + fshift) < (height - 1)
      endif
    endif

    ; Construct the weighted mean images
    ;
    ; Note that multiplying the unnormalized image (sdev * image)
    ; by the weight (1/sdev^2) is equivalent to dividing image/sdev

    if keyword_set(noweight) then begin
      if keyword_set(calibrated) then begin
        (*wsum[nsum]).image = (*wsum[nsum]).image + float( image*(1.0D*sdev) )
      endif else begin
        (*wsum[nsum]).image = (*wsum[nsum]).image + image
      endelse
    endif else begin
      (*wsum[nsum]).image = (*wsum[nsum]).image + float( image/(1.0D*sdev) )
    endelse

    ; Accumulate some parameters

    if keyword_set(noweight) then begin
      weight = 1.0D
      sumweights[nsum] = sumweights[nsum] + weight
      if keyword_set(calibrated) then $
           sumvariance[nsum] = sumvariance[nsum] + (1.0D*sdev)^2
    endif else begin
      weight = 1/(1.0D*sdev)^2
      sumweights[nsum] = sumweights[nsum] + weight
    endelse
    if nlooksOK[nsum] then newnlooks[nsum] = newnlooks[nsum] + nlooks
    if tauOK[nsum] then newtau[nsum] = newtau[nsum] + tau

    ; Work with two versions of the phase and the azimuth,
    ; using ranges [0,360) vs. [180,540)

    if phaseOK[nsum] then begin
      phase1 = phase - 360*floor(phase/360.0)
      phase2 = phase1 - 360*floor( (phase1 - 180)/360.0 )
    endif
    if azimOK[nsum] then begin
      az1 = azim - 360*floor(azim/360.0)
      az2 = az1 - 360*floor( (az1 - 180)/360.0 )
    endif

    ; Look for minimum or maximum values for some parameters

    if delbinOK[nsum] then begin
      newdelbin1[nsum] = newdelbin1[nsum] < delbin1
      newdelbin2[nsum] = newdelbin2[nsum] > delbin2
    endif
    if dopbinOK[nsum] then begin
      newdopbin1[nsum] = newdopbin1[nsum] < dopbin1
      newdopbin2[nsum] = newdopbin2[nsum] > dopbin2
    endif
    jd_in = getextrai('jdstart', stack=(n+1))
    jdstart[nsum] = (notnull(jd_in)) ? (jdstart[nsum] < jd_in) : -999999999L
    jd_in = getextrai('jdend', stack=(n+1))
    jdend[nsum] = (notnull(jd_in)) ? (jdend[nsum] > jd_in) : 999999999L
    if phaseOK[nsum] then begin
      phase1_min[nsum] = phase1_min[nsum] < phase1
      phase1_max[nsum] = phase1_max[nsum] > phase1
      phase2_min[nsum] = phase2_min[nsum] < phase2
      phase2_max[nsum] = phase2_max[nsum] > phase2
    endif
    if azimOK[nsum] then begin
      az1_min[nsum] = az1_min[nsum] < az1
      az1_max[nsum] = az1_max[nsum] > az1
      az2_min[nsum] = az2_min[nsum] < az2
      az2_max[nsum] = az2_max[nsum] > az2
    endif
    dist_in = getextrai('distmin', stack=(n+1))
    dist_min[nsum] = (notnull(dist_in)) ? (dist_min[nsum] < dist_in) : -999999.9
    dist_in = getextrai('distmax', stack=(n+1))
    dist_max[nsum] = (notnull(dist_in)) ? (dist_max[nsum] > dist_in) : 999999.9
    ra1_in = getextrai('ramin', stack=(n+1))
    if notnull(ra1_in) then begin
      ra1_min[nsum] = ra1_min[nsum] < ra1_in
      ra1_max[nsum] = ra1_max[nsum] > ra1_in
      ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
      ra2_min[nsum] = ra2_min[nsum] < ra2_in
      ra2_max[nsum] = ra2_max[nsum] > ra2_in
    endif else begin
      ra1_min[nsum] = -999.9
      ra1_max[nsum] = 999.9
      ra2_min[nsum] = -999.9
      ra2_max[nsum] = 999.9
    endelse
    ra1_in = getextrai('ramax', stack=(n+1))
    if notnull(ra1_in) then begin
      ra1_min[nsum] = ra1_min[nsum] < ra1_in
      ra1_max[nsum] = ra1_max[nsum] > ra1_in
      ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
      ra2_min[nsum] = ra2_min[nsum] < ra2_in
      ra2_max[nsum] = ra2_max[nsum] > ra2_in
    endif else begin
      ra1_min[nsum] = -999.9
      ra1_max[nsum] = 999.9
      ra2_min[nsum] = -999.9
      ra2_max[nsum] = 999.9
    endelse
    dec_in = getextrai('decmin', stack=(n+1))
    dec_min[nsum] = (notnull(dec_in)) ? (dec_min[nsum] < dec_in) : -999.9
    dec_in = getextrai('decmax', stack=(n+1))
    dec_max[nsum] = (notnull(dec_in)) ? (dec_max[nsum] > dec_in) : 999.9

    ; Construct weighted means of some parameters

    if elevOK[nsum] then newelev[nsum] = newelev[nsum] + weight*elev
    if rttOK[nsum] then newrtt[nsum] = newrtt[nsum] + weight*rtt
    if txpowerOK[nsum] then newtxpower[nsum] = newtxpower[nsum] + weight*txpower
    if tsysOK[nsum] then newtsys[nsum] = newtsys[nsum] + weight*tsys
    if gainOK[nsum] then newgain[nsum] = newgain[nsum] + weight*gain
    if phaseOK[nsum] then begin
      phase1_mean[nsum] = phase1_mean[nsum] + weight*phase1
      phase2_mean[nsum] = phase2_mean[nsum] + weight*phase2
    endif
    if azimOK[nsum] then begin
      az1_mean[nsum] = az1_mean[nsum] + weight*az1
      az2_mean[nsum] = az2_mean[nsum] + weight*az2
    endif
    if timezoneOK[nsum] then begin
      tempzone = getextrai('timezone',stack=(n+1))
      timezoneOK[nsum] = (tempzone eq timezone[nsum])
    endif
    if jdmeanOK[nsum] then begin
      jd_in = getextrai('jdmean', stack=(n+1))
      if notnull(jd_in) then begin
        jdmean[nsum] = jdmean[nsum] + weight*jd_in
      endif else begin
        jdmeanOK[nsum] = 0
      endelse
    endif
    if dist_meanOK[nsum] then begin
      dist_in = getextrai('distmean', stack=(n+1))
      if notnull(dist_in) then begin
        dist_mean[nsum] = dist_mean[nsum] + weight*dist_in
      endif else begin
        dist_meanOK[nsum] = 0
      endelse
    endif
    if ra_meanOK[nsum] then begin
      ra1_in = getextrai('ramean', stack=(n+1))
      if notnull(ra1_in) then begin
        ra1_mean[nsum] = ra1_mean[nsum] + weight*ra1_in
        ra2_in = ra1_in - 24*floor(ra1_in/24.0 - 0.5)
        ra2_mean[nsum] = ra2_mean[nsum] + weight[0]*ra2_in
      endif else begin
        ra_meanOK[nsum] = 0
      endelse
    endif
    if dec_meanOK[nsum] then begin
      dec_in = getextrai('decmean', stack=(n+1))
      if notnull(dec_in) then begin
        dec_mean[nsum] = dec_mean[nsum] + weight*dec_in
      endif else begin
        dec_meanOK[nsum] = 0
      endelse
    endif

    ; Go through the array of extra tags and find the ones which are present and have
    ; the same format and value for all images contributing to this sum: those are the
    ; ones we'll keep.  Similarly, keep a comment line only if every image has the same
    ; line.  Finally, if formats and values match for a given extra tag but trailing
    ; comments don't match, keep the extra tag but erase the trailing comment.
    ;
    ; Keep time tags in order to keep them grouped together within the extra tags:
    ; we'll change the values and comments later

    for k=newnextra[nsum]-1L,0,-1L do begin

      newextraformat = strlowcase(newextratags[nsum,k].format)
      newextraname = strlowcase(newextratags[nsum,k].name)
      newextravalue = strlowcase(newextratags[nsum,k].value)
      newextracomment = strlowcase(newextratags[nsum,k].comment)
      if newextraname eq 'year' or newextraname eq 'month'  or newextraname eq 'day'    $
                                or newextraname eq 'hour'   or newextraname eq 'minute' $
                                or newextraname eq 'second' then begin
        keepextra = 1
        keepcomment = 0
      endif else begin
        keepextra = 0
        keepcomment = 0
        extranum = where(extratags.name eq newextraname, extracount)
        if extracount gt 0 then begin
          for j=0L,extracount-1 do begin
            extraformat = strlowcase(extratags[extranum[j]].format)
            extraname = newextraname
            extravalue = strlowcase(extratags[extranum[j]].value)
            extracomment = strlowcase(extratags[extranum[j]].comment)
            match = ((notnull(extraname)) and $
                     (extraformat eq newextraformat) and (extravalue eq newextravalue)) or $
                    ((isnull(extraname)) and (extracomment eq newextracomment)) or $
                    ((extraname eq 'eph_row'))
            keepextra = keepextra or match
            keepcomment = keepcomment or (match and (extracomment eq newextracomment))
          endfor
        endif
      endelse

      if not keepextra then begin

        ; Don't create an error by trying to remove the last array element

        if newnextra[nsum] gt 1 then begin
          if k eq 0 then begin
            newextratags[nsum,0:newnextra[nsum]-2] = newextratags[nsum,1:newnextra[nsum]-1]
          endif else if k lt newnextra[nsum]-1 then begin
            newextratags[nsum,k:newnextra[nsum]-2] = newextratags[nsum,k+1:newnextra[nsum]-1]
          endif
        endif

        newnextra[nsum] = newnextra[nsum] - 1

      endif else if not keepcomment then begin
        newextratags[nsum,k].comment = ''
      endif

    endfor

    ; Increment the number of images contributing to this sum

    n_in[nsum] = n_in[nsum] + 1L

  endif  ;  count gt 0 and sumpol[ch-1] ?

; Look at the next image in the image stack

endfor

; All stacki images have been inspected

for nsum=0L,nsums-1 do begin

  if n_in[nsum] ge 2 then begin

    ; Normalize the weighted mean images and get the rms noise (sdev)
    ;
    ; Note that dividing by the sum of all weights, then normalizing
    ; by dividing by sdev = 1/sqrt(sum of all weights), is equivalent
    ; to multiplying by sdev

    if keyword_set(noweight) then begin
      if keyword_set(calibrated) then begin
        sdev = sqrt(sumvariance[nsum]) / n_in[nsum]
        (*wsum[nsum]).image = (*wsum[nsum]).image / sqrt(sumvariance[nsum])
      endif else begin
        (*wsum[nsum]).image = (*wsum[nsum]).image / sqrt(n_in[nsum])
      endelse
    endif else begin
      sdev = 1.0/sqrt(sumweights[nsum])
      (*wsum[nsum]).image = sdev * (*wsum[nsum]).image
    endelse

    ; Complete the weighted means for various tags

    if elevOK[nsum] then newelev[nsum] = newelev[nsum] / sumweights[nsum]
    if rttOK[nsum] then newrtt[nsum] = newrtt[nsum] / sumweights[nsum]
    if txpowerOK[nsum] then newtxpower[nsum] = newtxpower[nsum] / sumweights[nsum]
    if tsysOK[nsum] then newtsys[nsum] = newtsys[nsum] / sumweights[nsum]
    if gainOK[nsum] then newgain[nsum] = newgain[nsum] / sumweights[nsum]
    if jdmeanOK[nsum] then jdmean[nsum] = jdmean[nsum] / sumweights[nsum]
    if dist_meanOK[nsum] then dist_mean[nsum] = dist_mean[nsum] / sumweights[nsum]
    if ra_meanOK[nsum] then begin
      ra1_mean[nsum] = ra1_mean[nsum] / sumweights[nsum]
      ra2_mean[nsum] = ra2_mean[nsum] / sumweights[nsum]
    endif
    if dec_meanOK[nsum] then dec_mean[nsum] = dec_mean[nsum] / sumweights[nsum]

    ; Deal with rotation phase and azimuth:
    ;
    ; If the images contributing to this sum have phases which are more
    ; "concentrated" when considered over the range [180,540) rather than [0,360),
    ; then compute the mean phase using the former range (and then take mod 360).
    ; For example, a bunch of phases between 330-360 and 0-40 should yield a
    ; mean phase close to 360/0, not 180.  Ditto azimuth.

    if phaseOK[nsum] then begin
      phase1_range = phase1_max[nsum] - phase1_min[nsum]
      phase2_range = phase2_max[nsum] - phase2_min[nsum]
      if phase1_range gt phase2_range then begin
        newphase = (phase2_mean[nsum]/sumweights[nsum]) mod 360.0
      endif else begin
        newphase = phase1_mean[nsum]/sumweights[nsum]
      endelse
    endif
    if azimOK[nsum] then begin
      az1_range = az1_max[nsum] - az1_min[nsum]
      az2_range = az2_max[nsum] - az2_min[nsum]
      if az1_range gt az2_range then begin
        newazim = (az2_mean[nsum]/sumweights[nsum]) mod 360.0
      endif else begin
        newazim = az1_mean[nsum]/sumweights[nsum]
      endelse
    endif

    ; If all images contributing to this sum had mean Julian date extra tags,
    ; use the new mean Julian date to set some time/date tags

    if jdmeanOK[nsum] and timezoneOK[nsum] then begin
      jd_midnight = floor(jdmean[nsum] - 0.5) + 0.5D
      rctime = round( 86400*((jdmean[nsum] - 0.5D) mod 1) )  ;  nearest sec
      jduse = jd_midnight + rctime/86400.0
      caldat_roundsec,jdmean[nsum],mon,day,year,hr,min,sec
      calmean[nsum] = string(year,format='(i4)') + ' '     $
                      + monthnames[mon-1] + ' '            $
                      + string(day,format='(i2.2)') + ' '  $
                      + string(hr,format='(i2.2)') + ':'   $
                      + string(min,format='(i2.2)') + ':'  $
                      + string(sec,format='(i2.2)') + ' '  $
                      + timezone[nsum]
    endif else begin
      year = 1L
      mon = 1L  ;  so showstacki can call it 'Jan'
      day = 1L
      hr = 0L
      min = 0L
      sec = 0L
      if jdmeanOK[nsum] then begin
        print,'WARNING in sumim: No time information output for sum #', $
              nsum+1,format='(a,i0)'
        print,"                     due to discrepant 'timezone' extra tags"
      endif
    endelse

    ; Make sure there's at least one extra tag (so that the extratags
    ; array is defined): Add a dummy tag.

    if newnextra[nsum] eq 0 then begin
      newnextra[nsum] = 1L
      newextratags[nsum,0] = extratag
      newextratags[nsum,0].comment = 'This is a dummy comment'
    endif

    ; Put it all together and add it to the image stack

    ; Figure out what the sums' group numbers will be

    maxgroup = max(allgroups)
    newstart = 100*( 1L + maxgroup/100 ) + 1L
    newgroupnumber = lindgen(nsums) + newstart

    stackiStruc = {group:newgroupnumber[nsum], $
                   image:(*wsum[nsum]).image, $
                   extratags:newextratags[nsum,0:newnextra[nsum]-1], $ 
                   width:newwidth[nsum], $
                   height:newheight[nsum], $
                   nextra:newnextra[nsum], $
                   tname:newtname[nsum], $
                   pol:newpol[nsum]}
    stacki = [stacki, ptr_new(stackiStruc)]
    nstacki = nstacki + 1L
    print,'Sum #',nsum+1,' (',n_in[nsum],' ',newpol[nsum], $
          ' images) is stacki image #',nstacki,' (group #',newgroupnumber[nsum],')', $
          format='(a,i0,a,i0,3a,i0,a,i0,a)'

    ; Add various extra tags which have been changed

    setextrai,'i','year',string(year, format='(i0)'), $
              comment='year of mean RX date/time (UTC)',stack=nstacki
    setextrai,'i','month',string(mon, format='(i0)'), $
              comment='month of mean RX date/time (UTC)',stack=nstacki
    setextrai,'i','day',string(day, format='(i0)'), $
              comment='day of mean RX date/time (UTC)',stack=nstacki
    setextrai,'i','hour',string(hr, format='(i0)'), $
              comment='hour of mean RX date/time (UTC)',stack=nstacki
    setextrai,'i','minute',string(min, format='(i0)'), $
              comment='minute of mean RX date/time',stack=nstacki
    setextrai,'i','second',string(sec, format='(i0)'), $
              comment='second of mean RX date/time',stack=nstacki
    if year ge 1600 then begin
      doy = date2doy(mon, day, year)
      setextrai,'i','iday',string(doy, format='(i0)'), $
                comment='Day of Year for mean RX (UT January 1 = 1)',stack=nstacki
    endif else begin
      deleteextrai,'iday',/silent,stack=nstacki
    endelse
    if nlooksOK[nsum] then $
         setextrai,'i','nlooks',string(newnlooks[nsum], format='(i0)'), $
                   stack=nstacki
    if tauOK[nsum] then $
         setextrai,'f','int_time',newtau[nsum], $
                   comment='Exposure time included in image (s)',stack=nstacki
    setextrai,'f','eph_row',neweph_row[nsum], $
              comment='row in which ephemeris range lies, 0-based',stack=nstacki
    if elevOK[nsum] then $
         setextrai,'f','elevation',string(newelev[nsum], format='(f8.5)'), $
                   comment='mean elevation in deg',stack=nstacki
    if azimOK[nsum] then $
         setextrai,'f','azimuth',string(newazim, format='(f8.4)'), $
                   comment='mean azimuth in deg',stack=nstacki
    if rttOK[nsum] then $
         setextrai,'f','rtt',string(newrtt[nsum], format='(f10.5)'), $
                   comment='mean round-trip time in sec',stack=nstacki
    if txpowerOK[nsum] then begin
      deleteextrai,'tx_power',/silent,stack=nstacki
      setextrai,'i','tx_power',string(round(newtxpower[nsum]), format='(i0)'), $
                comment='mean TX power in kW',stack=nstacki
    endif
    if tsysOK[nsum] then $
         setextrai,'f','tsys',string(newtsys[nsum], format='(f8.4)'), $
                   comment='mean system temp in K',stack=nstacki
    if gainOK[nsum] then $
         setextrai,'f','gain_rxtx',string(newgain[nsum], format='(e11.5)'), $
                   comment='mean rx_gain * tx_gain',stack=nstacki
    if phaseOK[nsum] then $
         setextrai,'f','phase',string(newphase, format='(f8.4)'), $
                   comment='mean rotation phase in deg',stack=nstacki
    if not (keyword_set(noweight) and not keyword_set(calibrated)) then $
         setextrai,'f','sdev',string(sdev, format='(e10.4)'), $
                   comment='rms pixel noise in km^2',stack=nstacki
    if delbinOK[nsum] then begin
      setextrai,'i','delbin1',string(newdelbin1[nsum], format='(i0)'), $
                comment='min delay bin containing signal, 0-based',stack=nstacki
      setextrai,'i','delbin2',string(newdelbin2[nsum], format='(i0)'), $
                comment='max delay bin containing signal, 0-based',stack=nstacki
    endif
    if dopbinOK[nsum] then begin
      setextrai,'i','dopbin1',string(newdopbin1[nsum], format='(i0)'), $
                comment='min Doppler bin containing signal, 0-based',stack=nstacki
      setextrai,'i','dopbin2',string(newdopbin2[nsum], format='(i0)'), $
                comment='max Doppler bin containing signal, 0-based',stack=nstacki
    endif

    ; If all images contributing to this sum had extra tags for Julian dates,
    ; get rid of any old Julian date extra tags and create new ones
    ;
    ; Ditto for distance, RA, and dec

    if timezoneOK[nsum] then begin
      if jdstart[nsum] gt 0 then $
         setextrai,'d','jdstart',string(jdstart[nsum], format='(d13.5)'),stack=nstacki
      if jdmeanOK[nsum] then begin
         setextrai,'d','jdmean',string(jdmean[nsum], format='(d13.5)'),stack=nstacki
         setextrai,'s','calmean',calmean[nsum],stack=nstacki
      endif
      if jdend[nsum] lt 10000000L then $
         setextrai,'d','jdend',string(jdend[nsum], format='(d13.5)'),stack=nstacki
    endif
    if dist_min[nsum] gt 0 then $
       setextrai,'f','distmin',dist_min[nsum],stack=nstacki
    if dist_meanOK[nsum] then $
       setextrai,'f','distmean',dist_mean[nsum],stack=nstacki
    if dist_max[nsum] lt 100000.0 then $
       setextrai,'f','distmax',dist_max[nsum],stack=nstacki
    if ra1_min[nsum] ge  0 and ra1_max[nsum] lt 24 and $
       ra2_min[nsum] ge 12 and ra2_max[nsum] lt 36 and ra_meanOK[nsum] then begin
      ra1_range = ra1_max[nsum] - ra1_min[nsum]
      ra2_range = ra2_max[nsum] - ra2_min[nsum]
      if ra2_range lt ra1_range then begin
        setextrai,'f','ramin',(ra2_min[nsum] mod 24),stack=nstacki
        setextrai,'f','ramean',(ra2_mean[nsum] mod 24),stack=nstacki
        setextrai,'f','ramax',(ra2_max[nsum] mod 24),stack=nstacki
      endif else begin
        setextrai,'f','ramin',ra1_min[nsum],stack=nstacki
        setextrai,'f','ramean',ra1_mean[nsum],stack=nstacki
        setextrai,'f','ramax',ra1_max[nsum],stack=nstacki
      endelse
    endif
    if dec_min[nsum] ge -90.0 then $
       setextrai,'f','decmin',dec_min[nsum],stack=nstacki
    if dec_meanOK[nsum] then $
       setextrai,'f','decmean',dec_mean[nsum],stack=nstacki
    if dec_max[nsum] le 90.0 then $
       setextrai,'f','decmax',dec_max[nsum],stack=nstacki

    ; If the 'date-obs' (FITS) tag was present in the input images, there may be some
    ; standard comments which elaborate on it referring to rx start: delete them

    extranum = where((*stacki[nstacki-1]).extratags.comment eq $
                     '# The RX start time is a reasonable approximation to the experiment', count)
    if count gt 0 then for k=count-1L,0,-1 do deleteextrai,extranum[k]+1,stack=nstacki,/silent
    extranum = where((*stacki[nstacki-1]).extratags.comment eq $
                     '# mid-time. A more accurate mid-time is DATE-OBS + (EXPTIME-RTT)/2', count)
    if count gt 0 then for k=count-1L,0,-1 do deleteextrai,extranum[k]+1,stack=nstacki,/silent
    extranum = where((*stacki[nstacki-1]).extratags.comment eq $
                     '# which is typically about 3s earlier', count)
    if count gt 0 then for k=count-1L,0,-1 do deleteextrai,extranum[k]+1,stack=nstacki,/silent

  endif else begin

    ; Sums with zero or one image contributing: Just print a message

    print,'Sum #',nsum+1,' (',n_in[nsum],' ',newpol[nsum],' images) will not be output', $
          format='(a,i0,a,i0,3a)'

  endelse

endfor

; Clean up pointers

ptr_free,wsum

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showpeak,chan=chan,full=full,stack=st,help=help

; Display the maximum signal strength for the loaded pair;
; search only within the defined signal bins, unless /full is set.
;
; /stack does this instead for all pairs in the stack
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

if n_params() gt 0 or keyword_set(help) then begin
  print,' '
  print,'showpeak[,chan=1 or 2][,/full][,/stack][,/help]'
  print,' '
  print,'Display the peak signal strength for the loaded pair'
  print,'    -- look only within the defined signal bins (jsnr1-jsnr2),'
  print,'       unless /full is set, in which case consider all bins'
  print,' '
  print,'/stack displays maximum signal strength for all pairs'
  print,'     in the stack rather than the loaded pair'
  print,' '
  return
endif

cchanstrings = chanstrings
for i = 0, n_elements(chanstrings)-1 do cchanstrings[i] = chanstrings[i] + ': ' 

if not keyword_set(st) then begin

  ; Display peak signal strength for the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in showpeak: No pair is loaded'
    return
  endif
  nuse = 1L
  npol = n_elements((*loaded).tags)
endif else begin

  ; Display peak signal strength for every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in showpeak: The pair stack is empty'
    return
  endif
  nuse = nstack
  npol = n_elements((*stack[nuse]).tags)
endelse

; Decide which channel to use

showpol = intrr(npol)
if n_elements(chan) eq 0 then begin
  showpol = showpol + 1
endif else if chan ge 1 && chan le npol then begin
  showpol[chan-1] = 1
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
  return
endelse

print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this pair

  useStruc = (keyword_set(st)) ? *stack[n] : *loaded
  spec = useStruc.spec
  ndata = useStruc.ndata
  tags = useStruc.tags
  jsnr1 = tags.jsnr1
  jsnr2 = tags.jsnr2
  if keyword_set(st) then begin
    posfr = tags[0].posfr
    xjcen = tags[0].xjcen
    dfreq = tags[0].dfreq
    freq = posfr*dfreq*(findgen(ndata) - xjcen)
  endif else begin
    freq = useStruc.freq
  endelse

  ; Set the frequency range over which we want the peak signal

  if keyword_set(full) then begin
    firstbin = intarr(npol)
    lastbin = replicate(ndata-1, npol)
  endif else begin
    firstbin = jsnr1
    lastbin = jsnr2
  endelse

  ; Display the peak signal and the frequency at which it occurs

  if keyword_set(st) then begin
    printstring = 'Pair #' + string(n+1,format='(i0)') + ':  '
  endif else begin
    printstring = ''
  endelse
  for ch=1,npol do begin
    if showpol[ch-1] then begin
      peaksignal = max(spec[ch-1,firstbin[ch-1]:lastbin[ch-1]], maxbin)
      peakfreq = freq[maxbin+firstbin[ch-1]]
      printstring = printstring + cchanstrings[ch-1]  $
                    + string(peaksignal,format='(f8.2)') + ' at '    $
                    + string(peakfreq,format='(f8.2)') + ' Hz     '
    endif
  endfor
  print,printstring
  print,' '
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showpeak1,full=full,stack=st,help=help

; Display the maximum signal strength for the loaded single-channel spectrum;
; search only within the defined signal bins, unless /full is set.
;
; /stack does this instead for all spectra in the single-channel stack
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() gt 0 or keyword_set(help) then begin
  print,' '
  print,'showpeak1[,/full][,/stack][,/help]'
  print,'  Display the peak signal strength for the loaded single-channel spectrum'
  print,'  -- look only within the defined signal bins (jsnr1-jsnr2),'
  print,'     unless /full is set, in which case consider all bins'
  print,' '
  print,'/stack displays maximum signal strength for all spectra'
  print,'     in the single-channel stack rather than the loaded spectrum'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Display peak signal strength for the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in showpeak1: No single-channel spectrum is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Display peak signal strength for every spectrum in the single-chanel stack

  if nstack1 eq 0 then begin
    print,'ERROR in showpeak1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
endelse

print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack1[n] : *loaded1
  spec = useStruc.spec
  ndata = useStruc.ndata
  tags1 = useStruc.tags
  jsnr1 = tags1.jsnr1
  jsnr2 = tags1.jsnr2
  if keyword_set(st) then begin
    posfr = tags1.posfr
    xjcen = tags1.xjcen
    dfreq = tags1.dfreq
    freq = posfr*dfreq*(findgen(ndata) - xjcen)
  endif else begin
    freq = useStruc.freq
  endelse

  ; Set the frequency range over which we want the peak signal

  if keyword_set(full) then begin
    firstbin = 0
    lastbin = ndata - 1
  endif else begin
    firstbin = jsnr1
    lastbin = jsnr2
  endelse

  ; Display the peak signal and the frequency at which it occurs

  if keyword_set(st) then begin
    printstring = 'Spectrum #' + string(n+1,format='(i0)') + ':  '
  endif else begin
    printstring = ''
  endelse
  peaksignal = max(spec[firstbin:lastbin], maxbin)
  peakfreq = freq[maxbin+firstbin]
  printstring = printstring + string(peaksignal,format='(f8.2)') $
                + ' at ' + string(peakfreq,format='(f8.2)') + ' Hz'
  print,printstring
  print,' '
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showpeaki,full=full,stack=st,help=help

; Display the maximum signal strength for the loaded image;
; search only within the defined signal bins, unless /full is set.
;
; /stack does this instead for all images in the image stack

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() gt 0 or keyword_set(help) then begin
  print,' '
  print,'showpeaki[,/full][,/stack][,/help]'
  print,'  Display the peak signal strength for the loaded image'
  print,'  -- look only within the defined signal bins (jsnr1-jsnr2),'
  print,'     unless /full is set, in which case consider all bins'
  print,' '
  print,'/stack displays maximum signal strength for all images'
  print,'     in the image stack rather than the loaded image'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Display peak signal strength for the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,'ERROR in showpeaki: No image is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Display peak signal strength for every image in the stack

  if nstacki eq 0 then begin
    print,'ERROR in showpeaki: The image stack is empty'
    return
  endif
  nuse = nstacki
endelse

print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this image

  if keyword_set(st) then begin
    image = (*stacki[n]).image
    width = (*stacki[n]).width
    height = (*stacki[n]).height
    dfreq = getextrai('fres',stack=(n+1))
    eph_col = getextrai('eph_col',stack=(n+1))
    eph_row = getextrai('eph_row',stack=(n+1))
    delayunit = getextrai('delayunit',stack=(n+1))
    spb = getextrai('samples_per_baud',stack=(n+1))
    rows_per_baud = getextrai('rows_per_baud',stack=(n+1))
    baud = getextrai('baudlen',stack=(n+1))
    dopbin1 = getextrai('dopbin1',stack=(n+1))
    dopbin2 = getextrai('dopbin2',stack=(n+1))
    delbin1 = getextrai('delbin1',stack=(n+1))
    delbin2 = getextrai('delbin2',stack=(n+1))
  endif else begin
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
    dopbin1 = getextrai('dopbin1')
    dopbin2 = getextrai('dopbin2')
    delbin1 = getextrai('delbin1')
    delbin2 = getextrai('delbin2')
  endelse
  if isnull(delayunit) and notnull(baud) then begin
    if isnull(spb) then spb = 1L
    if isnull(rows_per_baud) then rows_per_baud = spb
    delayunit = baud/rows_per_baud
  endif
  if notnull(eph_col) and notnull(eph_row) and notnull(dfreq) and notnull(delayunit) then begin
    f = dfreq*(findgen(width) - eph_col)       ; Doppler (Hz)
    d = delayunit*(findgen(height) - eph_row)  ; delay (usec)
    good_deldop = 1
  endif else begin
    print,"Image #",n+1," needs the 'eph_col' and 'eph_row' and 'fres' tags,", $
          format='(a,i0,a)'
    print,"           plus either the 'delayunit' or 'baudlen' tag"
    good_deldop = 0
  endelse

  ; Set the delay range over which we want the peak signal

  if keyword_set(full) then begin
    firstrow = 0
    lastrow = height - 1
  endif else begin
    if isnull(delbin1) then begin
      delbin1 = 0L
      print,'delbin1 missing for image #',n+1,', set to ',delbin1, $
            format='(a,i0,a,i0)'
    endif else if delbin1 lt 0 then begin
      delbin1 = 0L
      print,'delbin1 < 0 for image #',n+1,', reset to ',delbin1, $
            format='(a,i0,a,i0)'
    endif
    if isnull(delbin2) then begin
      delbin2 = height - 1
      print,'delbin2 missing for image #',n+1,', set to ',delbin2, $
            format='(a,i0,a,i0)'
    endif else if delbin2 ge height then begin
      delbin2 = height - 1
      print,'delbin2 >= height for image #',n+1,', reset to ',delbin2, $
            format='(a,i0,a,i0)'
    endif
    firstrow = delbin1
    lastrow = delbin2
  endelse

  ; Set the frequency range over which we want the peak signal

  if keyword_set(full) then begin
    firstcol = 0
    lastcol = width - 1
  endif else begin
    if isnull(dopbin1) then begin
      dopbin1 = 0L
      print,'dopbin1 missing for image #',n+1,', set to ',dopbin1, $
            format='(a,i0,a,i0)'
    endif else if dopbin1 lt 0 then begin
      dopbin1 = 0L
      print,'dopbin1 < 0 for image #',n+1,', reset to ',dopbin1, $
            format='(a,i0,a,i0)'
    endif
    if isnull(dopbin2) then begin
      dopbin2 = width - 1
      print,'dopbin2 missing for image #',n+1,', set to ',dopbin2, $
            format='(a,i0,a,i0)'
    endif else if dopbin2 ge width then begin
      dopbin2 = width - 1
      print,'dopbin2 >= width for image #',n+1,', reset to ',dopbin2, $
            format='(a,i0,a,i0)'
    endif
    firstcol = dopbin1
    lastcol = dopbin2
  endelse

  ; Display the peak signal and the delay and frequency at which it occurs

  if delbin2 lt delbin1 then begin
    print,'Image #',n+1,' has delbin2 < delbin1: no computation possible', $
          format='(a,i0,a)'
  endif else if dopbin2 lt dopbin1 then begin
    print,'Image #',n+1,' has dopbin2 < dopbin1: no computation possible', $
          format='(a,i0,a)'
  endif else begin

    peaksignal = max(image[firstcol:lastcol,firstrow:lastrow], maxelement)
    width_use = lastcol - firstcol + 1
    maxcol = firstcol + (maxelement mod width_use)
    maxrow = firstrow + maxelement/width_use
    if good_deldop then begin
      print,'Image #',n+1,': ',peaksignal,' at ',d[maxrow],' usec, ', $
            f[maxcol],' Hz  (row ',maxrow,', column ',maxcol,')', $
            format='(a,i0,a,f8.2,a,f8.2,a,f8.3,a,i0,a,i0,a)'
    endif else begin
      print,'Image #',n+1,': ',peaksignal,' at row ',maxrow,', column ', $
            maxcol,format='(a,i0,a,f8.2,a,i0,a,i0)'
    endelse
    print,' '

  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showsnr,efbmin,efbmax,efbincr,n=n,help=help

; Display the signal-to-noise ratio (SNR) of the OC channel of the loaded pair.
; SNR is computed by smoothing the spectrum to effective resolutions from
; efbmin to efbmax in increments of efbincr (5 Hz if not specified), finding the
; peak signal (within the defined signal limits) for each resolution, and using the
; highest such peak.
;
; 2002 Nov 28: Allow user to specify power-law exponent n, which determines the
;              shape of the smoothing filter (see comments to procedure smoothf).
;              Note that n needn't be an integer.  Default value = 2.
 
common loadedBlock,loadedi,loaded1,loaded

if n_elements(n) eq 0 then n = 2.0

if n_params() lt 2 or n_params() gt 3 or keyword_set(help) then begin
  print,' '
  print,'showsnr,efbmin,efbmax[,efbincr][,n=n][,/help]'
  print,' '
  print,'Display the OC signal-to-noise ratio (SNR) of the loaded pair,'
  print,'     the peak strength of the optimally filtered OC spectrum'
  print,' '
  print,'     Smooth to trial effective resolutions between efbmin and efbmax,'
  print,'          at increments of efbincr (or 5 Hz if not specified)'
  print,' '
  print,'     n determines the shape of the smoothing filter; type "smoothf,/help"'
  print,"          for details.  Note that n needn't be an integer."
  print,'          (default = 2)'
  print,' '
  return
endif

if (*loaded).ndata le 2 then begin
  print,"ERROR in showsnr: No pair is loaded, so SNR can't be determined"
  return
endif else if efbmax lt efbmin then begin
  print,"ERROR in showsnr: Can't have efbmax < efbmin"
  return
endif else if n lt 0 then begin
  print,"ERROR in showsnr: Power-law exponent n can't be negative"
  return
endif

dfreq = (*loaded).tags[0].dfreq
if efbmax lt dfreq then begin
  print,"ERROR in showsnr: Can't have efbmax < raw resolution (",dfreq," Hz)", $
        format='(a,f7.3,a)'
  return
endif

; Set the effective resolution increment if not specified

if n_params() eq 2 then efbincr = 5

; Make sure that both min and max effective resolutions are at least double
; the raw resolution, so that there's no trouble with smoothing

if efbmin ge 2*dfreq then begin
  efbmin_use = 1.0*efbmin
  efbmax_use = 1.0*efbmax
endif else begin
  efbmin_use = 1.0*efbincr*ceil(2*dfreq/efbincr)
  print,'WARNING in showsnr: Reset minimum resolution to ',efbmin_use, $
        ' Hz, which is >= 2 x raw resolution',format='(a,f7.3,a)'
  if efbmax ge efbmin_use then begin
    efbmax_use = 1.0*efbmax
  endif else begin
    efbmax_use = 1.0*efbmin_use
    print,'WARNING in showsnr: Reset maximum resolution to ',efbmax_use, $
          ' Hz',format='(a,f7.3,a)'
  endelse
endelse

; Store the loaded pair so we can "refresh" it after each smoothing operation

storeloaded = ptr_new(*loaded)

; Compute the peak signal for each smoothed OC spectrum

lim1 = (*loaded).tags[0].jsnr1
lim2 = (*loaded).tags[0].jsnr2
SNR = 0.0
efbSNR = 0.0
for efb=efbmin_use,efbmax_use,efbincr do begin
  smoothf,efb,n=n,chan=1
  peakOCsignal = max((*loaded).spec[0,lim1:lim2])
  if peakOCsignal gt SNR then begin
    SNR = peakOCsignal
    efbSNR = efb
  endif
  *loaded = *storeloaded
endfor

; Display the results

print,' '
print,'OC SNR = ',SNR,' (at effective resolution ',efbSNR,' Hz)', $
      format='(a,f8.2,a,f8.2,a)'
print,' '

; Free up storage space

ptr_free,storeloaded

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showsnr1,efbmin,efbmax,efbincr,n=n,help=help

; Display the signal-to-noise ratio (SNR) of the loaded single-channel spectrum.
; SNR is computed by smoothing the spectrum to effective resolutions from
; efbmin to efbmax in increments of efbincr (5 Hz if not specified), finding the
; peak signal (within the defined signal limits) for each resolution, and using the
; highest such peak.
;
; 2002 Nov 28: Allow user to specify power-law exponent n, which determines the
;              shape of the smoothing filter (see comments to procedure smoothf).
;              Note that n needn't be an integer.  Default value = 2.
 
common loadedBlock,loadedi,loaded1,loaded

if n_elements(n) eq 0 then n = 2.0

if n_params() lt 2 or n_params() gt 3 or keyword_set(help) then begin
  print,' '
  print,'showsnr1,efbmin,efbmax[,efbincr][,n=n][,/help]'
  print,' '
  print,'Display the signal-to-noise ratio (SNR) of the loaded single-channel spectrum,'
  print,'     the peak strength of the optimally filtered spectrum'
  print,' '
  print,'     Smooth to trial effective resolutions between efbmin and efbmax,'
  print,'          at increments of efbincr (or 5 Hz if not specified)'
  print,' '
  print,'     n determines the shape of the smoothing filter; type "smoothf1,/help"'
  print,"          for details.  Note that n needn't be an integer."
  print,'          (default = 2)'
  print,' '
  return
endif

if (*loaded1).ndata le 2 then begin
  print,"ERROR in showsnr1: No single-channel spectrum is loaded, so SNR can't be determined"
  return
endif else if efbmax lt efbmin then begin
  print,"ERROR in showsnr1: Can't have efbmax < efbmin"
  return
endif else if n lt 0 then begin
  print,"ERROR in showsnr1: Power-law exponent n can't be negative"
  return
endif

dfreq = (*loaded1).tags.dfreq
if efbmax lt dfreq then begin
  print,"ERROR in showsnr1: Can't have efbmax < raw resolution (",dfreq," Hz)", $
        format='(a,f7.3,a)'
  return
endif

; Set the effective resolution increment if not specified

if n_params() eq 2 then efbincr = 5

; Make sure that both min and max effective resolutions are at least double
; the raw resolution, so that there's no trouble with smoothing

if efbmin ge 2*dfreq then begin
  efbmin_use = 1.0*efbmin
  efbmax_use = 1.0*efbmax
endif else begin
  efbmin_use = 1.0*efbincr*ceil(2*dfreq/efbincr)
  print,'WARNING in showsnr1: Reset minimum resolution to ',efbmin_use, $
        ' Hz, which is >= 2 x raw resolution',format='(a,f7.3,a)'
  if efbmax ge efbmin_use then begin
    efbmax_use = 1.0*efbmax
  endif else begin
    efbmax_use = 1.0*efbmin_use
    print,'WARNING in showsnr1: Reset maximum resolution to ',efbmax_use, $
          ' Hz',format='(a,f7.3,a)'
  endelse
endelse

; Store the loaded single-channel spectrum so we can "refresh" it
; after each smoothing operation

storeloaded1 = ptr_new(*loaded1)

; Compute the peak signal for each smoothed spectrum

lim1 = (*loaded1).tags.jsnr1
lim2 = (*loaded1).tags.jsnr2
SNR = 0.0
efbSNR = 0.0
for efb=efbmin_use,efbmax_use,efbincr do begin
  smoothf1,efb,n=n
  peakSignal = max((*loaded1).spec[lim1:lim2])
  if peakSignal gt SNR then begin
    SNR = peakSignal
    efbSNR = efb
  endif
  *loaded1 = *storeloaded1
endfor

; Display the results

print,' '
print,'SNR = ',SNR,' (at effective resolution ',efbSNR,' Hz)', $
      format='(a,f8.2,a,f8.2,a)'
print,' '

; Free up storage space

ptr_free,storeloaded1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showlim,chan=chan,stack=n,help=help

; Displays the signal limits for the loaded pair,
; or else for stack pair n if the stack keyword is set.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

if keyword_set(help) or n_params() gt 0 then begin
  print,'showlim[,chan=1 or 2][,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in showlim: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showlim: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in showlim: There are only ',nstack, $
          ' pairs in the stack', format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in showlim: No pair is loaded'
  return
endif

cchanstrings = chanstrings
for i=0, n_elements(chanstrings)-1 do cchanstrings[i] = chanstrings[i] + ': '

; Get the necessary elements of the structure
; whose parameters are to be displayed

dispStruc = (n_elements(n) gt 0) ? *stack[n-1] : *loaded
tags = dispStruc.tags
df = tags[0].dfreq
posfr = tags[0].posfr
bin0 = tags[0].xjcen
bin1 = tags.jsnr1
bin2 = tags.jsnr2
flim1 = posfr*df*(bin1 - 0.5 - bin0)
flim2 = posfr*df*(bin2 + 0.5 - bin0)
ndata = dispStruc.ndata
npol = n_elements(tags)

; Decide which channel to use

showpol = intarr(npol)
if n_elements(chan) eq 0 then begin
  showpol = showpol + 1
endif else if chan ge 1 && chan le npol then begin
  showpol[chan-1] = 1
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
  return
endelse
; Display the limits

printstring = ''
for ch=1,npol do begin
  if showpol[ch-1] then begin
    printstring = printstring + cchanstrings[ch-1] + 'freq. bins '     $
                  + string(bin1[ch-1],format='(i0)') + '-'            $
                  + string(bin2[ch-1],format='(i0)') + '  ('          $
                  + string(flim1[ch-1],format='(f8.2)') + ' Hz to '   $
                  + string(flim2[ch-1],format='(f8.2)') + ' Hz)     '
  endif
endfor
print,' '
print,printstring
if bin0 ne ndata/2 then $
    print,'    (Full range is bins 0-',ndata-1,', 0 Hz is bin ',bin0,')', $
          format='(a,i0,a,i0,a)'
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showlim1,stack=n,help=help

; Displays the signal limits for the loaded single-channel spectrum,
; or else for stack1 spectrum n if the stack keyword is set.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showlim1[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack1 eq 0 then begin
    print,'ERROR in showlim1: The single-channel stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showlim1: Must have n >= 1'
    return
  endif else if n gt nstack1 then begin
    print,'ERROR in showlim1: There are only ',nstack1, $
          ' spectra in the single-channel stack', format='(a,i0,a)'
    return
  endif
endif else if (*loaded1).ndata le 2 then begin
  print,'ERROR in showlim1: No single-channel spectrum is loaded'
  return
endif

; Get the necessary elements of the structure
; whose parameters are to be displayed

dispStruc = (n_elements(n) gt 0) ? *stack1[n-1] : *loaded1
tags1 = dispStruc.tags
df = tags1.dfreq
posfr = tags1.posfr
bin0 = tags1.xjcen
bin1 = tags1.jsnr1
bin2 = tags1.jsnr2
flim1 = posfr*df*(bin1 - 0.5 - bin0)
flim2 = posfr*df*(bin2 + 0.5 - bin0)
ndata = dispStruc.ndata

; Display the limits

print,' '
print,'Frequency bins ',bin1,'-',bin2,'  (',flim1,' Hz to ',  $
      flim2,' Hz)',format='(a,i0,a,i0,a,f8.2,a,f8.2,a)'
if bin0 ne ndata/2 then $
    print,'    (Full range is bins 0-',ndata-1,', 0 Hz is bin ',bin0,')', $
          format='(a,i0,a,i0,a)'
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showlimi,stack=n,help=help

; Displays the signal limits for the loaded image,
; or else for stacki image n if the stack keyword is set.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'showlimi[,stack=n][,/help])'
  return
endif else if n_elements(n) gt 0 then begin
  if nstacki eq 0 then begin
    print,'ERROR in showlimi: The image stack is empty'
    return
  endif else if n le 0 then begin
    print,'ERROR in showlimi: Must have n >= 1'
    return
  endif else if n gt nstacki then begin
    print,'ERROR in showlimi: There are only ',nstacki, $
          ' images in the image stack', format='(a,i0,a)'
    return
  endif
endif else if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in showlimi: No image is loaded'
  return
endif

; Get elements of the image whose parameters are to be displayed

if n_elements(n) gt 0 then begin
  width = (*stacki[n-1]).width
  height = (*stacki[n-1]).height
  dfreq = getextrai('fres',stack=n)
  eph_col = getextrai('eph_col',stack=n)
  eph_row = getextrai('eph_row',stack=n)
  delayunit = getextrai('delayunit',stack=n)
  spb = getextrai('samples_per_baud',stack=n)
  rows_per_baud = getextrai('rows_per_baud',stack=n)
  baud = getextrai('baudlen',stack=n)
  dopbin1 = getextrai('dopbin1',stack=n)
  dopbin2 = getextrai('dopbin2',stack=n)
  delbin1 = getextrai('delbin1',stack=n)
  delbin2 = getextrai('delbin2',stack=n)
endif else begin
  width = (*loadedi).width
  height = (*loadedi).height
  dfreq = getextrai('fres')
  eph_col = getextrai('eph_col')
  eph_row = getextrai('eph_row')
  delayunit = getextrai('delayunit')
  spb = getextrai('samples_per_baud')
  rows_per_baud = getextrai('rows_per_baud')
  baud = getextrai('baudlen')
  dopbin1 = getextrai('dopbin1')
  dopbin2 = getextrai('dopbin2')
  delbin1 = getextrai('delbin1')
  delbin2 = getextrai('delbin2')
endelse
if isnull(delayunit) and notnull(baud) then begin
  if isnull(spb) then spb = 1L
  if isnull(rows_per_baud) then rows_per_baud = spb
  delayunit = baud/rows_per_baud
endif
if notnull(eph_col) and notnull(eph_row) and notnull(dfreq) and notnull(delayunit) then begin
  f = dfreq*(findgen(width) - eph_col)       ; Doppler (Hz)
  d = delayunit*(findgen(height) - eph_row)  ; delay (usec)
endif else begin
  print,"ERROR in showlimi: Need the 'eph_col' and 'eph_row' and 'fres' tags,"
  print,"                   plus either the 'delayunit' or 'baudlen' tag"
  return
endelse

; Display the Doppler and delay limits

if notnull(dopbin1) and notnull(dopbin2) then begin
  flim1 = dfreq*(dopbin1 - 0.5 - eph_col)
  flim2 = dfreq*(dopbin2 + 0.5 - eph_col)
  printstring1 = 'Doppler bins ' + string(dopbin1,format='(i0)') + '-'  $
                 + string(dopbin2,format='(i0)')
endif else begin
  printstring1 = ''
endelse
if notnull(delbin1) and notnull(delbin2) then begin
  dlim1 = delayunit*(delbin1 - 0.5 - eph_row)
  dlim2 = delayunit*(delbin2 + 0.5 - eph_row)
  printstring2 = 'delay   bins ' + string(delbin1,format='(i0)') + '-'  $
                 + string(delbin2,format='(i0)')
endif else begin
  printstring2 = ''
endelse
if notnull(printstring1) and notnull(printstring2) then begin
  nskip = strlen(printstring1) - strlen(printstring2)
  if nskip gt 0 then begin
    for i=0L,nskip-1 do printstring2 = printstring2 + ' '
  endif else if nskip lt 0 then begin
    for i=0L,-nskip-1 do printstring1 = printstring1 + ' '
  endif
endif
if notnull(printstring1) then begin
  printstring1 = printstring1 + '  (' + string(flim1,format='(f8.2)')   $
                 + ' Hz to ' + string(flim2,format='(f8.2)') + ' Hz)'
endif else begin
  printstring1 = 'Doppler bins undefined'
endelse
if notnull(printstring2) then begin
  printstring2 = printstring2 + '  (' + string(dlim1,format='(f8.2)')   $
                 + ' us to ' + string(dlim2,format='(f8.2)') + ' us)'
endif else begin
  printstring2 = 'delay   bins undefined'
endelse

print,' '
print,printstring1
print,printstring2
print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro flip,chan=chan,stack=st,zeropivot=zeropivot,iqerror=iqerror,help=help

; Flip (reverse the direction of) the loaded pair.  If the stack keyword
; is set, do this instead for all pairs in the stack.
;
; By default, the flipping is done pivoting on the center of the spectrum,
; the xjcen tag is adjusted (if the center isn't zero Doppler), the posfr tag
; is negated, and (for the loaded pair) the frequency vector is also flipped
; (pivoting on its center).
;
; If zeropivot is set, flipping is done pivoting on zero Doppler (so xjcen
; stays the same).  posfr is negated, and (for the loaded pair) the frequency
; vector is also flipped (about zero Doppler).
;
; If iqerror is set, it's assumed that we're correcting a cabling error
; (switching I and Q), so the flipping is done pivoting on the middle of the
; spectrum, and xjcen, posfr, and the frequency vector are left unchanged.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if not keyword_set(zeropivot) then zeropivot = 0
if not keyword_set(iqerror) then iqerror = 0

if n_params() ne 0 or keyword_set(help) or (zeropivot and iqerror) or $
              (n_elements(chan) gt 0 and not keyword_set(st) and not iqerror) then begin
  print,' '
  print,'flip[,chan=1 or 2][,/stack][,/zeropivot OR ,/iqerror][,/help]'
  print,' '
  print,'    Reverse the direction of the loaded pair or (if /stack is set)'
  print,'        of every pair in the stack'
  print,' '
  print,'    If neither /zeropivot nor /iqerror is set, flipping is done about'
  print,'        the center of the spectrum, the xjcen tag is adjusted accordingly,'
  print,'        the posfr tag is negated, and (for the loaded pair) the frequency'
  print,'        vector is also flipped'
  print,' '
  print,'    If /zeropivot is set, flipping is done about zero Doppler, xjcen is'
  print,'        left unchanged, posfr is negated, and frequencies are flipped'
  print,' '
  print,'    /iqerror is used to correct a cabling error (switching I and Q):'
  print,'        flipping is done about the middle of the spectrum, and neither'
  print,'        xjcen, posfr, nor the frequency vector are changed'
  print,' '
  print,"    (You can't set both /zeropivot and /iqerror)"
  print,' '
  print,'    chan can only be specified if /stack or /iqerror is set, since both'
  print,'        channels of the loaded pair share one frequency vector.'
  print,'        If flipping one channel would leave the two channels of the pair'
  print,'        with different xjcen or posfr tags, flip will leave the channel'
  print,'        unchanged and will print a warning message.'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Flip the loaded pair

  if (*loaded).ndata le 2 then begin
    print,"ERROR in flip: There's no loaded pair to flip"
    return
  endif
  nuse = 1L
endif else begin

  ; Flip every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in flip: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; Loop through and process each spectrum

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack[n] : *loaded
  tags = useStruc.tags
  xjcen = tags[0].xjcen
  posfr = tags[0].posfr
  ndata = useStruc.ndata
  old_offset1 = tags.xjcen - tags.jsnr1
  old_offset2 = tags.jsnr2 - tags.xjcen
  npol = n_elements(tags)

; Check which channels to flip

  flippol = intarr(npol)
  if n_elements(chan) eq 0 then begin
    flippol = flippol + 1
  endif else if chan ge 1 && chan le npol then begin
    if (npol gt 2) then begin
      print, 'ERROR in flip: flipping individual channels in Stokes data'
      print, "doesn't work, as the cross terms don't match, sorry."
      return
    endif
    flippol[chan-1] = 1
  endif else begin
    print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
    return
  endelse
  
  ; Create temporary copies of flipped spectra and (perhaps) tags

  for ch=1,npol do begin
    if flippol[ch-1] then begin
      origspec = reform(useStruc.spec[ch-1,*])
      nshift = (zeropivot) ? (2*xjcen - ndata + 1L) : 0L
      flipspec = shift(reverse(origspec), nshift)
      useStruc.spec[ch-1,*] = flipspec

      ; If specified, set xjcen (the bin whose CENTER is zero Doppler)
      ; and posfr (+1 for frequency increasing rightward, -1 for leftward)

      if not zeropivot and not iqerror then useStruc.tags[ch-1].xjcen = ndata - xjcen - 1L
      if not iqerror then useStruc.tags[ch-1].posfr = -posfr

      ; Flip the signal limits

      useStruc.tags[ch-1].jsnr1 = $
          (useStruc.tags[ch-1].xjcen - old_offset2[ch-1]) > 0L
      useStruc.tags[ch-1].jsnr2 = $
          (useStruc.tags[ch-1].xjcen + old_offset1[ch-1]) < (ndata - 1)
    endif
  endfor

  ; Check that the OC and SC channels will still have the same
  ; xjcen and posfr tags; print a warning and refuse to change anything
  ; if this isn't the case. If it's Stokes, all bets are off so don't
  ; bother figuring it out.
  ;
  ; If everything is OK, flip the frequency vector if specified
  ; and replace the old structure with the new one

  xjcen_new = useStruc.tags.xjcen
  posfr_new = useStruc.tags.posfr
  if (xjcen_new[0] eq xjcen_new[1]) and (posfr_new[0] eq posfr_new[1]) then begin
    if not keyword_set(st) and not iqerror then begin
      dfreq = tags[0].dfreq
      f = posfr_new[0]*dfreq*(findgen(ndata) - xjcen_new[0])
      useStruc.freq = f
    endif
    if keyword_set(st) then *stack[n] = useStruc else *loaded = useStruc
  endif else begin
    print,'WARNING in flip: Pair #',n+1,' left unchanged',format='(a,i0,a)'
    print,'                 OC vs SC would have different xjcen (',     $
          xjcen_new[0],' vs ',xjcen_new[1],') or posfr (',              $
          posfr_new[0],' vs ',posfr_new[1],') tags',format='(4(a,i0),a)'
  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro flip1,stack=st,zeropivot=zeropivot,iqerror=iqerror,help=help

; Flip (reverse the direction of) the loaded single-channel spectrum.
; If the stack keyword is set, do this instead for all spectra in the
; single-channel stack.
;
; By default, the flipping is done pivoting on the center of the spectrum,
; the xjcen tag is adjusted (if the center isn't zero Doppler), the posfr tag
; is negated, and (for the loaded single-channel spectrum) the frequency
; vector is also flipped (pivoting on its center).
;
; If zeropivot is set, flipping is done pivoting on zero Doppler (so xjcen
; stays the same).  posfr is negated, and (for the loaded single-channel
; spectrum) the frequency vector is also flipped (about zero Doppler).
;
; If iqerror is set, it's assumed that we're correcting a cabling error
; (switching I and Q), so the flipping is done pivoting on the middle of the
; spectrum, and xjcen, posfr, and the frequency vector are left unchanged.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if not keyword_set(zeropivot) then zeropivot = 0
if not keyword_set(iqerror) then iqerror = 0

if n_params() ne 0 or keyword_set(help) or (zeropivot and iqerror) then begin
  print,' '
  print,'flip1[,/stack][,/zeropivot OR ,/iqerror][,/help]'
  print,' '
  print,'     Reverse the direction of the loaded single-channel spectrum or'
  print,'         (if /stack is set) of every spectrum in the single-channel stack'
  print,' '
  print,'     If neither /zeropivot nor /iqerror is set, flipping is done about'
  print,'         the center of the spectrum, the xjcen tag is adjusted accordingly,'
  print,'         the posfr tag is negated, and (for the loaded single-channel spectrum)'
  print,'         the frequency vector is also flipped'
  print,' '
  print,'     If /zeropivot is set, flipping is done about zero Doppler, xjcen is'
  print,'         left unchanged, posfr is negated, and frequencies are flipped'
  print,' '
  print,'     /iqerror is used to correct a cabling error (switching I and Q):'
  print,'         flipping is done about the middle of the spectrum, and neither'
  print,'         xjcen, posfr, nor the frequency vector are changed'
  print,' '
  print,"     (You can't set both /zeropivot and /iqerror)"
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Flip the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,"ERROR in flip1: There's no loaded single-channel spectrum to flip"
    return
  endif
  nuse = 1L
endif else begin

  ; Flip all spectra in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in flip1: The single-channel stack is empty'
    return
  endif

  nuse = nstack1
endelse

; Loop through and process each spectrum

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack1[n] : *loaded1
  tags1 = useStruc.tags
  xjcen = tags1.xjcen
  posfr = tags1.posfr
  ndata = useStruc.ndata
  old_offset1 = tags1.xjcen - tags1.jsnr1
  old_offset2 = tags1.jsnr2 - tags1.xjcen

  ; Create temporary copies of flipped spectrum and (perhaps) tags

  origspec = useStruc.spec
  nshift = (zeropivot) ? (2*xjcen - ndata + 1L) : 0
  flipspec = shift(reverse(origspec), nshift)
  useStruc.spec = flipspec

  ; If specified, set xjcen (the bin whose CENTER is zero Doppler)
  ; and posfr (+1 for frequency increasing rightward, -1 for leftward)

  if not zeropivot and not iqerror then useStruc.tags.xjcen = ndata - xjcen - 1L
  if not iqerror then useStruc.tags.posfr = -posfr

  ; Flip the signal limits

  useStruc.tags.jsnr1 = (useStruc.tags.xjcen - old_offset2) > 0L
  useStruc.tags.jsnr2 = (useStruc.tags.xjcen + old_offset1) < (ndata - 1)

  ; If this is the loaded single-channel spectrum and the frequency
  ; vector needs to be changed, do it now

  if not keyword_set(st) and not iqerror then begin
    xjcen_new = useStruc.tags.xjcen
    posfr_new = useStruc.tags.posfr
    dfreq = tags1.dfreq
    f = posfr_new*dfreq*(findgen(ndata) - xjcen_new)
    useStruc.freq = f
  endif

  ; Upload the revised structure

  if keyword_set(st) then *stack1[n] = useStruc else *loaded1 = useStruc

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro flipi,stack=st,zeropivot=zeropivot,iqerror=iqerror,help=help

; Flip (reverse the Doppler direction of) the loaded image.  If the stack keyword
; is set, do this instead for all images in the stack.
;
; By default, the flipping is done pivoting on the center of the spectrum
; and the eph_col tag is adjusted (if the center isn't zero Doppler)
;
; If zeropivot is set, flipping is done pivoting on the center of the zero-Doppler
; column (so eph_col stays almost the same).
;
; If iqerror is set, it's assumed that we're correcting a cabling error
; (switching I and Q), so the flipping is done pivoting on the middle of the
; spectrum, and eph_col is left unchanged.
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if not keyword_set(zeropivot) then zeropivot = 0
if not keyword_set(iqerror) then iqerror = 0

if n_params() ne 0 or keyword_set(help) or (zeropivot and iqerror) or $
              (n_elements(chan) gt 0 and not keyword_set(st) and not iqerror) then begin
  print,' '
  print,'flipi[,/stack][,/zeropivot OR ,/iqerror][,/help]'
  print,' '
  print,'    Reverse the Doppler direction of the loaded image or (if /stack is set)'
  print,'        of every image in the stack'
  print,' '
  print,'    If neither /zeropivot nor /iqerror is set, flipping is done about'
  print,'        the center of the image and the eph_col tag is adjusted accordingly'
  print,' '
  print,'    If /zeropivot is set, flipping is done about the center of the 0-Hz column,'
  print,'        so that eph_col changes by only a fraction of a column'
  print,' '
  print,'    /iqerror is used to correct a cabling error (switching I and Q):'
  print,'        flipping is done about the center of the image, and eph_col is'
  print,'        left unchanged'
  print,' '
  print,"    (You can't set both /zeropivot and /iqerror)"
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Flip the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,' '
    print,'ERROR in flipi: No image is loaded'
    print,' '
    return
  endif

  ; Get new 0-Hz column and work out the number of columns
  ; we'll need to shift the image after reversing it

  width = (*loadedi).width
  orig_eph_col = getextrai('eph_col')
  if keyword_set(iqerror) then begin
    flip_eph_col = orig_eph_col
    nshift = 0L
  endif else if keyword_set(zeropivot) then begin
    if isnull(orig_eph_col) then begin
      print,' '
      print,"ERROR in flipi: loaded image has no 'eph_col' tag"
      print,' '
      return
    endif
    flip_eph_col = 2*round(orig_eph_col) - orig_eph_col
    setextrai,'f','eph_col',flip_eph_col
    nshift = 2*round(orig_eph_col) - width + 1L
  endif else begin
    if isnull(orig_eph_col) then begin
      print,"WARNING in flipi: loaded image has no 'eph_col' tag"
      flip_eph_col = ''
    endif else begin
      flip_eph_col = width - 1 - orig_eph_col
      setextrai,'f','eph_col',flip_eph_col
    endelse
    nshift = 0L
  endelse

  ; Read the image, reverse it in Doppler, and shift it to put
  ; the center of the 0-Hz column in the desired place

  origimage = (*loadedi).image
  flipimage = shift(reverse(origimage, 1), [nshift,0])
  (*loadedi).image = flipimage

  ; Get new signal limits

  if not keyword_set(iqerror) then begin
    orig_dopbin1 = getextrai('dopbin1')
    orig_dopbin2 = getextrai('dopbin2')
    if notnull(orig_eph_col) and $
                   notnull(dopbin1) and notnull(dopbin2) then begin
      orig_offset1 = orig_eph_col - orig_dopbin1
      orig_offset2 = orig_dopbin2 - orig_eph_col
      flip_dopbin1 = round(flip_eph_col - orig_offset2) > 0L
      flip_dopbin2 = round(flip_eph_col + orig_offset1) < (width - 1)
      setextrai,'i','dopbin1',flip_dopbin1
      setextrai,'i','dopbin2',flip_dopbin2
    endif
  endif

  ; Adjust the DC column

  if not keyword_set(iqerror) then begin
    orig_dc_col = getextrai('dc_col')
    if notnull(orig_eph_col) and notnull(orig_dc_col) then begin
      orig_dc_offset = orig_dc_col - orig_eph_col
      flip_dc_col = flip_eph_col - orig_dc_offset
      setextrai,'f','dc_col',flip_dc_col
    endif
  endif

  ; Adjust the FITS "crpix1" tag

  if not keyword_set(iqerror) then begin
    orig_crpix1 = getextrai('crpix1')
    if notnull(orig_eph_col) and notnull(orig_crpix1) then begin
      flip_crpix1 = flip_eph_col + 1
      setextrai,'f','crpix1',flip_crpix1
    endif
  endif

  ; Add an extra tag to show that flipping has been performed

  if keyword_set(iqerror) then begin
    setextrai,'s','iqerror','corrected'
  endif else begin
    flipped = getextrai('flipped')
    if isnull(flipped) then begin
      setextrai,'s','flipped','true'
    endif else begin
      deleteextrai,'flipped',/silent
    endelse
  endelse

endif else begin

  ; Flip all images in the stack

  if nstacki eq 0 then begin
    print,'ERROR in flipi: The image stack is empty'
    return
  endif

  for n=0L,nstacki-1 do begin
    width = (*stacki[n]).width
    orig_eph_col = getextrai('eph_col',stack=(n+1))
    if keyword_set(iqerror) then begin
      flip_eph_col = orig_eph_col
      nshift = 0L
    endif else if keyword_set(zeropivot) then begin
      if isnull(orig_eph_col) then begin
        print,' '
        print,"ERROR in flipi: stacki image #",n+1, $
              " has no 'eph_col' tag",format='(a,i0,a)'
        print,' '
        return
      endif
      flip_eph_col = 2*round(orig_eph_col) - orig_eph_col
      setextrai,'f','eph_col',flip_eph_col,stack=(n+1)
      nshift = 2*round(eph_col) - width + 1L
    endif else begin
      if isnull(orig_eph_col) then begin
        print,"WARNING in flipi: stacki image #",n+1, $
              " has no 'eph_col' tag",format='(a,i0,a)'
        flip_eph_col = ''
      endif else begin
        flip_eph_col = width - 1 - orig_eph_col
        setextrai,'f','eph_col',flip_eph_col,stack=(n+1)
      endelse
      nshift = 0L
    endelse
    origimage = (*stacki[n]).image
    flipimage = shift(reverse(origimage, 1), [nshift,0])
    (*stacki[n]).image = flipimage
    if not keyword_set(iqerror) then begin
      orig_dopbin1 = getextrai('dopbin1',stack=(n+1))
      orig_dopbin2 = getextrai('dopbin2',stack=(n+1))
      if notnull(orig_eph_col) and $
                     notnull(dopbin1) and notnull(dopbin2) then begin
        orig_offset1 = orig_eph_col - orig_dopbin1
        orig_offset2 = orig_dopbin2 - orig_eph_col
        flip_dopbin1 = round(flip_eph_col - orig_offset2) > 0L
        flip_dopbin2 = round(flip_eph_col + orig_offset1) < (width - 1)
        setextrai,'i','dopbin1',flip_dopbin1,stack=(n+1)
        setextrai,'i','dopbin2',flip_dopbin2,stack=(n+1)
      endif
    endif
    if not keyword_set(iqerror) then begin
      orig_dc_col = getextrai('dc_col',stack=(n+1))
      if notnull(orig_eph_col) and notnull(orig_dc_col) then begin
        orig_dc_offset = orig_dc_col - orig_eph_col
        flip_dc_col = flip_eph_col - orig_dc_offset
        setextrai,'f','dc_col',flip_dc_col,stack=(n+1)
      endif
    endif
    if not keyword_set(iqerror) then begin
      orig_crpix1 = getextrai('crpix1',stack=(n+1))
      if notnull(orig_eph_col) and notnull(orig_crpix1) then begin
        flip_crpix1 = flip_eph_col + 1
        setextrai,'f','crpix1',flip_crpix1,stack=(n+1)
      endif
    endif
    if keyword_set(iqerror) then begin
      setextrai,'s','iqerror','corrected',stack=(n+1)
    endif else begin
      flipped = getextrai('flipped',stack=(n+1))
      if isnull(flipped) then begin
        setextrai,'s','flipped','true',stack=(n+1)
      endif else begin
        deleteextrai,'flipped',/silent,stack=(n+1)
      endelse
    endelse
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro renormalize,omitf=farray,bins=bins,chan=chan,stack=st,spikeremove=spikeremove, $
                srthresh=srthresh,setrmsm=setrmsm,nochange=nochange,help=help

; Renormalize the loaded pair to have zero mean and unit standard deviation.
; The signal region (jsnr1-jsnr2) is ignored; additional regions to ignore can
; be specified via the omitf keyword.
;
; /stack does this instead for all pairs in the stack
;
; If /spikeremove is set, turn on an algorithm to ignore spikes; see help message
; for function polymaskfit for details.
;
; srthresh specifies how strong (how many sigmas) such spikes must be to be ignored
; (default = 5)
;
; If /setrmsm is set, set the rmsm cw tag (ratio of actual to predicted rms scatter)
; to 1.0 after renormalization

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that the sets of frequency ranges to be omitted from the noise calculation
; have been correctly input (if input at all).  Define the number of such
; regions to be the number input plus one: one extra for the signal region
; specified via tags jsnr1 - jsnr2, which is always omitted from the calculation.

if n_elements(farray) eq 0 then begin
  n_omitregions = 1
  omitOK = 1
endif else begin
  omitsize = size(farray)
  omitdims = omitsize[0]
  dim1 = omitsize[1]
  if omitdims eq 0 or omitdims gt 2 or dim1 ne 2 then begin
    omitOK = 0
  endif else begin
    omitOK = 1
    n_omitregions = (omitdims eq 1) ? 2 : (1 + omitsize[2])
    for k=0L,n_omitregions-2 do begin
      if farray[0,k] gt farray[1,k] then omitOK = 0
    endfor
  endelse
endelse

removing_spikes = keyword_set(spikeremove)
if n_elements(srthresh) eq 0 then srthresh = 5.0
good_srthresh = srthresh gt 0

if n_params() ne 0 or not omitOK or keyword_set(help) then begin
  print,' '
  print,'renormalize[,omitf=farray[,/bins]][,chan=1 or 2][,/stack]  $'
  print,'           [,/spikeremove[,srthresh=srthresh]][,/setrmsm]  $'
  print,'           [,/nochange][,/help]'
  print,' '
  print,'Renormalize the loaded pair to have zero mean and unit standard deviation'
  print,'outside of the signal limits (tags jsnr1 - jsnr2)'
  print,' '
  print,'If desired you can specify additional regions to omit from the noise'
  print,"          calculation by giving those regions' Doppler limits in farray"
  print,' '
  print,'farray is an array of Doppler frequency pairs defining regions'
  print,'          (IN ADDITION to the signal region) which should be omitted from'
  print,'          the noise calculation.  Values are in Hz, unless /bins is set, in'
  print,'          which case they are bin numbers (counting from 0).  Values must'
  print,'          be in increasing order even if the spectrum has frequency'
  print,'          increasing leftward.  Examples of farray are [-150,-120]'
  print,'          or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  print,'/bins sets the units of farray entries as spectral bins, counting from 0;'
  print,'          if omitted, the units are Hz.'
  print,' '
  print,'/stack renormalizes all pairs in the stack rather than the loaded pair'
  print,' '
  print,'/spikeremove turns on an algorithm to ignore spikes;'
  print,'    call "dummy = polymaskfit(/help)" to see details'
  print,' '
  print,'srthresh specifies how strong (in sigmas) such spikes must be to be ignored'
  print,'   (default = 5)'
  print,' '
  print,'/nochange simply lists means and standard deviations without'
  print,'          actually changing the pair(s)'
  print,' '
  print,'/setrmsm sets the rmsm tag equal to 1.0, unless /nochange is set, in which'
  print,'          case it sets it to the calculated standard deviation'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Renormalize the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in renormalize: No pair is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Renormalize every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in renormalize: The stack is empty'
    return
  endif
  nuse = nstack
endelse

; Determine which polarization channels were specified

if n_elements(chan) eq 0 then begin
  dopol = [1,1]
  first_channel = 1
endif else if chan eq 1 then begin
  dopol = [1,0]
  first_channel = 1
endif else if chan eq 2 then begin
  dopol = [0,1]
  first_channel = 2
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse
polstring = ['    OC: ','    SC: ']

print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this pair

  useStruc = (keyword_set(st)) ? *stack[n] : *loaded
  pair = useStruc.spec
  tags = useStruc.tags
  posfr = tags[0].posfr
  xjcen = tags[0].xjcen
  dfreq = tags[0].dfreq
  jsnr1 = tags.jsnr1
  jsnr2 = tags.jsnr2
  ndata = useStruc.ndata
  if keyword_set(st) then begin
    freq = posfr*dfreq*(findgen(ndata) - xjcen)
  endif else begin
    freq = useStruc.freq
  endelse
  npol = n_elements(tags)
  if npol gt 2 then begin
    print, "ERROR in renormalize: Not yet implemented for Stokes data"
    return
  end

  ; Add the signal region to the list of regions to be omitted
  ; (keep track of the two polarization channels separately)

  if keyword_set(bins) then begin
    farray_use = lonarr(2,2,n_omitregions)   ; (chan, min/max, iregion)
    farray_use[*,0,0] = jsnr1
    farray_use[*,1,0] = jsnr2
    if n_omitregions gt 1 then begin
      farray_use[0,*,1:n_omitregions-1] = farray
      farray_use[1,*,1:n_omitregions-1] = farray
    endif
  endif else begin
    farray_use = fltarr(2,2,n_omitregions)
    if posfr eq 1 then begin
      farray_use[*,0,0] = (jsnr1 - xjcen)*dfreq
      farray_use[*,1,0] = (jsnr2 - xjcen)*dfreq
    endif else begin
      farray_use[*,0,0] = -(jsnr2 - xjcen)*dfreq
      farray_use[*,1,0] = -(jsnr1 - xjcen)*dfreq
    endelse
    if n_omitregions gt 1 then begin
      farray_use[0,*,1:n_omitregions-1] = farray
      farray_use[1,*,1:n_omitregions-1] = farray
    endif
  endelse

  ; Loop over channels

  printstring = 'Pair #' + string(n+1L,format='(i0)') + ':'
  for ch=1,2 do begin
    if dopol[ch-1] then begin
      spec = double(reform(pair[ch-1,*]))  ;  change back to float later

      ; Initialize a mask vector
  
      in_mask = intarr(ndata) + 1

      ; Mask out the frequency bin ranges to be omitted from the fit

      for k=0L,n_omitregions-1 do begin
        if keyword_set(bins) then begin
          omit_left = round(1.0*farray_use[ch-1,0,k])
          omit_right = round(1.0*farray_use[ch-1,1,k])
        endif else if posfr eq 1 then begin
          omit_left = round(xjcen + farray_use[ch-1,0,k]/dfreq)
          omit_right = round(xjcen + farray_use[ch-1,1,k]/dfreq)
        endif else begin
          omit_left = round(xjcen - farray_use[ch-1,1,k]/dfreq)
          omit_right = round(xjcen - farray_use[ch-1,0,k]/dfreq)
        endelse
        if omit_left lt 0 or omit_right ge ndata then begin 
          print,' '
          if ch eq first_channel then begin
            print,'WARNING for pair #',n+1,': Excluded region (bins ',omit_left,'-', $
                  omit_right,') extends beyond the spectrum (bins 0-',ndata-1,')',   $
                  format='(4(a,i0),a)'
          endif
        endif
        if omit_left lt ndata and omit_right ge 0 then begin
          omit_left = omit_left > 0
          omit_right = omit_right < (ndata - 1)
          in_mask[omit_left:omit_right] = 0
          if ch eq first_channel then begin
            print,'Pair #',n+1,':    excluding bins ',omit_left,'-', $
                  omit_right,' from mean and sdev calculation',format='(3(a,i0),a)'
          endif
        endif
      endfor

      ; Check that at least two bins have NOT been excluded

      if total(in_mask) lt 2 then begin
        print,'Pair #',n+1,' channel ',ch,' has < 2 signal bins left: no changes made', $
              format='(a,i0,a,i0,a)'
      endif else begin

        ; Compute the mean over non-signal bins

        specmean = polymaskfit(freq,spec,0,in_mask=in_mask,spikeremove=removing_spikes, $
                               srthresh=srthresh,yerror=specrms)
        if keyword_set(nochange) then begin
          printstring = printstring + polstring[ch-1] + 'mean = '  $
                        + string(specmean,format='(f6.3)') + ' and rms = '  $
                        + string(specrms,format='(f5.3)')
        endif else begin
          printstring = printstring + polstring[ch-1] + 'subtract mean = '  $
                        + string(specmean,format='(f6.3)') + ', divide by rms = '  $
                        + string(specrms,format='(f5.3)')
        endelse

        ; Save the normalized spectrum if requested

        if keyword_set(nochange) then begin
          if keyword_set(st) then begin
            if keyword_set(setrmsm) then (*stack[n]).tags[ch-1].rmsm = float(specrms)
          endif else begin
            if keyword_set(setrmsm) then (*loaded).tags[ch-1].rmsm = float(specrms)
          endelse
        endif else begin
          spec = float((spec - specmean[0])/specrms)
          if keyword_set(st) then begin
            (*stack[n]).spec[ch-1,*] = spec
            if keyword_set(setrmsm) then (*stack[n]).tags[ch-1].rmsm = 1.0
          endif else begin
            (*loaded).spec[ch-1,*] = spec
            if keyword_set(setrmsm) then (*loaded).tags[ch-1].rmsm = 1.0
          endelse
        endelse

      endelse
    endif
  endfor
  print,printstring

endfor

print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro renormalize1,omitf=farray,bins=bins,stack=st,spikeremove=spikeremove, $
                 srthresh=srthresh,setrmsm=setrmsm,nochange=nochange,help=help

; Renormalize the loaded single-channel spectrum to have zero mean and unit
; standard deviation.  The signal region (jsnr1-jsnr2) is ignored; additional
; regions to ignore can be specified via the omitf keyword.
;
; /stack does this instead for all spectra in the single-channel stack
;
; If /spikeremove is set, turn on an algorithm to ignore spikes; see help message
; for function polymaskfit for details.
;
; srthresh specifies how strong (how many sigmas) such spikes must be to be ignored
; (default = 5)
;
; If /setrmsm is set, set the rmsm cw tag (ratio of actual to predicted rms scatter)
; to 1.0 after renormalization

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that the sets of frequency ranges to be omitted from the noise calculation
; have been correctly input (if input at all).  Define the number of such
; regions to be the number input plus one: one extra for the signal region
; specified via tags jsnr1 - jsnr2, which is always omitted from the calculation.

if n_elements(farray) eq 0 then begin
  n_omitregions = 1
  omitOK = 1
endif else begin
  omitsize = size(farray)
  omitdims = omitsize[0]
  dim1 = omitsize[1]
  if omitdims eq 0 or omitdims gt 2 or dim1 ne 2 then begin
    omitOK = 0
  endif else begin
    omitOK = 1
    n_omitregions = (omitdims eq 1) ? 2 : (1 + omitsize[2])
    for k=0L,n_omitregions-2 do begin
      if farray[0,k] gt farray[1,k] then omitOK = 0
    endfor
  endelse
endelse

removing_spikes = keyword_set(spikeremove)
if n_elements(srthresh) eq 0 then srthresh = 5.0
good_srthresh = srthresh gt 0

if n_params() ne 0 or not omitOK or keyword_set(help) then begin
  print,' '
  print,'renormalize1[,omitf=farray[,/bins]][,/stack][,/spikeremove[,srthresh=srthresh]]  $'
  print,'            [,/setrmsm][,/nochange][,/help]'
  print,' '
  print,'Renormalize the loaded single-channel spectrum to have zero mean and unit'
  print,'          standard deviation outside of the signal limits (tags jsnr1 - jsnr2)'
  print,' '
  print,'If desired you can specify additional regions to omit from the noise'
  print,"          calculation by giving those regions' Doppler limits in farray"
  print,' '
  print,'farray is an array of Doppler frequency pairs defining regions'
  print,'          (IN ADDITION to the signal region) which should be omitted from'
  print,'          the noise calculation.  Values are in Hz, unless /bins is set, in'
  print,'          which case they are bin numbers (counting from 0).  Values must'
  print,'          be in increasing order even if the spectrum has frequency'
  print,'          increasing leftward.  Examples of farray are [-150,-120]'
  print,'          or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  print,'/bins sets the units of farray entries as spectral bins, counting from 0;'
  print,'          if omitted, the units are Hz.'
  print,' '
  print,'/stack renormalizes all spectra in the single-channel stack rather than the'
  print,'          loaded spectrum'
  print,' '
  print,'/spikeremove turns on an algorithm to ignore spikes;'
  print,'    call "dummy = polymaskfit(/help)" to see details'
  print,' '
  print,'srthresh specifies how strong (in sigmas) such spikes must be to be ignored'
  print,'   (default = 5)'
  print,' '
  print,'/nochange simply lists means and standard deviations without'
  print,'          actually changing the spectrum/spectra'
  print,' '
  print,'/setrmsm sets the rmsm tag equal to 1.0, unless /nochange is set, in which'
  print,'          case it sets it to the calculated standard deviation'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Renormalize the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in renormalize1: No single-channel spectrum is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Renormalize every spectrum in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in renormalize1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
endelse

print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack1[n] : *loaded1
  spec = useStruc.spec
  tags1 = useStruc.tags
  posfr = tags1.posfr
  xjcen = tags1.xjcen
  dfreq = tags1.dfreq
  jsnr1 = tags1.jsnr1
  jsnr2 = tags1.jsnr2
  ndata = useStruc.ndata
  if keyword_set(st) then begin
    freq = posfr*dfreq*(findgen(ndata) - xjcen)
  endif else begin
    freq = useStruc.freq
  endelse

  ; Add the signal region to the list of regions to be omitted

  if keyword_set(bins) then begin
    farray_use = lonarr(2,n_omitregions)   ; (min/max, iregion)
    farray_use[0,0] = jsnr1
    farray_use[1,0] = jsnr2
    if n_omitregions gt 1 then farray_use[*,1:n_omitregions-1] = farray
  endif else begin
    farray_use = fltarr(2,n_omitregions)
    if posfr eq 1 then begin
      farray_use[0,0] = (jsnr1 - xjcen)*dfreq
      farray_use[1,0] = (jsnr2 - xjcen)*dfreq
    endif else begin
      farray_use[0,0] = -(jsnr2 - xjcen)*dfreq
      farray_use[1,0] = -(jsnr1 - xjcen)*dfreq
    endelse
    if n_omitregions gt 1 then farray_use[*,1:n_omitregions-1] = farray
  endelse

  spec = double(spec)  ;  change back to float later

  ; Initialize a mask vector
  
  in_mask = intarr(ndata) + 1

  ; Mask out the frequency bin ranges to be omitted from the fit

  for k=0L,n_omitregions-1 do begin
    if keyword_set(bins) then begin
      omit_left = round(1.0*farray_use[0,k])
      omit_right = round(1.0*farray_use[1,k])
    endif else if posfr eq 1 then begin
      omit_left = round(xjcen + farray_use[0,k]/dfreq)
      omit_right = round(xjcen + farray_use[1,k]/dfreq)
    endif else begin
      omit_left = round(xjcen - farray_use[1,k]/dfreq)
      omit_right = round(xjcen - farray_use[0,k]/dfreq)
    endelse
    if omit_left lt 0 or omit_right ge ndata then begin 
      print,' '
      print,'WARNING for spectrum #',n+1,': Excluded region (bins ',omit_left,'-', $
            omit_right,') extends beyond the spectrum (bins 0-',ndata-1,')',   $
            format='(4(a,i0),a)'
    endif
    if omit_left lt ndata and omit_right ge 0 then begin
      omit_left = omit_left > 0
      omit_right = omit_right < (ndata - 1)
      in_mask[omit_left:omit_right] = 0
      print,'Spectrum #',n+1,':    excluding bins ',omit_left,'-', $
            omit_right,' from mean and sdev calculation',format='(3(a,i0),a)'
    endif
  endfor

  ; Check that at least two bins have NOT been excluded

  if total(in_mask) lt 2 then begin
    print,'Spectrum #',n+1,' has < 2 signal bins left: no changes made', $
          format='(a,i0,a)'
  endif else begin

    ; Compute the mean over non-signal bins

    specmean = polymaskfit(freq,spec,0,in_mask=in_mask,spikeremove=removing_spikes, $
                           srthresh=srthresh,yerror=specrms)
    printstring = 'Spectrum #' + string(n+1L,format='(i0)') + ':    '
    if keyword_set(nochange) then begin
      printstring = printstring + 'mean = ' + string(specmean,format='(f6.3)')  $
                    + ' and rms = ' + string(specrms,format='(f5.3)')
    endif else begin
      printstring = printstring + 'subtract mean = '  $
                    + string(specmean,format='(f6.3)') + ', divide by rms = '  $
                    + string(specrms,format='(f5.3)')
    endelse
    print,printstring

    ; Save the normalized spectrum if requested

    if keyword_set(nochange) then begin
      if keyword_set(st) then begin
        if keyword_set(setrmsm) then (*stack1[n]).tags.rmsm = float(specrms)
      endif else begin
        if keyword_set(setrmsm) then (*loaded1).tags.rmsm = float(specrms)
      endelse
    endif else begin
      spec = float((spec - specmean[0])/specrms)
      if keyword_set(st) then begin
        (*stack1[n]).spec = spec
        if keyword_set(setrmsm) then (*stack1[n]).tags.rmsm = 1.0
      endif else begin
        (*loaded1).spec = spec
        if keyword_set(setrmsm) then (*loaded1).tags.rmsm = 1.0
      endelse
    endelse

  endelse

endfor

print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro renormalizei,delomit=delarray,dopomit=doparray,bins=bins,stack=st,nochange=nochange, $
                 help=help

; Renormalize the loaded image to have zero mean and unit standard deviation for the
; region outside of the signal limits (extra tags dopbin1, dopbin2, delbin1, delbin2)
;
; /stack does this instead for all images in the image stack
;
; Modified 6/26/06 to add the delomit=delarray and dopomit=doparray keyword options
; for omitting delay Doppler ranges from the noise calculation (in addition to the
; signal region); this was spurred by the need to mask out the secondary component
; of a binary asteroid

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that the sets of delay and Doppler ranges to be omitted from the noise
; calculation have been correctly input (if input at all).  Define the number of such
; delay regions to be the number input plus one: one extra for the signal region
; specified via extra tags delbin1 - delbin2, which is always omitted from the
; calculation.  Similarly for Doppler regions, since signal region dopbin1 - dopbin2
; is always omitted.

if n_elements(delarray) eq 0 then begin
  n_omitdelregions = 1
  omitdelOK = 1
endif else begin
  omitsize = size(delarray)
  omitdims = omitsize[0]
  dim1 = omitsize[1]
  if omitdims eq 0 or omitdims gt 2 or dim1 ne 2 then begin
    omitdelOK = 0
  endif else begin
    omitdelOK = 1
    n_omitdelregions = (omitdims eq 1) ? 2 : (1 + omitsize[2])
    for k=0L,n_omitdelregions-2 do begin
      if delarray[0,k] gt delarray[1,k] then omitdelOK = 0
    endfor
  endelse
endelse

if n_elements(doparray) eq 0 then begin
  n_omitdopregions = 1
  omitdopOK = 1
endif else begin
  omitsize = size(doparray)
  omitdims = omitsize[0]
  dim1 = omitsize[1]
  if omitdims eq 0 or omitdims gt 2 or dim1 ne 2 then begin
    omitdopOK = 0
  endif else begin
    omitdopOK = 1
    n_omitdopregions = (omitdims eq 1) ? 2 : (1 + omitsize[2])
    for k=0L,n_omitdopregions-2 do begin
      if doparray[0,k] gt doparray[1,k] then omitdopOK = 0
    endfor
  endelse
endelse

if n_params() ne 0 or not (omitdelOK and omitdopOK) $
                   or n_omitdelregions ne n_omitdopregions or keyword_set(help) then begin
  print,' '
  print,'renormalizei[,delomit=delarray,dopomit=doparray[,/bins]][,/stack][,/nochange][,/help]'
  print,' '
  print,'Renormalize the loaded image to have zero mean and unit standard deviation'
  print,'outside of the signal limits (extra tags dopbin1, dopbin2, delbin1, delbin2)'
  print,' '
  print,'If desired you can specify additional rectangular regions to omit from the noise'
  print,"          calculation by giving those regions' delay limits in delarray and"
  print,'           their Doppler limits in doparray'
  print,' '
  print,'delarray is an array of delay pairs defining regions (IN ADDITION to the signal region)'
  print,'          which should be omitted from the noise calculation.  Values are in usec,'
  print,'          unless /bins is set, in which case they are row numbers (counting from 0).'
  print,'          Examples of delarray are [-40,-5] or [[-40,-5], [57,123.2], [306,422]].'
  print,' '
  print,'doparray is an array of Doppler frequency pairs defining regions'
  print,'          (IN ADDITION to the signal region) which should be omitted from'
  print,'          the noise calculation.  Values are in Hz, unless /bins is set, in'
  print,'          which case they are column numbers (counting from 0).  Examples of'
  print,'          doparray are [-150,-120] or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  print,'NOTE: delarray and doparray must have the SAME number of pairs'
  print,' '
  print,'/bins sets the units of delarray and doparray entries as image rows and columns,'
  print,'          counting from 0; if omitted, the units are usec and Hz.'
  print,' '
  print,'/stack renormalizes all images in the image stack rather than the loaded image'
  print,' '
  print,'/nochange simply lists means and standard deviations without'
  print,'          actually changing the image(s)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Set signal limits for the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,'ERROR in renormalizei: No image is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Set signal limits for every image in the stack

  if nstacki eq 0 then begin
    print,'ERROR in renormalizei: The image stack is empty'
    return
  endif
  nuse = nstacki
endelse

print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this image

  if keyword_set(st) then begin
    image = (*stacki[n]).image
    width = (*stacki[n]).width
    height = (*stacki[n]).height
    dfreq = getextrai('fres',stack=(n+1))
    eph_col = getextrai('eph_col',stack=(n+1))
    eph_row = getextrai('eph_row',stack=(n+1))
    delayunit = getextrai('delayunit',stack=(n+1))
    spb = getextrai('samples_per_baud',stack=(n+1))
    rows_per_baud = getextrai('rows_per_baud',stack=(n+1))
    baud = getextrai('baudlen',stack=(n+1))
    dopbin1 = getextrai('dopbin1',stack=(n+1))
    dopbin2 = getextrai('dopbin2',stack=(n+1))
    delbin1 = getextrai('delbin1',stack=(n+1))
    delbin2 = getextrai('delbin2',stack=(n+1))
    rowsPerBaud = getextrai('rows_per_baud',stack=(n+1))
  endif else begin
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
    dopbin1 = getextrai('dopbin1')
    dopbin2 = getextrai('dopbin2')
    delbin1 = getextrai('delbin1')
    delbin2 = getextrai('delbin2')
    rowsPerBaud = getextrai('rows_per_baud')
  endelse
  if isnull(delayunit) and notnull(baud) then begin
    if isnull(spb) then spb = 1L
    if isnull(rows_per_baud) then rows_per_baud = spb
    delayunit = baud/rows_per_baud
  endif
  if isnull(dopbin1) then begin
    dopbin1 = 0L
    print,'dopbin1 missing for image #',n+1,', set to ',dopbin1, $
          format='(a,i0,a,i0)'
  endif else if dopbin1 lt 0 then begin
    dopbin1 = 0L
    print,'dopbin1 < 0 for image #',n+1,', reset to ',dopbin1, $
          format='(a,i0,a,i0)'
  endif
  if isnull(dopbin2) then begin
    dopbin2 = width - 1
    print,'dopbin2 missing for image #',n+1,', set to ',dopbin2, $
          format='(a,i0,a,i0)'
  endif else if dopbin2 ge width then begin
    dopbin2 = width - 1
    print,'dopbin2 >= width for image #',n+1,', reset to ',dopbin2, $
          format='(a,i0,a,i0)'
  endif
  if isnull(delbin1) then begin
    delbin1 = 0L
    print,'delbin1 missing for image #',n+1,', set to ',delbin1, $
          format='(a,i0,a,i0)'
  endif else if delbin1 lt 0 then begin
    delbin1 = 0L
    print,'delbin1 < 0 for image #',n+1,', reset to ',delbin1, $
          format='(a,i0,a,i0)'
  endif
  if isnull(delbin2) then begin
    delbin2 = height - 1
    print,'delbin2 missing for image #',n+1,', set to ',delbin2, $
          format='(a,i0,a,i0)'
  endif else if delbin2 ge height then begin
    delbin2 = height - 1
    print,'delbin2 >= height for image #',n+1,', reset to ',delbin2, $
          format='(a,i0,a,i0)'
  endif
  if isnull(rowsPerBaud) then rowsPerBaud = 1L
  nbaud = height/rowsPerBaud

  ; If there's any non-signal region, do the renormalization

  if (dopbin1 eq 0) and (dopbin2 eq width - 1) and (delbin1 eq 0)  $
                    and (delbin2 eq height - 1) then begin
    print,'Image #',n+1,' is all signal: no changes made',format='(a,i0,a)'
  endif else if dopbin2 lt dopbin1 then begin
    print,'Image #',n+1,' has dopbin2 < dopbin1: no changes made', $
          format='(a,i0,a)'
  endif else if delbin2 lt delbin1 then begin
    print,'Image #',n+1,' has delbin2 < delbin1: no changes made', $
          format='(a,i0,a)'
  endif else if nbaud*rowsPerBaud ne height then begin
    print,'Image #',n+1,' has a noninteger number of bauds: no changes made', $
          format='(a,i0,a)'
  endif else begin

    ; Add the signal region to the list of regions to be omitted

    if keyword_set(bins) then begin
      delarray_use = lonarr(2,n_omitdelregions)
      delarray_use[*,0] = [delbin1,delbin2]
      if n_omitdelregions gt 1 then delarray_use[*,1:n_omitdelregions-1] = delarray
      doparray_use = lonarr(2,n_omitdopregions)
      doparray_use[*,0] = [dopbin1,dopbin2]
      if n_omitdopregions gt 1 then doparray_use[*,1:n_omitdopregions-1] = doparray
    endif else begin
      delarray_use = fltarr(2,n_omitdelregions)
      delarray_use[*,0] = ([delbin1,delbin2] - eph_row)*delayunit
      if n_omitdelregions gt 1 then delarray_use[*,1:n_omitdelregions-1] = delarray
      doparray_use = fltarr(2,n_omitdopregions)
      doparray_use[*,0] = ([dopbin1,dopbin2] - eph_col)*dfreq
      if n_omitdopregions gt 1 then doparray_use[*,1:n_omitdopregions-1] = doparray
    endelse

    ; Initialize a mask array

    mask = intarr(width, height) + 1

    ; Mask out the frequency bin ranges to be omitted from the fit

    for k=0L,n_omitdelregions-1 do begin
      if keyword_set(bins) then begin
        omit_del1 = round(1.0*delarray_use[0,k])
        omit_del2 = round(1.0*delarray_use[1,k])
        omit_dop1 = round(1.0*doparray_use[0,k])
        omit_dop2 = round(1.0*doparray_use[1,k])
      endif else begin
        omit_del1 = round(eph_row + delarray_use[0,k]/delayunit)
        omit_del2 = round(eph_row + delarray_use[1,k]/delayunit)
        omit_dop1 = round(eph_col + doparray_use[0,k]/dfreq)
        omit_dop2 = round(eph_col + doparray_use[1,k]/dfreq)
      endelse
      if omit_del1 lt 0 or omit_del2 ge height $
                        or omit_dop1 lt 0 or omit_dop2 ge width then begin 
        print,' '
        print,'WARNING for image #',n+1,': Excluded region (cols ', $
              omit_dop1,'-',omit_dop2,', rows ',omit_del1,'-',omit_del2, $
              ') extends beyond the image (cols 0-',width-1,', rows 0-', $
              height-1,')',format='(7(a,i0),a)'
      endif
      if omit_del1 lt height and omit_del2 ge 0 $
                             and omit_dop1 lt width and omit_dop2 ge 0 then begin
        omit_del1 = omit_del1 > 0
        omit_del2 = omit_del2 < (height - 1)
        omit_dop1 = omit_dop1 > 0
        omit_dop2 = omit_dop2 < (width - 1)
        mask[omit_dop1:omit_dop2,omit_del1:omit_del2] = 0
        print,'Image #',n+1,': excluding rectangular region (cols ', $
              omit_dop1,'-',omit_dop2,', rows ',omit_del1,'-', $
              omit_del2,') from mean and sdev calculation', $
              format='(a,i0,a,i0,a,i0,a,i0,a,i0,a)'
      endif
    endfor

    ; Check that at least two pixels have NOT been excluded

    if total(mask) lt 2 then begin
      print,'Image #',n+1,' has < 2 signal pixels left: no changes made', $
            format='(a,i0,a)'
    endif else begin

      ; Rearrange the image by spb; also convert it to double precision
      ; (but output float at the end)
      ;
      ; Also rearrange the mask

      image = 1.0D*reform(temporary(image), width, rowsPerBaud, nbaud)
      mask = reform(temporary(mask), width, rowsPerBaud, nbaud)

      ; Compute the mean over non-signal pixels for each row per baud

      actionstring = (keyword_set(nochange)) ? '' : 'subtract '
      for i = 0L,rowsPerBaud-1 do begin
        image_mean_1rpb = total(image[*,i,*]*mask[*,i,*]) / total(mask[*,i,*])
        image[*,i,*] = image[*,i,*] - image_mean_1rpb
        print,'Image #',n+1,': ',actionstring,'mean  =',image_mean_1rpb, $
              ' for spb ',i+1,format='(a,i0,3a,e13.5,a,i0)'
      endfor

      ; Scale the image to the standard deviation of the non-signal region

      actionstring = (keyword_set(nochange)) ? '' : 'divide by '
      image_sdev = sqrt( total((image*mask)^2) / (total(mask) - 1) )
      image = temporary(image) / image_sdev
      print,'Image #',n+1,': ',actionstring,'sdev =',image_sdev, $
            format='(a,i0,3a,e13.5)'
      print,' '

      ; Rearrange image so that different rows per baud are interleaved again,
      ; and convert it back to single precision

      image = float(reform(temporary(image), width, rowsPerBaud*nbaud))

      ; Save the normalized image

      if not keyword_set(nochange) then begin
        if keyword_set(st) then begin
          (*stacki[n]).image = image
        endif else begin
          (*loadedi).image = image
        endelse
      endif

    endelse

  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro meanrms,help=help,_extra=_ext

; Compute the mean value (bias) and rms scatter about that bias for the loaded pair,
; unless /stack is set, in which case this is done for every pair in the pair stack.
; The signal region (jsnr1-jsnr2) is ignored; additional regions to ignore can
; be specified via the omitf keyword.
;
; /stack does this instead for all pairs in the stack
;
; If /spikeremove is set, turn on an algorithm to ignore spikes; see help message
; for function polymaskfit for details.
;
; srthresh specifies how strong (how many sigmas) such spikes must be to be ignored
; (default = 5)
;
; If /setrmsm is set, set the rmsm cw tag (ratio of actual to predicted rms scatter)
; equal to the computed rms scatter
;
; 2006 Jun 25: frequency refers to center of bin, not left edge
; 2006 Jun 29: move this code into the "renormalize" routine and use that routine
;              to do the work for meanrms

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,' '
  print,'meanrms[,omitf=farray[,/bins]][,chan=1 or 2][,/stack]  $'
  print,'       [,/spikeremove[,srthresh=srthresh]][,/setrmsm][,/help]'
  print,' '
  print,'Compute the mean spectral value and rms scatter about the mean'
  print,'outside of the signal limits (tags jsnr1 - jsnr2) for the loaded pair'
  print,' '
  print,'If desired you can specify additional regions to omit from the noise'
  print,"          calculation by giving those regions' Doppler limits in farray"
  print,' '
  print,'farray is an array of Doppler frequency pairs defining regions'
  print,'          (IN ADDITION to the signal region) which should be omitted from'
  print,'          the noise calculation.  Values are in Hz, unless /bins is set, in'
  print,'          which case they are bin numbers (counting from 0).  Values must'
  print,'          be in increasing order even if the spectrum has frequency'
  print,'          increasing leftward.  Examples of farray are [-150,-120]'
  print,'          or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  print,'/bins sets the units of farray entries as spectral bins, counting from 0;'
  print,'          if omitted, the units are Hz.'
  print,' '
  print,'/stack renormalizes all pairs in the stack rather than the loaded pair'
  print,' '
  print,'/spikeremove turns on an algorithm to ignore spikes;'
  print,'    call "dummy = polymaskfit(/help)" to see details'
  print,' '
  print,'srthresh specifies how strong (in sigmas) such spikes must be to be ignored'
  print,'   (default = 5)'
  print,' '
  print,'/setrmsm sets the rmsm tag equal to the computed rms scatter'
  print,' '
  return
endif

renormalize,/nochange,_extra=_ext

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro meanrms1,help=help,_extra=_ext

; Compute the mean value (bias) and rms scatter about that bias for the loaded,
; single-channel spectrum, unless /stack is set, in which case this is done for
; every spectrum in the single-channel stack.  The signal region (jsnr1-jsnr2) is
; ignored; additional regions to ignore can be specified via the omitf keyword.
;
; /stack does this instead for all spectra in the single-channel stack
;
; If /spikeremove is set, turn on an algorithm to ignore spikes; see help message
; for function polymaskfit for details.
;
; srthresh specifies how strong (how many sigmas) such spikes must be to be ignored
; (default = 5)
;
; If /setrmsm is set, set the rmsm cw tag (ratio of actual to predicted rms scatter)
; equal to the computed rms scatter
;
; 2006 Jun 25: frequency refers to center of bin, not left edge
; 2006 Jun 29: move this code into the "renormalize1" routine and use that routine
;              to do the work for meanrms1

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,' '
  print,'meanrms1[,omitf=farray[,/bins]][,/stack][,/spikeremove[,srthresh=srthresh]]  $'
  print,'        [,/setrmsm][,/help]'
  print,' '
  print,'Compute the mean spectral value and rms scatter about the mean outside of the'
  print,'signal limits (tags jsnr1 - jsnr2) for the loaded single-channel spectrum'
  print,' '
  print,'If desired you can specify additional regions to omit from the noise'
  print,"          calculation by giving those regions' Doppler limits in farray"
  print,' '
  print,'farray is an array of Doppler frequency pairs defining regions'
  print,'          (IN ADDITION to the signal region) which should be omitted from'
  print,'          the noise calculation.  Values are in Hz, unless /bins is set, in'
  print,'          which case they are bin numbers (counting from 0).  Values must'
  print,'          be in increasing order even if the spectrum has frequency'
  print,'          increasing leftward.  Examples of farray are [-150,-120]'
  print,'          or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  print,'/bins sets the units of farray entries as spectral bins, counting from 0;'
  print,'          if omitted, the units are Hz.'
  print,' '
  print,'/stack renormalizes all spectra in the single-channel stack rather than the'
  print,'          loaded spectrum'
  print,' '
  print,'/spikeremove turns on an algorithm to ignore spikes;'
  print,'    call "dummy = polymaskfit(/help)" to see details'
  print,' '
  print,'srthresh specifies how strong (in sigmas) such spikes must be to be ignored'
  print,'   (default = 5)'
  print,' '
  print,'/setrmsm sets the rmsm tag equal to the computed rms scatter'
  print,' '
  return
endif

renormalize1,/nochange,_extra=_ext

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro meanrmsi,help=help,_extra=_ext

; Compute the mean value (bias) and rms scatter about that bias for the loaded
; image, unless /stack is set, in which case this is done for every image in
; the image stack.  The signal region -- a rectangle defined by columns
; dopbin1-dopbin2 and rows delbin1-delbin2, counting from zero -- is ignored.
; Other regions to ignore can be specified via the delomit and dopomit keywords.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,' '
  print,'meanrmsi[,delomit=delarray,dopomit=doparray][,/bins][,/stack][,/help]'
  print,' '
  print,'Compute the mean pixel value (bias) and rms scatter about the mean'
  print,'    for the loaded image, unless /stack is set, in which case'
  print,'    this is done for every image in the image stack'
  print,' '
  print,'delarray is an array of delay pairs defining regions (IN ADDITION to the signal region)'
  print,'          which should be omitted from the noise calculation.  Values are in usec,'
  print,'          unless /bins is set, in which case it is in rows (counting from 0).  Examples'
  print,'          of delarray are [-40,-5] or [[-40,-5], [57,123.2], [306,422]].'
  print,' '
  print,'doparray is an array of Doppler frequency pairs defining regions'
  print,'          (IN ADDITION to the signal region) which should be omitted from'
  print,'          the noise calculation.  Values are in Hz, unless /bins is set, in'
  print,'          which case it is in columns (counting from 0).  Examples of doparray'
  print,'          are [-150,-120] or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  print,'NOTE: delarray and doparray must have the SAME number of pairs'
  print,' '
  print,'/bins sets the units of delarray and doparray entries as image rows and columns,'
  print,'          counting from 0; if omitted, the units are usec and Hz.'
  print,' '
  return
endif

renormalizei,/nochange,_extra=_ext

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro resetphase,period=period,jd0=jd0,phase0=phase0,stack=st,help=help

; Reset the rotation phase tag for the loaded pair, or else for every
; pair in the stack if the /stack keyword is set.
;
; The rotation period, the epoch (jd0) at which phase = phase0, and the
; value of phase0 can be given as keywords if different from the value(s)
; already in the extra tags; the extra tags are changed accordingly.
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

keywordsOK = 1
if n_elements(period) gt 0 then begin
  keywordsOK = keywordsOK and (period gt 0)
endif else begin
  period = -9.99
endelse
if n_elements(jd0) gt 0 then begin
  keywordsOK = keywordsOK and (jd0 gt 0) and (size(jd0, /type) eq 5)
endif else begin
  jd0 = -9.99
endelse
if n_elements(phase0) gt 0 then begin
  keywordsOK = keywordsOK and (phase0 ge 0) and (phase0 lt 360)
endif else begin
  phase0 = -9.99
endelse

if n_params() gt 0 or not keywordsOK or keyword_set(help) then begin
  print,' '
  print,'resetphase[,period=period][,jd0=jd0][,phase0=phase0][,/stack][,/help]'
  print,' '
  print,'Reset the rotation phase tag of the loaded pair, or else of every'
  print,'    pair in the stack if /stack is set'
  print,' '
  print,'Using one or more of the following will cause the corresponding'
  print,'    extra tags to be changed before the new phase is computed:'
  print,' '
  print,'    period = new rotation period (hr)'
  print,'    jd0    = epoch at which rotation phase = phase0 (Julian date)'
  print,'    phase0 = rotation phase when epoch = jd0 (deg)'
  print,' '
  print,'    (jd0 must be double-precision; must have 0 <= phase0 < 360)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Reset rotation phase for the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in resetphase: No pair is loaded'
    return
  endif
  nuse = 1L
  period_use = fltarr(nuse)
  jd0_use = dblarr(nuse)
  phase0_use = fltarr(nuse)
  jdmean_use = dblarr(nuse)

  if period lt 0 then begin
    tempperiod = getextra('period')
    if isnull(tempperiod) then begin
      print,"ERROR in resetphase: The loaded pair has no 'period' extra tag"
      return
    endif
    period_use[0] = tempperiod
  endif else begin
    period_use[0] = period
  endelse
  if jd0 lt 0 then begin
    tempjd0 = getextra('jd0')
    if isnull(tempjd0) then begin
      print,"ERROR in resetphase: The loaded pair has no 'jd0' extra tag"
      return
    endif
    jd0_use[0] = tempjd0
  endif else begin
    jd0_use[0] = jd0
  endelse
  if phase0 lt 0 then begin
    tempphase0 = getextra('phase0')
    if isnull(tempphase0) then begin
      print,"ERROR in resetphase: The loaded pair has no 'phase0' extra tag"
      return
    endif
    phase0_use[0] = tempphase0
  endif else begin
    phase0_use[0] = phase0
  endelse

  tempjdmean = getextra('jdmean')
  if isnull(tempjdmean) then begin
    print,"ERROR in resetphase: The loaded pair has no 'jdmean' extra tag"
    return
  endif
  jdmean_use[0] = tempjdmean

  ; Compute the new rotation phase

  phase = (jdmean_use - jd0_use)*(24.0/period_use)*360 + phase0_use
  phase = phase - 360*floor(phase/360)

  ; Assign the new rotation phase to the tags;
  ; also assign the new period, jd0, and phase0 to the extra tags
  ; if any of these were specified as keywords

  settag,'phase',phase[0],/silent
  if period gt 0 then setextra,'f','period',strtrim(period, 2)
  if jd0 gt 0 then setextra,'d','jd0',string(jd0,format='(d13.5)')
  if phase0 ge 0 then setextra,'f','phase0',strtrim(phase0, 2)

endif else begin

  ; Reset rotation phase for every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in resetphase: The pair stack is empty'
    return
  endif
  nuse = nstack
  period_use = fltarr(nuse)
  jd0_use = dblarr(nuse)
  phase0_use = fltarr(nuse)
  jdmean_use = dblarr(nuse)

  for n=1,nuse do begin
    if period lt 0 then begin
      tempperiod = getextra('period',stack=n)
      if isnull(tempperiod) then begin
        print,"ERROR in resetphase: stack pair #",n," has no 'period' extra tag", $
              format='(a,i0,a)'
        return
      endif
      period_use[n-1] = tempperiod
    endif else begin
      period_use[n-1] = period
    endelse
    if jd0 lt 0 then begin
      tempjd0 = getextra('jd0',stack=n)
      if isnull(tempjd0) then begin
        print,"ERROR in resetphase: stack pair #",n," has no 'jd0' extra tag", $
              format='(a,i0,a)'
        return
      endif
      jd0_use[n-1] = tempjd0
    endif else begin
      jd0_use[n-1] = jd0
    endelse
    if phase0 lt 0 then begin
      tempphase0 = getextra('phase0',stack=n)
      if isnull(tempphase0) then begin
        print,"ERROR in resetphase: stack pair #",n," has no 'phase0' extra tag", $
              format='(a,i0,a)'
        return
      endif
      phase0_use[n-1] = tempphase0
    endif else begin
      phase0_use[n-1] = phase0
    endelse

    tempjdmean = getextra('jdmean',stack=n)
    if isnull(tempjdmean) then begin
      print,"ERROR in resetphase: stack pair #",n," has no 'jdmean' extra tag", $
            format='(a,i0,a)'
      return
    endif
    jdmean_use[n-1] = tempjdmean

  endfor

  phase = (jdmean_use - jd0_use)*(24.0/period_use)*360 + phase0_use
  phase = phase - 360*floor(phase/360)

  for n=1,nuse do begin
    settag,'phase',phase[n-1],/silent,stack=n
    if period gt 0 then setextra,'f','period',strtrim(period, 2),stack=n
    if jd0 gt 0 then setextra,'d','jd0',string(jd0,format='(d13.5)'),stack=n
    if phase0 ge 0 then setextra,'f','phase0',strtrim(phase0, 2),stack=n
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro resetphase1,period=period,jd0=jd0,phase0=phase0,stack=st,help=help

; Reset the rotation phase tag for the loaded single-channel spectrum,
; or else for every spectrum in the single-channel stack if the /stack
; keyword is set.
;
; The rotation period, the epoch (jd0) at which phase = phase0, and the
; value of phase0 can be given as keywords if different from the value(s)
; already in the extra tags; the extra tags are changed accordingly.
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

keywordsOK = 1
if n_elements(period) gt 0 then begin
  keywordsOK = keywordsOK and (period gt 0)
endif else begin
  period = -9.99
endelse
if n_elements(jd0) gt 0 then begin
  keywordsOK = keywordsOK and (jd0 gt 0) and (size(jd0, /type) eq 5)
endif else begin
  jd0 = -9.99
endelse
if n_elements(phase0) gt 0 then begin
  keywordsOK = keywordsOK and (phase0 ge 0) and (phase0 lt 360)
endif else begin
  phase0 = -9.99
endelse

if n_params() gt 0 or not keywordsOK or keyword_set(help) then begin
  print,' '
  print,'resetphase1[,period=period][,jd0=jd0][,phase0=phase0][,/stack][,/help]'
  print,' '
  print,'Reset the rotation phase tag of the loaded single-channel spectrum,'
  print,'    or else of every spectrum in the single-channel stack if /stack is set'
  print,' '
  print,'Using one or more of the following will cause the corresponding'
  print,'    extra tags to be changed before the new phase is computed:'
  print,' '
  print,'    period = new rotation period (hr)'
  print,'    jd0    = epoch at which rotation phase = phase0 (Julian date)'
  print,'    phase0 = rotation phase when epoch = jd0 (deg)'
  print,' '
  print,'    (jd0 must be double-precision; must have 0 <= phase0 < 360)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Reset rotation phase for the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in resetphase1: No single-channel spectrum is loaded'
    return
  endif
  nuse = 1L
  period_use = fltarr(nuse)
  jd0_use = dblarr(nuse)
  phase0_use = fltarr(nuse)
  jdmean_use = dblarr(nuse)

  if period lt 0 then begin
    tempperiod = getextra1('period')
    if isnull(tempperiod) then begin
      print,"ERROR in resetphase1: The loaded single-channel spectrum has no 'period' extra tag"
      return
    endif
    period_use[0] = tempperiod
  endif else begin
    period_use[0] = period
  endelse
  if jd0 lt 0 then begin
    tempjd0 = getextra1('jd0')
    if isnull(tempjd0) then begin
      print,"ERROR in resetphase1: The loaded single-channel spectrum has no 'jd0' extra tag"
      return
    endif
    jd0_use[0] = tempjd0
  endif else begin
    jd0_use[0] = jd0
  endelse
  if phase0 lt 0 then begin
    tempphase0 = getextra1('phase0')
    if isnull(tempphase0) then begin
      print,"ERROR in resetphase1: The loaded single-channel spectrum has no 'phase0' extra tag"
      return
    endif
    phase0_use[0] = tempphase0
  endif else begin
    phase0_use[0] = phase0
  endelse

  tempjdmean = getextra1('jdmean')
  if isnull(tempjdmean) then begin
    print,"ERROR in resetphase1: The loaded single-channel spectrum has no 'jdmean' extra tag"
    return
  endif
  jdmean_use[0] = tempjdmean

  ; Compute the new rotation phase

  phase = (jdmean_use - jd0_use)*(24.0/period_use)*360 + phase0_use
  phase = phase - 360*floor(phase/360)

  ; Assign the new rotation phase to the tags;
  ; also assign the new period, jd0, and phase0 to the extra tags
  ; if any of these were specified as keywords

  settag1,'phase',phase[0],/silent
  if period gt 0 then setextra1,'f','period',strtrim(period, 2)
  if jd0 gt 0 then setextra1,'d','jd0',string(jd0,format='(d13.5)')
  if phase0 ge 0 then setextra1,'f','phase0',strtrim(phase0, 2)

endif else begin

  ; Reset rotation phase for every spectrum in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in resetphase1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
  period_use = fltarr(nuse)
  jd0_use = dblarr(nuse)
  phase0_use = fltarr(nuse)
  jdmean_use = dblarr(nuse)

  for n=1,nuse do begin
    if period lt 0 then begin
      tempperiod = getextra1('period',stack=n)
      if isnull(tempperiod) then begin
        print,"ERROR in resetphase1: stack1 spectrum #",n," has no 'period' extra tag", $
              format='(a,i0,a)'
        return
      endif
      period_use[n-1] = tempperiod
    endif else begin
      period_use[n-1] = period
    endelse
    if jd0 lt 0 then begin
      tempjd0 = getextra1('jd0',stack=n)
      if isnull(tempjd0) then begin
        print,"ERROR in resetphase1: stack1 spectrum #",n," has no 'jd0' extra tag", $
              format='(a,i0,a)'
        return
      endif
      jd0_use[n-1] = tempjd0
    endif else begin
      jd0_use[n-1] = jd0
    endelse
    if phase0 lt 0 then begin
      tempphase0 = getextra1('phase0',stack=n)
      if isnull(tempphase0) then begin
        print,"ERROR in resetphase1: stack1 spectrum #",n," has no 'phase0' extra tag", $
              format='(a,i0,a)'
        return
      endif
      phase0_use[n-1] = tempphase0
    endif else begin
      phase0_use[n-1] = phase0
    endelse

    tempjdmean = getextra1('jdmean',stack=n)
    if isnull(tempjdmean) then begin
      print,"ERROR in resetphase1: stack1 spectrum #",n," has no 'jdmean' extra tag", $
            format='(a,i0,a)'
      return
    endif
    jdmean_use[n-1] = tempjdmean

  endfor

  phase = (jdmean_use - jd0_use)*(24.0/period_use)*360 + phase0_use
  phase = phase - 360*floor(phase/360)

  for n=1,nuse do begin
    settag1,'phase',phase[n-1],/silent,stack=n
    if period gt 0 then setextra1,'f','period',strtrim(period, 2),stack=n
    if jd0 gt 0 then setextra1,'d','jd0',string(jd0,format='(d13.5)'),stack=n
    if phase0 ge 0 then setextra1,'f','phase0',strtrim(phase0, 2),stack=n
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro resetphasei,period=period,jd0=jd0,phase0=phase0,stack=st,help=help

; Reset the rotation phase extra tag for the loaded image, or else for
; every image in the image stack if the /stack keyword is set.
;
; The rotation period, the epoch (jd0) at which phase = phase0, and the
; value of phase0 can be given as keywords if different from the value(s)
; already in the extra tags; the extra tags are changed accordingly.
 
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

keywordsOK = 1
if n_elements(period) gt 0 then begin
  keywordsOK = keywordsOK and (period gt 0)
endif else begin
  period = -9.99
endelse
if n_elements(jd0) gt 0 then begin
  keywordsOK = keywordsOK and (jd0 gt 0) and (size(jd0, /type) eq 5)
endif else begin
  jd0 = -9.99
endelse
if n_elements(phase0) gt 0 then begin
  keywordsOK = keywordsOK and (phase0 ge 0) and (phase0 lt 360)
endif else begin
  phase0 = -9.99
endelse

if n_params() gt 0 or not keywordsOK or keyword_set(help) then begin
  print,' '
  print,'resetphasei[,period=period][,jd0=jd0][,phase0=phase0][,/stack][,/help]'
  print,' '
  print,'Reset the rotation phase extra tag of the loaded image, or else of every image'
  print,'    in the image stack if /stack is set'
  print,' '
  print,'Using one or more of the following will cause the corresponding'
  print,'    extra tags to be changed before the new phase is computed:'
  print,' '
  print,'    period = new rotation period (hr)'
  print,'    jd0    = epoch at which rotation phase = phase0 (Julian date)'
  print,'    phase0 = rotation phase when epoch = jd0 (deg)'
  print,' '
  print,'    (jd0 must be double-precision; must have 0 <= phase0 < 360)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Reset rotation phase for the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,'ERROR in resetphasei: No image is loaded'
    return
  endif
  nuse = 1L
  period_use = fltarr(nuse)
  jd0_use = dblarr(nuse)
  phase0_use = fltarr(nuse)
  jdmean_use = dblarr(nuse)

  if period lt 0 then begin
    tempperiod = getextrai('period')
    if isnull(tempperiod) then begin
      print,"ERROR in resetphasei: The loaded image has no 'period' extra tag"
      return
    endif
    period_use[0] = tempperiod
  endif else begin
    period_use[0] = period
  endelse
  if jd0 lt 0 then begin
    tempjd0 = getextrai('jd0')
    if isnull(tempjd0) then begin
      print,"ERROR in resetphasei: The loaded image has no 'jd0' extra tag"
      return
    endif
    jd0_use[0] = tempjd0
  endif else begin
    jd0_use[0] = jd0
  endelse
  if phase0 lt 0 then begin
    tempphase0 = getextrai('phase0')
    if isnull(tempphase0) then begin
      print,"ERROR in resetphasei: The loaded image has no 'phase0' extra tag"
      return
    endif
    phase0_use[0] = tempphase0
  endif else begin
    phase0_use[0] = phase0
  endelse

  tempjdmean = getextrai('jdmean')
  if isnull(tempjdmean) then begin
    print,"ERROR in resetphasei: The loaded image has no 'jdmean' extra tag"
    return
  endif
  jdmean_use[0] = tempjdmean

  ; Compute the new rotation phase

  phase = (jdmean_use - jd0_use)*(24.0/period_use)*360 + phase0_use
  phase = phase - 360*floor(phase/360)

  ; Assign the new rotation phase to the extra tags;
  ; also assign the new period, jd0, and phase0 to the extra tags
  ; if any of these were specified as keywords

  setextrai,'f','phase',string(phase[0],format='(f8.4)'), $
            comment='mean rotation phase in deg'
  if period gt 0 then setextrai,'f','period',strtrim(period, 2)
  if jd0 gt 0 then setextrai,'d','jd0',string(jd0,format='(d13.5)')
  if phase0 ge 0 then setextrai,'f','phase0',strtrim(phase0, 2)

endif else begin

  ; Reset rotation phase for every image in the image stack

  if nstacki eq 0 then begin
    print,'ERROR in resetphasei: The image stack is empty'
    return
  endif
  nuse = nstacki
  period_use = fltarr(nuse)
  jd0_use = dblarr(nuse)
  phase0_use = fltarr(nuse)
  jdmean_use = dblarr(nuse)

  for n=1,nuse do begin
    if period lt 0 then begin
      tempperiod = getextrai('period',stack=n)
      if isnull(tempperiod) then begin
        print,"ERROR in resetphasei: stacki image #",n," has no 'period' extra tag", $
              format='(a,i0,a)'
        return
      endif
      period_use[n-1] = tempperiod
    endif else begin
      period_use[n-1] = period
    endelse
    if jd0 lt 0 then begin
      tempjd0 = getextrai('jd0',stack=n)
      if isnull(tempjd0) then begin
        print,"ERROR in resetphasei: stacki image #",n," has no 'jd0' extra tag", $
              format='(a,i0,a)'
        return
      endif
      jd0_use[n-1] = tempjd0
    endif else begin
      jd0_use[n-1] = jd0
    endelse
    if phase0 lt 0 then begin
      tempphase0 = getextrai('phase0',stack=n)
      if isnull(tempphase0) then begin
        print,"ERROR in resetphasei: stacki image #",n," has no 'phase0' extra tag", $
              format='(a,i0,a)'
        return
      endif
      phase0_use[n-1] = tempphase0
    endif else begin
      phase0_use[n-1] = phase0
    endelse

    tempjdmean = getextrai('jdmean',stack=n)
    if isnull(tempjdmean) then begin
      print,"ERROR in resetphasei: stacki image #",n," has no 'jdmean' extra tag", $
            format='(a,i0,a)'
      return
    endif
    jdmean_use[n-1] = tempjdmean

  endfor

  phase = (jdmean_use - jd0_use)*(24.0/period_use)*360 + phase0_use
  phase = phase - 360*floor(phase/360)

  for n=1,nuse do begin
    setextrai,'f','phase',string(phase[n-1],format='(f8.4)'), $
              comment='mean rotation phase in deg',stack=n
    if period gt 0 then setextrai,'f','period',strtrim(period, 2),stack=n
    if jd0 gt 0 then setextrai,'d','jd0',string(jd0,format='(d13.5)'),stack=n
    if phase0 ge 0 then setextrai,'f','phase0',strtrim(phase0, 2),stack=n
  endfor

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro extend,mode,lim1,lim2,normal=normal,stack=st,help=help

; Extend the loaded pair by adding extra bins with zero signal
;
; If keyword normal is set, use pseudorandom normal variates
;     (zero mean, unit variance) rather than zeros
;
; If keyword st is set, extend every pair in the stack instead
;
; 2006 Jun 25: frequency refers to center of bin, not left edge
; 2008 Mar 23: add "normal" keyword

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then mode = -1 ; just so it's defined

if keyword_set(help) or n_params() eq 0 or n_params() gt 3 or $
     mode lt 0 or mode gt 3 or $
     (mode eq 0 and (n_params() ne 1 or not keyword_set(st))) or $
     (mode eq 1 and n_params() eq 1) or $
     (mode eq 2 and n_params() lt 3) or $
     (mode eq 3 and n_params() eq 1) then begin
  print,' '
  print,'extend,mode[,lim1[,lim2]][,/stack][,/help]'
  print,' '
  print,'Extend the loaded pair by adding extra bins with zero signal, based on:'
  print,' '
  print,'  mode = 0: max freq range covered by the various stack pairs', $
        ' (must set /stack)',format='(2a)'
  print,'  mode = 1: frequencies lim1, lim2, or +/- lim1)'
  print,'  mode = 2: frequency bins lim1, lim2 (0-based)'
  print,'  mode = 3: frequency bin margins lim1, lim2 or lim1, lim1'
  print,' '
  print,'/normal extends the pair with pseudorandom normally distributed'
  print,'     zero-mean, unit-variance deviates rather than with all zeros'
  print,' '
  print,'/stack extends every stack pair instead of the loaded pair'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Extend the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in extend: No pair is loaded'
    return
  endif
  nuse = 1
endif else begin

  ; Extend every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in extend: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; Set the min and max frequencies (or frequency bins)

if n_params() eq 1 then begin

  ; mode = 0: Find the maximum frequency range covered by the various stack pairs
  
  lim1_use = 9.99e20
  lim2_use = -9.99e20
  for n=0L,nuse-1 do begin
    tags = (*stack[n]).tags
    dfreq = tags[0].dfreq
    posfr = tags[0].posfr
    xjcen = tags[0].xjcen
    ndata = (*stack[n]).ndata
    if posfr eq 1 then begin
      fmin = -dfreq*(xjcen + 0.5)
      fmax = dfreq*(ndata - xjcen - 0.5)
    endif else begin
      fmin = -dfreq*(ndata - xjcen - 0.5)
      fmax = dfreq*(xjcen + 0.5)
    endelse
    lim1_use = lim1_use < fmin
    lim2_use = lim2_use > fmax
  endfor
endif else if n_params() eq 2 then begin
  if (lim1 le 0) then begin
    print,'ERROR in extend: Must have lim1 > 0 if lim2 not specified'
    return
  endif
  lim1_use = (mode eq 1) ? -lim1 : lim1
  lim2_use = lim1
endif else begin
  if mode ne 3 and lim2 le lim1 then begin
    print,'ERROR in extend: Must have lim2 > lim1 for mode = 1 or 2'
    return
  endif else if mode eq 3 and (lim1 lt 0 or lim2 lt 0) then begin
    print,'ERROR in extend: Must have lim1, lim2 >= 0 for mode = 3'
    return
  endif
  lim1_use = lim1
  lim2_use = lim2
endelse

; Now process the pair(s)

for n=0L,nuse-1 do begin

  ; Get a few elements of this pair

  useStruc = (keyword_set(st)) ? *stack[n] : *loaded
  pair = useStruc.spec
  tags = useStruc.tags
  dfreq = tags[0].dfreq
  posfr = tags[0].posfr
  xjcen = tags[0].xjcen
  jsnr1 = tags.jsnr1
  jsnr2 = tags.jsnr2
  ndata = useStruc.ndata
  npol = n_elements(tags)
  f = (keyword_set(st)) ? posfr*dfreq*(findgen(ndata) - xjcen) : useStruc.freq
  
  ; Assign the left and right limits according to the
  ; method specified via the mode parameter

  if mode le 1 then begin

    ; Convert frequency limits to frequency bins

    if posfr eq 1 then begin
      bin1 = xjcen + round(lim1_use/dfreq)
      bin2 = xjcen + round(lim2_use/dfreq)
    endif else begin
      bin1 = xjcen - round(lim2_use/dfreq)
      bin2 = xjcen - round(lim1_use/dfreq)
    endelse
  endif else if mode eq 2 then begin

    ; Directly specify frequency bins

    bin1 = round(1.0*lim1_use)
    bin2 = round(1.0*lim2_use)
  endif else begin

    ; Extend frequency bins according to specified margins

    bin1 = -round(1.0*lim1_use)
    bin2 = ndata - 1 + round(1.0*lim2_use)
  endelse

  ; Make sure that the new limits don't fall within the spectrum
  ; and that more than two frequency bins remain after extending

  if bin1 gt 0 or bin2 lt (ndata-1) then begin
    print,'ERROR on pair #',n+1,': Requested frequency bins [', $
          bin1,'-',bin2,'] omit part of the spectrum [0-',ndata-1,']', $
          format='(4(a,i0),a)'
    print,'       Pair left unchanged'
  endif else if bin2 le bin1+1 then begin
    print,'ERROR on pair #',n+1,': Requested frequency bins [', $
          bin1,'-',bin2,'] would leave too few points in extended spectrum', $
          format='(3(a,i0),a)'
    print,'       Pair left unchanged'
  endif else begin
    newndata = bin2 - bin1 + 1
    newxjcen = xjcen - bin1
    newf = posfr*dfreq*(findgen(newndata) - newxjcen)
    tags.xjcen = newxjcen
    tags.jsnr1 = (jsnr1 - bin1) > 0L
    tags.jsnr2 = (jsnr2 - bin1) < (newndata - 1)
    if keyword_set(normal) then begin
      newpair = float(randomn(seed, npol, newndata, /double, /normal)) * $
        (tags.rmsm # replicate(1.0, newndata))
    endif else begin
      newpair = fltarr(npol, newndata)
    endelse
    newpair[*, (-bin1):(ndata-1-bin1)] = pair
    print,'Extending pair #',n+1,' to ',newndata,' points, frequency bins ',  $
          bin1,' - ',bin2,'  (',newf[0]-0.5*posfr*dfreq,' Hz to ',            $
          newf[newndata-1]+0.5*posfr*dfreq,' Hz)', $
          format='(a,i0,a,i0,a,i0,a,i0,a,f9.2,a,f9.2,a)'
    if keyword_set(st) then begin
      newStruc = {group:useStruc.group, spec:newpair, tags:tags, $
                  extratags:useStruc.extratags, ndata:newndata,  $
                  ntags:useStruc.ntags, nextra:useStruc.nextra, tname:useStruc.tname}
    endif else begin
      newStruc = {freq:newf, spec:newpair, tags:tags, extratags:useStruc.extratags, $
                  ndata:newndata, ntags:useStruc.ntags, nextra:useStruc.nextra,     $
                  tname:useStruc.tname}
    endelse

    if keyword_set(st) then *stack[n] = newStruc else *loaded = newStruc

  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro extend1,mode,lim1,lim2,normal=normal,stack=st,help=help

; Extend the loaded single-channel spectrum by adding extra bins with zero signal
;
; If keyword normal is set, use pseudorandom normal variates
;     (zero mean, unit variance) rather than zeros
;
; If keyword st is set, extend every spectrum in the single-channel stack instead
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then mode = -1 ; just so it's defined

if keyword_set(help) or n_params() eq 0 or n_params() gt 3 or $
     mode lt 0 or mode gt 3 or $
     (mode eq 0 and (n_params() ne 1 or not keyword_set(st))) or $
     (mode eq 1 and n_params() eq 1) or $
     (mode eq 2 and n_params() lt 3) or $
     (mode eq 3 and n_params() eq 1) then begin
  print,' '
  print,'extend1,mode[,lim1[,lim2]][,/stack][,/help]'
  print,' '
  print,'Extend the loaded single-channel spectrum by adding extra bins'
  print,'     with zero signal, based on:'
  print,' '
  print,'  mode = 0: max freq range covered by the various stack1 spectra', $
        ' (must set /stack)',format='(2a)'
  print,'  mode = 1: frequencies lim1, lim2, or +/- lim1'
  print,'  mode = 2: frequency bins lim1, lim2 (0-based)'
  print,'  mode = 3: frequency bin margins lim1, lim2 or lim1, lim1'
  print,' '
  print,'/normal extends the spectrum with pseudorandom normally distributed'
  print,'     zero-mean, unit-variance deviates rather than with all zeros'
  print,' '
  print,'/stack extends every stack1 spectrum instead of the ', $
        'loaded single-channel spectrum',format='(2a)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Extend the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in extend1: No single-channel spectrum is loaded'
    return
  endif
  nuse = 1
endif else begin

  ; Extend every spectrum in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in extend1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
endelse

; Set the min and max frequencies (or frequency bins)

if n_params() eq 1 then begin

  ; mode = 0: Find the maximum frequency range covered by the various stack1 spectra
  
  lim1_use = 9.99e20
  lim2_use = -9.99e20
  for n=0L,nuse-1 do begin
    tags = (*stack1[n]).tags
    dfreq = tags1.dfreq
    posfr = tags1.posfr
    xjcen = tags1.xjcen
    ndata = (*stack1[n]).ndata
    if posfr eq 1 then begin
      fmin = -dfreq*(xjcen + 0.5)
      fmax = dfreq*(ndata - xjcen - 0.5)
    endif else begin
      fmin = -dfreq*(ndata - xjcen - 0.5)
      fmax = dfreq*(xjcen + 0.5)
    endelse
    lim1_use = lim1_use < fmin
    lim2_use = lim2_use > fmax
  endfor
endif else if n_params() eq 2 then begin
  if (lim1 le 0) then begin
    print,'ERROR in extend1: Must have lim1 > 0 if lim2 not specified'
    return
  endif
  lim1_use = (mode eq 1) ? -lim1 : lim1
  lim2_use = lim1
endif else begin
  if mode ne 3 and lim2 le lim1 then begin
    print,'ERROR in extend1: Must have lim2 > lim1 for mode = 1 or 2'
    return
  endif else if mode eq 3 and (lim1 lt 0 or lim2 lt 0) then begin
    print,'ERROR in extend1: Must have lim1, lim2 >= 0 for mode = 3'
    return
  endif
  lim1_use = lim1
  lim2_use = lim2
endelse

; Now process the spectrum (or spectra)

for n=0L,nuse-1 do begin

  ; Get a few elements of this pair

  useStruc = (keyword_set(st)) ? *stack1[n] : *loaded1
  spec = useStruc.spec
  tags1 = useStruc.tags
  dfreq = tags1.dfreq
  posfr = tags1.posfr
  xjcen = tags1.xjcen
  jsnr1 = tags1.jsnr1
  jsnr2 = tags1.jsnr2
  ndata = useStruc.ndata
  f = (keyword_set(st)) ? posfr*dfreq*(findgen(ndata) - xjcen) : useStruc.freq
  
  ; Assign the left and right limits according to the
  ; method specified via the mode parameter

  if mode le 1 then begin

    ; Convert frequency limits to frequency bins

    if posfr eq 1 then begin
      bin1 = xjcen + round(lim1_use/dfreq)
      bin2 = xjcen + round(lim2_use/dfreq)
    endif else begin
      bin1 = xjcen - round(lim2_use/dfreq)
      bin2 = xjcen - round(lim1_use/dfreq)
    endelse
  endif else if mode eq 2 then begin

    ; Directly specify frequency bins

    bin1 = round(1.0*lim1_use)
    bin2 = round(1.0*lim2_use)
  endif else begin

    ; Extend frequency bins according to specified margins

    bin1 = -round(1.0*lim1_use)
    bin2 = ndata - 1 + round(1.0*lim2_use)
  endelse

  ; Make sure that the new limits don't fall within the spectrum
  ; and that more than two frequency bins remain after extending

  if bin1 gt 0 or bin2 lt (ndata-1) then begin
    print,'ERROR on spectrum #',n+1,': Requested frequency bins [', $
          bin1,'-',bin2,'] omit part of the spectrum [0-',ndata-1,']', $
          format='(4(a,i0),a)'
    print,'       Pair left unchanged'
  endif else if bin2 le bin1+1 then begin
    print,'ERROR on spectrum #',n+1,': Requested frequency bins [', $
          bin1,'-',bin2,'] would leave too few points in extended spectrum', $
          format='(3(a,i0),a)'
    print,'       Spectrum left unchanged'
  endif else begin
    newndata = bin2 - bin1 + 1
    newxjcen = xjcen - bin1
    newf = posfr*dfreq*(findgen(newndata) - newxjcen)
    tags1.xjcen = newxjcen
    tags1.jsnr1 = (jsnr1 - bin1) > 0L
    tags1.jsnr2 = (jsnr2 - bin1) < (newndata - 1)
    if keyword_set(normal) then begin
      newspec = float(randomn(seed, newndata, /double, /normal))
    endif else begin
      newspec = fltarr(newndata)
    endelse
    newspec[(-bin1):(ndata-1-bin1)] = spec
    print,'Extending spectrum #',n+1,' to ',newndata,' points, frequency bins ', $
          bin1,' - ',bin2,'  (',newf[0]-0.5*posfr*dfreq,' Hz to ', $
          newf[newndata-1]+0.5*posfr*dfreq,' Hz)', $
          format='(a,i0,a,i0,a,i0,a,i0,a,f9.2,a,f9.2,a)'
    if keyword_set(st) then begin
      newStruc = {group:useStruc.group, spec:newspec, tags:tags1, $
                  extratags:useStruc.extratags, ndata:newndata,   $
                  ntags:useStruc.ntags, nextra:useStruc.nextra,   $
                  tname:useStruc.tname, pol:useStruc.pol}
    endif else begin
      newStruc = {freq:newf, spec:newspec, tags:tags1, extratags:useStruc.extratags, $
                  ndata:newndata, ntags:useStruc.ntags, nextra:useStruc.nextra,      $
                  tname:useStruc.tname, pol:useStruc.pol}
    endelse

    if keyword_set(st) then *stack1[n] = newStruc else *loaded1 = newStruc

  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro extendi,mode,doprange=doprange,delrange=delrange,normal=normal,stack=st,help=help

; Extend the loaded image by adding extra rows and/or columns with zero signal
;
; If keyword normal is set, use pseudorandom normal variates
;     (zero mean, unit variance) rather than zeros
;
; If keyword st is set, extend every image in the image stack instead

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then mode = -1 ; just so it's defined

if keyword_set(help) or n_params() ne 1 or mode lt 0 or mode gt 3                   $
        or (mode eq 0 and (n_elements(delrange) gt 0 or n_elements(doprange) gt 0   $
                                                 or not keyword_set(st)  ))         $
        or (mode ne 0 and n_elements(delrange) eq 0 and n_elements(doprange) eq 0)  $
        then begin
  print,' '
  print,'extendi,mode[,doprange=doprange][,delrange=delrange][,/stack][,/help]'
  print,' '
  print,'Extend the loaded image by adding extra rows and/or columns'
  print,'     with zero signal, based on:'
  print,' '
  print,'  mode = 0: max delay-Doppler range covered by the various stacki images'
  print,'            (must set /stack and not use delrange or doprange in this mode)'
  print,' '
  print,'  mode = 1: delay range in usec and/or Doppler range in Hz'
  print,'            (doprange=100 is the same as doprange=[-100,100])'
  print,' '
  print,'  mode = 2: delay rows and/or Doppler columns (0-based)'
  print,' '
  print,'  mode = 3: delay row margins and/or Doppler column margins'
  print,'            (delrange=10 is the same as delrange=[10,10])'
  print,'            (doprange=15 is the same as doprange=[15,15])'
  print,' '
  print,'/normal extends the image with pseudorandom normally distributed'
  print,'     zero-mean, unit-variance deviates rather than with all zeros'
  print,' '
  print,'/stack extends every stacki image instead of the loaded image'
  print,' '
  return
endif

; Decide which images we're extending, check that all necessary extra tags
; are present, and store some information for later use

if not keyword_set(st) then begin

  ; Extend the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,'ERROR in extendi: No image is loaded'
    return
  endif
  nuse = 1L

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
  if isnull(eph_col) or isnull(eph_row) or isnull(dfreq) or isnull(delayunit) then begin
    print,"ERROR in extendi: The loaded image", $
          " needs the 'eph_col' and 'eph_row' and 'fres' tags,",format='(2a)'
    print,"                    plus either the 'delayunit' or 'baudlen' tag"
    return
  endif

endif else begin

  ; Extend every image in the image stack

  if nstacki eq 0 then begin
    print,'ERROR in extendi: The image stack is empty'
    return
  endif
  nuse = nstacki
  width = fltarr(nuse)
  height = fltarr(nuse)
  dfreq = fltarr(nuse)
  delayunit = fltarr(nuse)
  eph_col = fltarr(nuse)
  eph_row = fltarr(nuse)

  for n=0L,nuse-1 do begin
    width[n] = (*stacki[n]).width
    height[n] = (*stacki[n]).height
    dfreq[n] = getextrai('fres',stack=(n+1))
    eph_col[n] = getextrai('eph_col',stack=(n+1))
    eph_row[n] = getextrai('eph_row',stack=(n+1))
    delayunit[n] = getextrai('delayunit',stack=(n+1))
    spb = getextrai('samples_per_baud',stack=(n+1))
    rows_per_baud = getextrai('rows_per_baud',stack=(n+1))
    baud = getextrai('baudlen',stack=(n+1))
    if isnull(delayunit[n]) and notnull(baud) then begin
      if isnull(spb) then spb = 1L
      if isnull(rows_per_baud) then rows_per_baud = spb
      delayunit[n] = baud/rows_per_baud
    endif
    if isnull(eph_col) or isnull(eph_row) or isnull(dfreq) or isnull(delayunit) then begin
      print,'ERROR in extendi: stacki image #',n+1, $
            " needs the 'eph_col' and 'eph_row' and 'fres' tags,",format='(a,i0,a)'
      print,"                    plus either the 'delayunit' or 'baudlen' tag"
      return
    endif
  endfor

endelse

; Set the min and max delay-Doppler values (or image bins)

if mode eq 0 then begin

  ; Find the maximum delay-Doppler range covered by the various stacki images
  
  doplim1 = 9.99e20
  doplim2 = -9.99e20
  dellim1 = 9.99e20
  dellim2 = -9.99e20
  for n=0L,nuse-1 do begin
    fmin = -dfreq[n]*(eph_col[n] + 0.5)                 ; Doppler (Hz)
    fmax = dfreq[n]*(width[n] - eph_col[n] - 0.5)
    dmin = -delayunit[n]*(eph_row[n] + 0.5)             ; delay (usec)
    dmax = delayunit[n]*(height[n] - eph_row[n] - 0.5)
    doplim1 = doplim1 < fmin
    doplim2 = doplim2 > fmax
    dellim1 = dellim1 < dmin
    dellim2 = dellim2 > dmax
  endfor

endif else begin

  ; Get delay and/or Doppler limits from the keywords

  if n_elements(doprange) gt 0 then begin
    if n_elements(doprange) eq 2 then begin
      doplim1 = doprange[0]
      doplim2 = doprange[1]
      if mode ne 3 and doplim2 lt doplim1 then begin
        print,'ERROR in extendi: Must have doprange[1] >= doprange[0] for mode = 1 or 2'
        return
      endif else if mode eq 3 and (doprange[0] lt 0 or doprange[1] lt 0) then begin
        print,'ERROR in extendi: Must have doprange[0], doprange[1] >= 0 for mode = 3'
        return
      endif
    endif else if n_elements(doprange) eq 1 then begin
      if mode eq 1 then begin
        doplim1 = -doprange[0]
        doplim2 = doprange[0]
      endif else if mode eq 3 then begin
        doplim1 = doprange[0]
        doplim2 = doprange[0]
      endif else begin
        print,'ERROR in extendi: Must specify doprange=[lim1,lim2] for mode = 2'
        return
      endelse
      if doplim2 lt 0 then begin
        print,'ERROR in extendi: Must have scalar doprange >= 0'
        return
      endif
    endif else begin
      print,'ERROR in extendi: doprange must be a scalar or a 2-element vector'
      return
    endelse
  endif

  if n_elements(delrange) gt 0 then begin
    if n_elements(delrange) eq 2 then begin
      dellim1 = delrange[0]
      dellim2 = delrange[1]
      if mode ne 3 and dellim2 lt dellim1 then begin
        print,'ERROR in extendi: Must have delrange[1] >= delrange[0] for mode = 1 or 2'
        return
      endif else if mode eq 3 and (delrange[0] lt 0 or delrange[1] lt 0) then begin
        print,'ERROR in extendi: Must have delrange[0], delrange[1] >= 0 for mode = 3'
        return
      endif
    endif else if n_elements(delrange) eq 1 then begin
      if mode eq 3 then begin
        dellim1 = delrange[0]
        dellim2 = delrange[0]
      endif else begin
        print,'ERROR in extendi: Must specify delrange=[lim1,lim2] for mode = 1 or 2'
        return
      endelse
      if dellim2 lt 0 then begin
        print,'ERROR in extendi: Must have scalar delrange >= 0'
        return
      endif
    endif else begin
      print,'ERROR in extendi: delrange must be a 2-element vector'
      return
    endelse
  endif

endelse

; Now process the image(s)

change_delay = mode eq 0 or n_elements(delrange) gt 0
change_doppler = mode eq 0 or n_elements(doprange) gt 0

print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this image

  if keyword_set(st) then begin
    image = (*stacki[n]).image
    dopbin1 = getextrai('dopbin1',stack=(n+1))
    dopbin2 = getextrai('dopbin2',stack=(n+1))
    delbin1 = getextrai('delbin1',stack=(n+1))
    delbin2 = getextrai('delbin2',stack=(n+1))
  endif else begin
    image = (*loadedi).image
    dopbin1 = getextrai('dopbin1')
    dopbin2 = getextrai('dopbin2')
    delbin1 = getextrai('delbin1')
    delbin2 = getextrai('delbin2')
  endelse
  f = dfreq[n]*(findgen(width[n]) - eph_col[n])       ; Doppler (Hz)
  d = delayunit[n]*(findgen(height[n]) - eph_row[n])  ; delay (usec)
  
  ; Assign the left and right limits according to the
  ; method specified via the mode parameter

  if mode le 1 then begin

    ; Convert delay-Doppler limits to image bins

    if change_delay then begin
      deledgebin1 = round(eph_row[n] + dellim1/delayunit[n])
      deledgebin2 = round(eph_row[n] + dellim2/delayunit[n])
    endif else begin
      deledgebin1 = 0L
      deledgebin2 = height[n] - 1
    endelse

    if change_doppler then begin
      dopedgebin1 = round(eph_col[n] + doplim1/dfreq[n])
      dopedgebin2 = round(eph_col[n] + doplim2/dfreq[n])
    endif else begin
      dopedgebin1 = 0L
      dopedgebin2 = width[n] - 1
    endelse

  endif else if mode eq 2 then begin

    ; Directly specify image bins

    if change_delay then begin
      deledgebin1 = round(1.0*dellim1)
      deledgebin2 = round(1.0*dellim2)
    endif else begin
      deledgebin1 = 0L
      deledgebin2 = height[n] - 1
    endelse

    if change_doppler then begin
      dopedgebin1 = round(1.0*doplim1)
      dopedgebin2 = round(1.0*doplim2)
    endif else begin
      dopedgebin1 = 0L
      dopedgebin2 = width[n] - 1
    endelse

  endif else begin

    ; Extend image according to specified margins

    if change_delay then begin
      deledgebin1 = -round(1.0*dellim1)
      deledgebin2 = height[n] - 1 + round(1.0*dellim2)
    endif else begin
      deledgebin1 = 0L
      deledgebin2 = height[n] - 1
    endelse

    if change_doppler then begin
      dopedgebin1 = -round(1.0*doplim1)
      dopedgebin2 = width[n] - 1 + round(1.0*doplim2)
    endif else begin
      dopedgebin1 = 0L
      dopedgebin2 = width[n] - 1
    endelse

  endelse

  ; Make sure that the new limits don't fall within the image and that
  ; more than two bins remain in each dimension after extending

  if dopedgebin1 gt 0 or dopedgebin2 lt (width[n]-1) then begin
    print,'ERROR on image #',n+1,': Requested Doppler bins [', $
          dopedgebin1,'-',dopedgebin2,'] omit part of the image [0-', $
          width[n]-1,']',format='(4(a,i0),a)'
    print,'       Image left unchanged'
  endif else if dopedgebin2 le dopedgebin1+1 then begin
    print,'ERROR on image #',n+1,': Requested Doppler bins [', $
          dopedgebin1,'-',dopedgebin2, $
          '] would leave too few points in extended image', $
          format='(3(a,i0),a)'
    print,'       Image left unchanged'
  endif else if deledgebin1 gt 0 or deledgebin2 lt (height[n]-1) then begin
    print,'ERROR on image #',n+1,': Requested delay bins [', $
          deledgebin1,'-',deledgebin2,'] omit part of the image [0-', $
          height[n]-1,']',format='(4(a,i0),a)'
    print,'       Image left unchanged'
  endif else if deledgebin2 le deledgebin1+1 then begin
    print,'ERROR on image #',n+1,': Requested delay bins [', $
          deledgebin1,'-',deledgebin2, $
          '] would leave too few points in extended image', $
          format='(3(a,i0),a)'
    print,'       Image left unchanged'
  endif else begin

    ; Everything's OK, so change the image and display the results
    ;
    ; Note that you have to replace the entire image structure: IDL won't
    ; let you change the dimensions of the existing image array.

    newwidth = dopedgebin2 - dopedgebin1 + 1
    newheight = deledgebin2 - deledgebin1 + 1
    neweph_col = eph_col[n] - dopedgebin1
    neweph_row = eph_row[n] - deledgebin1
    if keyword_set(normal) then begin
      newimage = float(randomn(seed, newwidth, newheight, /double, /normal))
    endif else begin
      newimage = fltarr(newwidth, newheight)
    endelse
    newimage[(-dopedgebin1):(width[n]-1-dopedgebin1), (-deledgebin1):(height[n]-1-deledgebin1)] = image
    newf = dfreq[n]*(findgen(newwidth) - neweph_col)
    newd = delayunit[n]*(findgen(newheight) - neweph_row)
    if notnull(dopbin1) then dopbin1 = (dopbin1 - dopedgebin1) > 0L
    if notnull(dopbin2) then dopbin2 = (dopbin2 - dopedgebin1) < (newwidth - 1)
    if notnull(delbin1) then delbin1 = (delbin1 - deledgebin1) > 0L
    if notnull(delbin2) then delbin2 = (delbin2 - deledgebin1) < (newheight - 1)

    if keyword_set(st) then begin
      setextrai,'f','eph_col',neweph_col,stack=(n+1)
      setextrai,'f','eph_row',neweph_row,stack=(n+1)
      if notnull(dopbin1) then setextrai,'i','dopbin1',dopbin1,stack=(n+1)
      if notnull(dopbin2) then setextrai,'i','dopbin2',dopbin2,stack=(n+1)
      if notnull(delbin1) then setextrai,'i','delbin1',delbin1,stack=(n+1)
      if notnull(delbin2) then setextrai,'i','delbin2',delbin2,stack=(n+1)
      oldStruc = *stacki[n]
      newStruc = {group:oldStruc.group,                          $
                  image:newimage, extratags:oldStruc.extratags,  $
                  width:newwidth, height:newheight,              $
                  nextra:oldStruc.nextra, tname:oldStruc.tname, pol:oldStruc.pol}
      *stacki[n] = newStruc
    endif else begin
      setextrai,'f','eph_col',neweph_col
      setextrai,'f','eph_row',neweph_row
      if notnull(dopbin1) then setextrai,'i','dopbin1',dopbin1
      if notnull(dopbin2) then setextrai,'i','dopbin2',dopbin2
      if notnull(delbin1) then setextrai,'i','delbin1',delbin1
      if notnull(delbin2) then setextrai,'i','delbin2',delbin2
      oldStruc = *loadedi
      newStruc = {image:newimage, extratags:oldStruc.extratags,  $
                  width:newwidth, height:newheight,              $
                  nextra:oldStruc.nextra, tname:oldStruc.tname, pol:oldStruc.pol}
      *loadedi = newStruc
    endelse

    printstring1 = 'Extending image #' + string(n+1,format='(i0)') + ' to '  $
                   + string(newwidth,format='(i0)') + 'x'                     $
                   + string(newheight,format='(i0)') + ' pixels, '
    nskip1 = strlen(printstring1)
    printstring1 = printstring1 + 'Doppler bins '               $
                   + string(dopedgebin1,format='(i0)') + ' - '  $
                   + string(dopedgebin2,format='(i0)')
    formatstring = '(' + string(nskip1,format='(i0)') + 'x,a)'
    printstring2 = string('delay   bins ',format=formatstring)    $
                   + string(deledgebin1,format='(i0)')  + ' - '   $
                   + string(deledgebin2,format='(i0)')
    nskip2 = strlen(printstring1) - strlen(printstring2)
    if nskip2 gt 0 then begin
      for i=0L,nskip2-1 do printstring2 = printstring2 + ' '
    endif else if nskip2 lt 0 then begin
      for i=0L,-nskip2-1 do printstring1 = printstring1 + ' '
    endif
    printstring1 = printstring1 + '  ('                                             $
                   + string(newf[0]-0.5*dfreq[n],format='(f8.2)') + ' Hz to '       $
                   + string(newf[newwidth-1]+0.5*dfreq[n],format='(f8.2)') + ' Hz)'
    printstring2 = printstring2 + '  ('                                          $
                   + string(newd[0]-0.5*delayunit[n],format='(f8.2)')            $
                   + ' us to '                                                   $
                   + string(newd[newheight-1]+0.5*delayunit[n],format='(f8.2)')  $
                   + ' us)'
    print,printstring1
    print,printstring2

  endelse

endfor

print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro decimate,nbins,sample=sampnum,stack=st,help=help

; Decimate the loaded pair by averaging nbins bins at a time
;
; If keyword st is set, decimate every pair in the stack instead
;
; This routine is based in part on Mike Nolan's decimate.pro.
;
; 2013 Jun 12: fix computation of new sdev tag

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

use_sampling = (n_elements(sampnum) gt 0)
if not use_sampling then sampnum = 0       ; just to define it

if keyword_set(help) or n_params() ne 1 then begin
  print,' '
  print,'decimate,n[,sample=sampnum][,/stack][,/help]'
  print,' '
  print,'Decimate the loaded pair by averaging adjacent groups of n bins'
  print,' '
  print,'Groups are chosen such that the 0-Hz bin (xjcen) is the leftmost'
  print,"bin in its group; 'leftover' bins at either end are ignored"
  print,' '
  print,'If the sample keyword is used, instead replace each group of n bins'
  print,'with the sampnum-th bin in the group (sampnum = 1 or 2 or ... or n)'
  print,' '
  print,'/stack decimates every stack pair instead of the loaded pair'
  print,' '
  return
endif else if use_sampling and (sampnum le 0 or sampnum gt nbins) then begin
  print,' '
  print,'decimate,n[,sample=sampnum][,/stack][,/help]'
  print,' '
  print,'ERROR: must have 1 <= sampnum <= ',nbins,format='(a,i0)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Decimate the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in decimate: No pair is loaded'
    return
  endif
  nuse = 1
endif else begin

  ; Decimate every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in decimate: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; Process the pair(s)

for n=0L,nuse-1 do begin

  ; Get a few elements of this pair

  useStruc = (keyword_set(st)) ? *stack[n] : *loaded
  pair = useStruc.spec
  tags = useStruc.tags
  dfreq = tags[0].dfreq
  posfr = tags[0].posfr
  xjcen = tags[0].xjcen
  jsnr1 = tags.jsnr1
  jsnr2 = tags.jsnr2
  ndata = useStruc.ndata
  f = (keyword_set(st)) ? posfr*dfreq*(findgen(ndata) - xjcen) : useStruc.freq

  ; Compute the leftmost and rightmost bins in the original pair which will
  ; contribute to the decimated pair, leaving off any "leftover" bins on the
  ; edges which aren't divided evenly by nbins; then compute the length of
  ; the decimated pair.

  leftbin = xjcen - nbins*(xjcen/nbins)
  rightbin = xjcen + nbins*((ndata - xjcen)/nbins) - 1
  ndatanew = (rightbin - leftbin + 1)/nbins

  ; Get our starting version of the tags, to be edited later

  tagsnew = tags

  ; Perform decimation and create an appropriate comment string.  If we are averaging,
  ; renormalize the spectra to unit noise standard deviation, and reduce the rmsc tag
  ; to indicate reduced r.m.s. noise.
  ;
  ; If we are sampling every nbins-th spectral bin, we increase the sdev tag by a factor
  ; of nbins in order to keep the cross section about the same as before.  If we are
  ; instead averaging nbins spectral bins at a time, there is an additional noise reduction
  ; factor of sqrt(nbins), so we only increase sdev by a factor of sqrt(nbins).

  if use_sampling then begin
    pairnew = pair[*, leftbin + nbins*lindgen(ndatanew) + sampnum - 1]
    tagsnew.sdev = tags.sdev * nbins
    commentstring = 'decimated (sampled) by a factor of ' + string(nbins,format='(i0)')  $
                    + ' (sampnum = ' + string(sampnum,format='(i0)') + ')'
  endif else begin
    pairnew = total( reform(pair[*,leftbin:rightbin], 2, nbins, ndatanew), 2 ) / sqrt(nbins)
    tagsnew.sdev = tags.sdev * sqrt(nbins)
    tagsnew.rmsc = tags.rmsc / sqrt(nbins)
    commentstring = 'decimated (averaged) by a factor of ' + string(nbins,format='(i0)')
  endelse

  ; Get the new frequency vector and frequency resolution

  fnew = f[leftbin + nbins*lindgen(ndatanew)]
  tagsnew.dfreq = dfreq*nbins

  ; Get the new 0-Hz bin and the new signal limits; for the latter, err on the
  ; side of making the limits wider rather than narrower

  xjcennew = (xjcen - leftbin)/nbins
  tagsnew.jsnr1 = xjcennew - ceil((xjcen - 1.0*jsnr1)/nbins)
  tagsnew.jsnr2 = xjcennew + ceil((jsnr2 - 1.0*xjcen)/nbins)
  tagsnew.xjcen = xjcennew

  ; Replace the old pair with the new one and add the comment as an extra tag

  if keyword_set(st) then begin
    newStruc = {group:useStruc.group, spec:pairnew, tags:tagsnew, $
                extratags:useStruc.extratags, ndata:ndatanew,     $
                ntags:useStruc.ntags, nextra:useStruc.nextra, tname:useStruc.tname}
    *stack[n] = newStruc
    addcomment,commentstring,stack=(n+1)
  endif else begin
    newStruc = {freq:fnew, spec:pairnew, tags:tagsnew, $
                extratags:useStruc.extratags, ndata:ndatanew, ntags:useStruc.ntags, $
                nextra:useStruc.nextra, tname:useStruc.tname}
    *loaded = newStruc
    addcomment,commentstring
  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro decimate1,nbins,sample=sampnum,stack=st,help=help

; Decimate the loaded single-channel spectrum by averaging nbins bins at a time
;
; If keyword st is set, decimate every spectrum in the single-channel stack instead
;
; 2013 Jun 12: fix computation of new sdev tag
;              forgot to include polarization in new structure

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

use_sampling = (n_elements(sampnum) gt 0)
if not use_sampling then sampnum = 0       ; just to define it

if keyword_set(help) or n_params() ne 1 then begin
  print,' '
  print,'decimate1,n[,sample=sampnum][,/stack][,/help]'
  print,' '
  print,'Decimate the loaded single-channel spectrum  by averaging adjacent groups of n bins'
  print,' '
  print,'Groups are chosen such that the 0-Hz bin (xjcen) is the leftmost'
  print,"bin in its group; 'leftover' bins at either end are ignored"
  print,' '
  print,'If the sample keyword is used, instead replace each group of n bins'
  print,'with the sampnum-th bin in the group (sampnum = 1 or 2 or ... or n)'
  print,' '
  print,'/stack decimates every stack1 spectrum instead of the loaded spectrum'
  print,' '
  return
endif else if use_sampling and (sampnum le 0 or sampnum gt nbins) then begin
  print,' '
  print,'decimate1,n[,sample=sampnum][,/stack][,/help]'
  print,' '
  print,'ERROR: must have 1 <= sampnum <= ',nbins,format='(a,i0)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Decimate the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in decimate1: No single-channel spectrum is loaded'
    return
  endif
  nuse = 1
endif else begin

  ; Decimate every spectrum in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in decimate1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
endelse

; Process the spectra

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack1[n] : *loaded1
  spec = useStruc.spec
  tags1 = useStruc.tags
  dfreq = tags1.dfreq
  posfr = tags1.posfr
  xjcen = tags1.xjcen
  jsnr1 = tags1.jsnr1
  jsnr2 = tags1.jsnr2
  ndata = useStruc.ndata
  f = (keyword_set(st)) ? posfr*dfreq*(findgen(ndata) - xjcen) : useStruc.freq

  ; Compute the leftmost and rightmost bins in the original spectrum that will
  ; contribute to the decimated spectrum, leaving off any "leftover" bins on the
  ; edges which aren't divided evenly by nbins; then compute the length of
  ; the decimated spectrum.

  leftbin = xjcen - nbins*(xjcen/nbins)
  rightbin = xjcen + nbins*((ndata - xjcen)/nbins) - 1
  ndatanew = (rightbin - leftbin + 1)/nbins

  ; Get our starting version of the tags, to be edited later

  tags1new = tags1

  ; Perform decimation and create an appropriate comment string.  If we are averaging,
  ; renormalize the spectrum to unit noise standard deviation, and reduce the rmsc tag
  ; to indicate reduced r.m.s. noise.
  ;
  ; If we are sampling every nbins-th spectral bin, we increase the sdev tag by a factor
  ; of nbins in order to keep the cross section about the same as before.  If we are
  ; instead averaging nbins spectral bins at a time, there is an additional noise reduction
  ; factor of sqrt(nbins), so we only increase sdev by a factor of sqrt(nbins).

  if use_sampling then begin
    specnew = spec[leftbin + nbins*lindgen(ndatanew) + sampnum - 1]
    tags1new.sdev = tags1.sdev * nbins
    commentstring = 'decimated (sampled) by a factor of ' + string(nbins,format='(i0)')  $
                    + ' (sampnum = ' + string(sampnum,format='(i0)') + ')'
  endif else begin
    specnew = total( reform(spec[leftbin:rightbin], nbins, ndatanew), 1 ) / sqrt(nbins)
    tags1new.sdev = tags1.sdev * sqrt(nbins)
    tags1new.rmsc = tags1.rmsc / sqrt(nbins)
    commentstring = 'decimated (averaged) by a factor of ' + string(nbins,format='(i0)')
  endelse

  ; Get the new frequency vector and frequency resolution

  fnew = f[leftbin + nbins*lindgen(ndatanew)]
  tags1new.dfreq = dfreq*nbins

  ; Get the new 0-Hz bin and the new signal limits; for the latter, err on the
  ; side of making the limits wider rather than narrower
  
  xjcennew = (xjcen - leftbin)/nbins
  tags1new.jsnr1 = xjcennew - ceil((xjcen - 1.0*jsnr1)/nbins)
  tags1new.jsnr2 = xjcennew + ceil((jsnr2 - 1.0*xjcen)/nbins)
  tags1new.xjcen = xjcennew

  ; Replace the old spectrum with the new one and add the comment as an extra tag

  if keyword_set(st) then begin
    newStruc = {group:useStruc.group, spec:specnew, tags:tags1new, $
                extratags:useStruc.extratags, ndata:ndatanew,      $
                ntags:useStruc.ntags, nextra:useStruc.nextra,      $
                tname:useStruc.tname, pol:useStruc.pol}
    *stack1[n] = newStruc
    addcomment1,commentstring,stack=(n+1)
  endif else begin
    newStruc = {freq:fnew, spec:specnew, tags:tags1new, $
                extratags:useStruc.extratags, ndata:ndatanew, ntags:useStruc.ntags, $
                nextra:useStruc.nextra, tname:useStruc.tname, pol:useStruc.pol}
    *loaded1 = newStruc
    addcomment1,commentstring
  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro decimatei,delpix=nrows,doppix=ncols,rowsample=rowsampnum,colsample=colsampnum, $
              stack=st,help=help

; Decimate the loaded image by averaging nrows delay rows and/or ncols Doppler columns at a time
;
; If keyword st is set, decimate every image in the image stack instead
;
; 2013 Jun 12: fix computation of new sdev extra tag
;              switch from a single "sampnum" keyword to "rowsampnum" and "colsampnum" keywords

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

decimate_delay = (n_elements(nrows) gt 0)
decimate_doppler = (n_elements(ncols) gt 0)
use_delay_sampling = (n_elements(rowsampnum) gt 0)
use_doppler_sampling = (n_elements(colsampnum) gt 0)
if not decimate_delay then nrows = 0       ; just to define it
if not decimate_doppler then ncols = 0
if not use_delay_sampling then rowsampnum = 0
if not use_doppler_sampling then colsampnum = 0

if keyword_set(help) or n_params() ne 0 then begin
  print,' '
  print,'decimatei,delpix=nrows,doppix=ncols                                        $'
  print,'          [,rowsample=rowsampnum][,colsample=colsampnum][,/stack][,/help]'
  print,' '
  print,'Decimate the loaded image by averaging adjacent groups of'
  print,'nrows delay rows, ncols Doppler columns, or both'
  print,' '
  print,'If the rowsample keyword is used, instead of averaging rows, replace each'
  print,'group of nrows delay rows with the rowsampnum-th row in the group'
  print,'(rowsampnum = 1 or 2 or ... or nrows)'
  print,' '
  print,'If the colsample keyword is used, instead of averaging columns, replace each'
  print,'group of ncols Doppler columns with the colsampnum-th column in the group'
  print,'(colsampnum = 1 or 2 or ... or ncols)'
  print,' '
  print,'Groups of rows are chosen such that the 0-usec row (eph_row) is the'
  print,"bottommost row in its group (rowsampnum = 1); 'leftover' rows at the"
  print,'top and bottom edges of the image are ignored.'
  print,' '
  print,'Groups of columns are chosen such that the 0-Hz row (eph_col) is the'
  print,"leftmost row in its group (colsampnum = 1); 'leftover' columns at the"
  print,"left and right edges of the image are ignored."
  print,' '
  print,'/stack decimates every stacki image instead of the loaded image'
  print,' '
  return
endif else if not (decimate_delay or decimate_doppler) then begin
  print,' '
  print,'decimatei,delpix=nrows,doppix=ncols                                        $'
  print,'          [,rowsample=rowsampnum][,colsample=colsampnum][,/stack][,/help]'
  print,' '
  print,'ERROR: must specify nrows, ncols, or both'
  print,' '
  return
endif else if use_delay_sampling and $
              (rowsampnum le 0 or (decimate_delay and rowsampnum gt nrows)) then begin
  print,' '
  print,'decimatei,delpix=nrows,doppix=ncols                                        $'
  print,'          [,rowsample=rowsampnum][,colsample=colsampnum][,/stack][,/help]'
  print,' '
  print,'ERROR: must have 1 <= rowsampnum <= ',nrows,format='(a,i0)'
  print,' '
  return
endif else if use_doppler_sampling and $
              (colsampnum le 0 or (decimate_doppler and colsampnum gt ncols)) then begin
  print,' '
  print,'decimatei,delpix=nrows,doppix=ncols                                        $'
  print,'          [,rowsample=rowsampnum][,colsample=colsampnum][,/stack][,/help]'
  print,' '
  print,'ERROR: must have 1 <= colsampnum <= ',ncols,format='(a,i0)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Decimate the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,'ERROR in decimatei: No image is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Decimate every image in the image stack

  if nstacki eq 0 then begin
    print,'ERROR in decimatei: The image stack is empty'
    return
  endif
  nuse = nstacki
endelse

; Process the images

for n=0L,nuse-1 do begin

  ; Get a few elements of this image

  if keyword_set(st) then begin
    useStruc = *stacki[n]
    image = useStruc.image
    width = useStruc.width
    height = useStruc.height
    dfreq = getextrai('fres',stack=(n+1))
    eph_col = getextrai('eph_col',stack=(n+1))
    eph_row = getextrai('eph_row',stack=(n+1))
    delayunit = getextrai('delayunit',stack=(n+1))
    rangeunit = getextrai('rangeunit',stack=(n+1))
    spb = getextrai('samples_per_baud',stack=(n+1))
    rowsPerBaud = getextrai('rows_per_baud',stack=(n+1))
    baud = getextrai('baudlen',stack=(n+1))
    codemethod = getextrai('codemethod',stack=(n+1))
    dopbin1 = getextrai('dopbin1',stack=(n+1))
    dopbin2 = getextrai('dopbin2',stack=(n+1))
    delbin1 = getextrai('delbin1',stack=(n+1))
    delbin2 = getextrai('delbin2',stack=(n+1))
    sdev = getextrai('sdev',stack=(n+1))
    rows_decimated = getextrai('rows_decimated',stack=(n+1))
    cols_decimated = getextrai('cols_decimated',stack=(n+1))
  endif else begin
    useStruc = *loadedi
    image = useStruc.image
    width = useStruc.width
    height = useStruc.height
    dfreq = getextrai('fres')
    eph_col = getextrai('eph_col')
    eph_row = getextrai('eph_row')
    delayunit = getextrai('delayunit')
    rangeunit = getextrai('rangeunit')
    spb = getextrai('samples_per_baud')
    rowsPerBaud = getextrai('rows_per_baud')
    baud = getextrai('baudlen')
    codemethod = getextrai('codemethod')
    dopbin1 = getextrai('dopbin1')
    dopbin2 = getextrai('dopbin2')
    delbin1 = getextrai('delbin1')
    delbin2 = getextrai('delbin2')
    sdev = getextrai('sdev')
    rows_decimated = getextrai('rows_decimated')
    cols_decimated = getextrai('cols_decimated')
  endelse
  if isnull(spb) then spb = 1L
  if isnull(rowsPerBaud) then rowsPerBaud = spb
  stride = spb/rowsPerBaud
  if isnull(delayunit) and notnull(baud) then delayunit = baud/rowsPerBaud
  if isnull(dfreq) or isnull(delayunit) or isnull(eph_col) or isnull(eph_row) then begin
    print,"ERROR in decimatei: Need the 'fres' and 'eph_col' and 'eph_row' tags plus either"
    print,"                    the 'delayunit' or 'baudlen' tag for image #",n+1,format='(a,i0)'
    return
  endif
  if notnull(sdev) and isnull(codemethod) then begin
    print,"WARNING in decimatei: Assuming codemethod = 'short' for image #",n+1,format='(a,i0)'
    codemethod = 'short'
  endif

  ; Decimate in delay if specified

  if decimate_delay then begin

    ; Compute the bottommost and topmost rows in the original image that will
    ; contribute to the decimated image, leaving off any "leftover" rows on the
    ; edges that aren't divided evenly by nrows; then compute the height of the
    ; decimated image.  The groups are chosen such that the row that contains
    ; eph_row is the bottommmost row within its group.

    bottomrow = floor( eph_row - nrows*floor(eph_row/nrows) )
    toprow = bottomrow + nrows*( (height - bottomrow)/nrows ) - 1
    heightnew = (toprow - bottomrow + 1)/nrows

    ; Perform decimation and create an appropriate comment string.  Also
    ; compute bottomoffset, the (floating-point) row number in the old image
    ; that corresponds to the center of the bottommost row in the new image.
    ;
    ; If we are sampling every nrows-th row, we increase the sdev extra tag by a factor
    ; of nrows in order to keep the cross section about the same as before.  If we are
    ; instead averaging nrows rows at a time, there is an additional noise reduction
    ; factor -- a factor less than sqrt(nrows) due to correlated noise in adjacent
    ; image rows -- so we increase sdev by a smaller factor than in the sampling case.

    if use_delay_sampling then begin
      image2 = image[*, bottomrow + nrows*lindgen(heightnew) + rowsampnum - 1]
      delcommentstring = 'decimated (sampled) in delay by a factor of ' $
                         + string(nrows,format='(i0)') + ' (rowsampnum = '    $
                         + string(rowsampnum,format='(i0)') + ')'
      bottomoffset = bottomrow + rowsampnum - 1.0
      if notnull(sdev) then sdev2 = sdev * nrows
    endif else begin

      ; Must account for correlations between adjacent delay rows:
      ; For an explanation of the following steps, see comments
      ; in routine xseci in operations.pro

      n_covar = (codemethod eq 'long_orig') ? spb : 2*spb - 1
      j = indgen(n_covar)
      if codemethod eq 'long_orig' then begin
        covar = ( 1 - j/(1.0*spb) )^2
      endif else begin
        indexmask = intarr(n_covar)
        indexmask[0:spb-1] = 1
        covar = ( ( (2*spb-j-1)*(2*spb-j)*(2*spb-j+1)              $
                     - indexmask*4*(spb-j-1)*(spb-j)*(spb-j+1) )   $
                  / (6.0*spb^3) )^2
      endelse
      corr = covar/covar[0]
      w = where(j mod stride eq 0, rowcount)
      corr_row = corr[w]
      delfactor_row = (nrows - indgen(rowcount)) > 0
      corrsum = 2*total(delfactor_row*corr_row) - nrows   ; includes k < 0
      sdevfactor = sqrt(corrsum)
      image2 = total( reform(image[*, bottomrow:toprow], width, nrows, heightnew), 2 ) $
               / sdevfactor
      if notnull(sdev) then sdev2 = sdev*sdevfactor
      delcommentstring = 'decimated (averaged) in delay by a factor of ' $
                         + string(nrows,format='(i0)')
      bottomoffset = bottomrow + (nrows - 1)/2.0
    endelse

    ; Get the new delay resolution

    delayunitnew = delayunit*nrows
    if notnull(rangeunit) then rangeunitnew = rangeunit*nrows

    ; Get the new eph_row tag and the new signal limits; for the latter, err on the
    ; side of making the limits wider rather than narrower
  
    eph_rownew = (eph_row - bottomoffset) / nrows
    if notnull(delbin1) then delbin1new = floor((delbin1 - bottomoffset)/nrows)
    if notnull(delbin2) then delbin2new = ceil((delbin2 - bottomoffset)/nrows)

    ; Compute the overall delay decimation factor

    rows_decimatednew = (isnull(rows_decimated)) ? nrows : rows_decimated*nrows

  endif else begin
    image2 = image
    heightnew = height
    if notnull(sdev) then sdev2 = sdev
    delayunitnew = delayunit
    if notnull(rangeunit) then rangeunitnew = rangeunit
    eph_rownew = eph_row
    if notnull(delbin1) then delbin1new = delbin1
    if notnull(delbin2) then delbin2new = delbin2
  endelse

  ; Decimate in Doppler if specified

  if decimate_doppler then begin

    ; Compute the leftmost and rightmost columns in the original image which will
    ; contribute to the decimated image, leaving off any "leftover" columns on the
    ; edges that aren't divided evenly by ncols; then compute the width of the
    ; decimated image.  The groups are chosen such that the column that contains
    ; eph_col is the leftmost column within its group.

    leftcol = floor( eph_col - ncols*floor(eph_col/ncols) )
    rightcol = leftcol + ncols*( (width - leftcol)/ncols ) - 1
    widthnew = (rightcol - leftcol + 1)/ncols

    ; Perform decimation and create an appropriate comment string.  Also
    ; compute leftoffset, the (floating-point) column number in the old image
    ; that corresponds to the center of the leftmost column in the new image.
    ;
    ; If we are sampling every ncols-th column, we increase the sdev extra tag by a factor
    ; of ncols in order to keep the cross section about the same as before.  If we are
    ; instead averaging ncols columns at a time, there is an additional noise reduction
    ; factor of sqrt(ncols), so we only increase sdev by a factor of sqrt(ncols).

    if use_doppler_sampling then begin
      imagenew = image2[leftcol + ncols*lindgen(widthnew) + colsampnum - 1, *]
      dopcommentstring = 'decimated (sampled) in Doppler by a factor of ' $
                         + string(ncols,format='(i0)') + ' (colsampnum = '    $
                         + string(colsampnum,format='(i0)') + ')'
      leftoffset = leftcol + colsampnum - 1.0
      if notnull(sdev) then sdevnew = sdev2 * ncols
    endif else begin
      imagenew = total( reform(image2[leftcol:rightcol, *], ncols, widthnew, heightnew), 1 ) $
                 / sqrt(ncols)
      dopcommentstring = 'decimated (averaged) in Doppler by a factor of ' $
                         + string(ncols,format='(i0)')
      leftoffset = leftcol + (ncols - 1)/2.0
      if notnull(sdev) then sdevnew = sdev2 * sqrt(ncols)
    endelse

    ; Get the new frequency resolution

    dfreqnew = dfreq*ncols

    ; Get the new eph_col tag and signal limits; for the latter, err on the
    ; side of making the limits wider rather than narrower
  
    eph_colnew = (eph_col - leftoffset)/ncols
    if notnull(dopbin1) then dopbin1new = floor((dopbin1 - leftoffset)/ncols)
    if notnull(dopbin2) then dopbin2new = ceil((dopbin2 - leftoffset)/ncols)

    ; Compute the overall Doppler decimation factor

    cols_decimatednew = (isnull(cols_decimated)) ? ncols : cols_decimated*ncols

  endif else begin
    imagenew = image2
    widthnew = width
    if notnull(sdev) then sdevnew = sdev2
    dfreqnew = dfreq
    eph_colnew = eph_col
    if notnull(dopbin1) then dopbin1new = dopbin1
    if notnull(dopbin2) then dopbin2new = dopbin2
  endelse

  ; Replace the old image with the new one and add the comment(s) as an extra tag(s)

  if keyword_set(st) then begin
    newStruc = {group:useStruc.group,                          $
                image:imagenew, extratags:useStruc.extratags,  $
                width:widthnew, height:heightnew,              $
                nextra:useStruc.nextra, tname:useStruc.tname, pol:useStruc.pol}
    *stacki[n] = newStruc
    setextrai,'f','delayunit',delayunitnew,stack=(n+1)
    if notnull(rangeunit) then setextrai,'f','rangeunit',rangeunitnew,stack=(n+1)
    setextrai,'f','eph_row',eph_rownew,stack=(n+1)
    if notnull(delbin1) then setextrai,'i','delbin1',delbin1new,stack=(n+1)
    if notnull(delbin2) then setextrai,'i','delbin2',delbin2new,stack=(n+1)
    setextrai,'f','fres',dfreqnew,stack=(n+1)
    setextrai,'f','eph_col',eph_colnew,stack=(n+1)
    if notnull(dopbin1) then setextrai,'i','dopbin1',dopbin1new,stack=(n+1)
    if notnull(dopbin2) then setextrai,'f','dopbin2',dopbin2new,stack=(n+1)
    if notnull(sdev) then setextrai,'f','sdev',string(sdevnew,format='(e10.4)'),stack=(n+1)
    if decimate_delay then begin
      setextrai,'i','rows_decimated',rows_decimatednew,stack=(n+1)
      addcommenti,delcommentstring,stack=(n+1)
    endif
    if decimate_doppler then begin
      setextrai,'i','cols_decimated',cols_decimatednew,stack=(n+1)
      addcommenti,dopcommentstring,stack=(n+1)
    endif
  endif else begin
    newStruc = {image:imagenew, extratags:useStruc.extratags,  $
                width:widthnew, height:heightnew,              $
                nextra:useStruc.nextra, tname:useStruc.tname, pol:useStruc.pol}
    *loadedi = newStruc
    setextrai,'f','delayunit',delayunitnew
    if notnull(eph_row) then setextrai,'f','eph_row',eph_rownew
    if notnull(delbin1) then setextrai,'i','delbin1',delbin1new
    if notnull(delbin2) then setextrai,'i','delbin2',delbin2new
    setextrai,'f','fres',dfreqnew
    if notnull(eph_col) then setextrai,'f','eph_col',eph_colnew
    if notnull(dopbin1) then setextrai,'i','dopbin1',dopbin1new
    if notnull(dopbin2) then setextrai,'f','dopbin2',dopbin2new
    if notnull(sdev) then setextrai,'f','sdev',string(sdevnew,format='(e10.4)')
    if decimate_delay then begin
      setextrai,'i','rows_decimated',rows_decimatednew
      addcommenti,delcommentstring
    endif
    if decimate_doppler then begin
      setextrai,'i','cols_decimated',cols_decimatednew
      addcommenti,dopcommentstring
    endif
  endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro totalpower,stack=st,help=help

; Create a total-power (OC + SC) sum from the loaded pair and push it onto
; the single-channel stack.
;
; The spectrum's sdev tag will be updated but all other tags and extra tags
; will be taken from the OC half of the pair.  In particular, the jcp tag
; will have a value of l, indicating an OC spectrum; but a comment will be
; added indicating that it's actually a total-power spectrum.
;
; If /stack is set then do this for each pair in the stack
; rather than for the loaded pair

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'totalpower[,/stack][,/help]'
  print,' '
  print,'Create a total-power (OC + SC) spectral sum from the loaded pair and'
  print,'push it onto the single-channel stack.  This sum will have the correct'
  print,'sdev tag but will keep the OC values for all other tags and extra tags.'
  print,' '
  print,'In particular it will keep the OC value (= 1) for the jcp tag,'
  print,'but a comment will be added noting that it is actually total power.'
  print,' '
  print,'If /stack is set, create a total-power sum for each pair in the stack'
  print,'rather than for the loaded pair.'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Sum polarizations for the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in totalpower: No pair is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Sum polarizations for every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in totalpower: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; Store the loaded pair so we can use that "space" and reload the pair later on

storeloaded = ptr_new(*loaded)

; Loop through the pair(s) to be summed

for n=1L,nuse do begin

  ; If the stack keyword is being used, load stack pair n
  ; so that the loaded pair contains the OC tags we need

  if keyword_set(st) then load,n

  ; Get a few elements of this pair

  pair = (*loaded).spec
  tags = (*loaded).tags

  ; Compute the total-power sdev and spectrum

  oc_sdev = tags[0].sdev
  sc_sdev = tags[1].sdev
  tp_sdev = sqrt(oc_sdev^2 + sc_sdev^2)
  oc = pair[0,*]*oc_sdev
  sc = pair[1,*]*sc_sdev
  tp = oc + sc

  ; Put the normalized total-power spectrum into the OC half of the
  ; loaded pair and the total-power sdev into the OC half of that tag,
  ; then add a comment noting that this is actually a total-power spectrum

  (*loaded).spec[0,*] = tp/tp_sdev
  (*loaded).tags[0].sdev = tp_sdev
  addcomment,'This is a total-power spectrum'

  ; Push the total-power spectrum onto the single-channel stack

  split,chan=1,/silent

endfor

; Reload the pair that was there at the start

*loaded = *storeloaded
ptr_free,storeloaded

; Tell the user what was done

if keyword_set(st) then begin
  print,'Total-power spectra for stack pairs 1 - ',nstack, $
        ' added to the single-channel stack (spectra ',nstack1-nstack+1, $
        ' - ',nstack1,')',format='(a,i0,a,i0,a,i0,a)'
endif else begin
  print,'Total-power spectrum for the loaded pair', $
        ' added to the single-channel stack (spectrum ', $
        nstack1,')',format='(2a,i0,a)'
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro stokes,stack=st,txpol=txpol,rot=rot,help=help

; Create Stokes sectra from the loaded pair and push the new (wide) spectrum
; on to the stack.
;
; The spectrum's sdev tag will be updated but all other tags and extra tags
; will be taken from the OC half of the pair.  In particular, the jcp tag
; will have a value of 1, indicating an OC spectrum; but a comment will be
; added indicating that it's actually a total-power spectrum.
;
; If /stack is set then do this for each pair in the stack
; rather than for the loaded pair

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

nstokes = 11 ; How many we are adding.

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'stokes[,txpol=txpol][,/stack][,/help]'
  print,' '
  print, 'Add Stokes (S1 .. S4, dp, dlp) polarizations to the loaded spectrum.'
  print,'push it onto the single-channel stack.  These channels will update the'
  print,'sdev tag but will keep the OC values for all other tags and extra tags.'
  print,' '
  print,"txpol = 'lcp' or 'rcp' default is Arecibo with txpol=lcp, so OC is rcp."
  print," linear tx not yet implemented"
  print,'It will have tags of 5..8 for jcp.'
  print,' '
  print,'If /stack is set, create a total-power sum for each pair in the stack'
  print,'rather than for the loaded pair.'
  print,' '
  return
endif

if n_elements(txpol) eq 0 then txp = 1 else begin
  if (tolower(txpol))[0] eq 'l' then begin
    txp = 1
  endif else if (tolower(txpol))[0] eq 'r' then begin
    txp = -1
  endif else begin
    print, "txpol must be 'lcp' or 'rcp'"
    return
  endelse ; bad
endelse

if not keyword_set(st) then begin

  ; Sum polarizations for the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in stokes: No pair is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Sum polarizations for every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in stokes: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; if we're doing the stack, Store the loaded pair so we can use that "space" and reload the pair later on

if keyword_set(st) then storeloaded = ptr_new(*loaded)

; Loop through the pair(s) to be summed

for n=1L,nuse do begin

  ; If the stack keyword is being used, load stack pair n
  ; so that the loaded pair contains the OC tags we need

  if keyword_set(st) then load,n

  ; Get a few elements of this pair

  pair = (*loaded).spec
  tags = (*loaded).tags
  nspec =(*loaded).ndata
  npol = n_elements(tags)
  sdev = tags.sdev
  if npol lt 4 then begin
    print, "ERROR in stokes: Need complex input channels to compute Stokes"
    return
  end

  ; Compute the total-power sdev and spectrum

  rcsspec = pair[0:3,*] * (sdev[0:3] # replicate(1, nspec))
  s = fltarr(nstokes,nspec)
  stags = replicate(tags[0], nstokes) ; same as OC except as changed
  stags.sdev = 1.0

; s1 = I
  s[0,*] = rcsspec[0,*] + rcsspec[1,*]
  stags[0].jcp = 5

; s2 = Q
  s[1,*] = 2 * rcsspec[2,*]
  stags[1].jcp = 6

; S3 = U
  s[2,*] = -2 * rcsspec[3,*]
  stags[2].jcp = 7
  
; S4 = V
  s[3,*] = (rcsspec[0,*] - rcsspec[1,*]) * txp
  stags[3].jcp = 8

; mu = SC/OC Note: NaN is a legitimate result
  s[4,*] = rcsspec[1,*] / rcsspec[0,*]
  stags[4].jcp = 9

; m = degree of polarization  DP
  rmsp = sqrt(s[1,*]^2+s[2,*]^2+s[3,*]^2)
  s[5,*] = rmsp / s[0,*]
  stags[5].jcp = 10

; m_l = degree of linear polarization DL
  s[6,*] = sqrt(s[1,*]^2+s[2,*]^2) / s[0,*]
  stags[6].jcp = 11

; chi = 1/2 asin (s4 / m s1) CH
  s[7,*] = 0.5 * asin(s[3,*] / rmsp)
  stags[7].jcp = 12

; red = (m s1 ( 1 + sin 2chi)/2) ^{1/2} RD
;     = (rms/2 (1 + S4 / rms)) ^{1/s}
  s[8,*] = sqrt(rmsp * (1 + s[3,*]/rmsp) * 0.5)
  stags[8].jcp = 13

; green = (s1 (1 -m) ) ^{1/2} GN
  s[9,*] = sqrt(s[0,*] - rmsp)
  stags[9].jcp = 14

; blue = (m s1 (1 - 2 chi) /2) ^{1/2} BL
;      = (rms/2 (1 - S4 / rms)) ^{1/s}

  s[10,*] = sqrt(rmsp * (1 - s[3,*]/rmsp) * 0.5)
  stags[10].jcp = 15

; Create a new structure that may (probably is) be a different size from the old one.

  newStruct = {freq:(*loaded).freq, spec:[(*loaded).spec[0:3,*], s], $
            tags:[(*loaded).tags[0:3],stags], extratags:(*loaded).extratags, $
            ndata:(*loaded).ndata, ntags:(*loaded).ntags, $
            nextra:(*loaded).nextra, tname:(*loaded).tname}

  *loaded = newStruct
  if keyword_set(st) then begin
    ptr_free,st[n-1]
    *st[n-1] = *loaded
  endif

endfor

if keyword_set(st) then begin
; Reload the pair that was there at the start

  *loaded = *storeloaded
  ptr_free,storeloaded
endif

; Tell the user what was done

if keyword_set(st) then begin
  print,' spectra for stack pairs 1 - '+strtrim(string(nstack),2), +$
        ' add Stokes channels'
endif else begin
  print,'Loaded pair added Stokes channels'
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro shiftspec,mode,fshift,roundshift=roundshift,full=full,chan=chan,stack=st, $
              silent=silent,help=help

; Shift the loaded pair in frequency
;
; Modified 2009 Sep 8 by CM:
;     Correct for reduced r.m.s. noise that results from interpolation

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common channelBlock, chanstrings, maxchan

if n_params() eq 0 then mode = -1  ; just so it's defined

if keyword_set(help) or n_params() eq 0 or n_params() gt 2 or $
     mode lt 0 or mode gt 2 or (mode eq 0 and n_params() ne 1) or $
     (mode ne 0 and n_params() ne 2) then begin
  print,' '
  print,'shiftspec,mode[,fshift][,/full][,/roundshift][,chan=1 or 2][,/stack]  $'
  print,'                       [,/silent][,/help]'
  print,' '
  print,'Shift the loaded pair in frequency:'
  print,'    shift OC peak to 0 Hz             (mode = 0)'
  print,'    fshift = frequency shift in Hz    (mode = 1)'
  print,'    fshift = frequency shift in bins  (mode = 2)'
  print,' '
  print,'When mode = 0, the OC peak over the defined signal range'
  print,'    (tags jsnr1-jsnr2) is used, unless the /full keyword is set,'
  print,'    in which case shiftspec uses the OC peak for the full spectrum'
  print,' '
  print,"Interpolation (IDL's 'interpolate' procedure with cubic convolution)"
  print,'    is used; for shifts by an integer number of bins, this is'
  print,"    equivalent to a simple shift (IDL's 'shift' procedure), except"
  print,'    that instead of a circular shift the spectra are zero-filled'
  print,' '
  print,'The signal range (jsnr1 and jsnr2 tags) is shifted as well,'
  print,'    but the 0-Hz bin (xjcen tag) is left unchanged'
  print,' '
  print,'/full is relevant when mode = 0: it tells shiftspec to search for the'
  print,'    peak OC signal over the full spectrum rather than just over the'
  print,'    defined signal range (tags jsnr1-jsnr2)'
  print,' '
  print,'/roundshift rounds off fshift so that spectra are shifted by an'
  print,'    integer number of bins; otherwise no rounding is carried out'
  print,' '
  print,'/stack shifts every stack pair rather than the loaded pair'
  print,' '
  print,'/silent suppresses screen output'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Shift the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in shiftspec: No pair is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Shift every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in shiftspec: The pair stack is empty'
    return
  endif
  nuse = nstack
endelse

; Shift spectra using the method specified via the mode parameter

if not keyword_set(silent) then print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this pair

  useStruc = (keyword_set(st)) ? *stack[n] : *loaded
  pair = useStruc.spec
  tags = useStruc.tags
  df = tags[0].dfreq
  posfr = tags[0].posfr
  xjcen = tags[0].xjcen
  sdev = tags.sdev
  jsnr1 = tags.jsnr1
  jsnr2 = tags.jsnr2
  ndata = useStruc.ndata
  npol = n_elements(tags)

; Determine which polarization channels were specified

  shiftpol = intarr(npol)
  if n_elements(chan) eq 0 then begin
    shiftpol = shiftpol + 1
    chanstring = ''
  endif else if chan ge 1 and chan le npol then begin
    shiftpol[chan-1] = 1
    chanstring = chanstrings[chan-1] + ' '
  endif else begin
    print,'Must use chan = 1 (OC) or 2 (SC) .. npol or omit (all)'
    return
  endelse

  ; Figure out how many (floating-point) bins to shift the spectra

  if mode eq 0 then begin
    if keyword_set(full) then begin
      firstbin = 0L
      lastbin  = ndata - 1
    endif else begin
      firstbin = jsnr1[0]
      lastbin  = jsnr2[0]
    endelse
    OCspec = reform(pair[0,*])
    peaksignal = max(OCspec[firstbin:lastbin], maxbin)
    maxbin = maxbin + firstbin
    binshift = -1.0*(maxbin - xjcen)
  endif else if mode eq 1 then begin
    binshift = posfr*fshift/df
  endif else begin
    binshift = 1.0*fshift
  endelse

  ; Round to an integer number of bins if desired

  if keyword_set(roundshift) then binshift = round(binshift)

  ; Shift the spectra via interpolation; avoid using the "shift" function,
  ; even for shifts by integer numbers of bins, because it produces
  ; circular shifts and I prefer filling in the shifted region with zeroes.
  ; Then correct the spectra (and the "sdev" tags) for the fact that
  ; interpolation by a fractional number of bins reduces the r.m.s. noise.

  coords = findgen(ndata) - binshift

  for ch=1,npol do begin
    if shiftpol[ch-1] then begin
      spec = reform(pair[ch-1,*])
      spec = interpolate(spec, coords, cubic=-1.0, missing=0.0)
      sdevfactor = sdevfactor_cubicconv_cw(binshift, cubic=-1.0)
      sdev = sdev*sdevfactor  ; 0.0 < sdevfactor <= 1.0
      spec = spec/sdevfactor
      jsnr1[ch-1] = jsnr1[ch-1] + floor(binshift)
      jsnr2[ch-1] = jsnr2[ch-1] + ceil(binshift)
      if (jsnr1[ch-1] lt 0 and jsnr2[ch-1] ge 0) $
             or (jsnr1[ch-1] lt ndata and jsnr2[ch-1] ge ndata) then begin
        print,'WARNING for pair #',n+1,' chan #',ch, $
              ': signal range partly extends beyond the spectrum', $
              format='(a,i0,a,i0,a)'
        print,'        -- signal range truncated'
        jsnr1[ch-1] = (jsnr1[ch-1] > 0L) < (ndata - 1)
        jsnr2[ch-1] = (jsnr2[ch-1] > 0L) < (ndata - 1)
      endif else if (jsnr1[ch-1] ge ndata or jsnr2[ch-1] lt 0) then begin
        print,'WARNING for pair #',n+1,' chan #',ch, $
              ': signal range is entirely outside the spectrum', $
              format='(a,i0,a,i0,a)'
        print,'        -- signal range reset to the full spectrum'
        jsnr1[ch-1] = 0L
        jsnr2[ch-1] = ndata - 1
      endif
      if keyword_set(st) then begin
        (*stack[n]).spec[ch-1,*] = spec
        (*stack[n]).tags[ch-1].sdev = sdev[ch-1]
        (*stack[n]).tags[ch-1].jsnr1 = jsnr1[ch-1]
        (*stack[n]).tags[ch-1].jsnr2 = jsnr2[ch-1]
      endif else begin
        (*loaded).spec[ch-1,*] = spec
        (*loaded).tags[ch-1].sdev = sdev[ch-1]
        (*loaded).tags[ch-1].jsnr1 = jsnr1[ch-1]
        (*loaded).tags[ch-1].jsnr2 = jsnr2[ch-1]
      endelse
    endif
  endfor

  ; Display the results

  if keyword_set(st) then begin
    startstring = 'Stack #' + string(n+1, format='(i0)') + ': Shifting '
  endif else begin
    startstring = 'Shifting '
  endelse
  if keyword_set(roundshift) then begin
    shiftstring = 'by ' + string(binshift,format='(i+0)') + ' bins  ('  $
                  + string(posfr*df*binshift,format='(f+0)') + ' Hz)'
  endif else begin
    shiftstring = 'by ' + string(binshift,format='(f+0)') + ' bins  ('  $
                  + string(posfr*df*binshift,format='(f+0)') + ' Hz)'
  endelse
  if not keyword_set(silent) then print,startstring,chanstring,shiftstring

endfor

if not keyword_set(silent) then print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro shiftspec1,mode,fshift,roundshift=roundshift,stack=st, $
               silent=silent,help=help

; Shift the loaded single-channel spectrum in frequency
;
; Modified 2009 Sep 8 by CM:
;     Correct for reduced r.m.s. noise that results from interpolation

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() eq 0 then mode = -1  ; just so it's defined

if keyword_set(help) or n_params() eq 0 or n_params() gt 2 or $
     mode lt 0 or mode gt 2 or (mode eq 0 and n_params() ne 1) or $
     (mode ne 0 and n_params() ne 2) then begin
  print,' '
  print,'shiftspec1,mode[,fshift][,/full][,/roundshift][,/stack]  $'
  print,'                        [,/silent][,/help]'
  print,' '
  print,'Shift the loaded single-channel spectrum in frequency:'
  print,'    shift peak to 0 Hz                (mode = 0)'
  print,'    fshift = frequency shift in Hz    (mode = 1)'
  print,'    fshift = frequency shift in bins  (mode = 2)'
  print,' '
  print,'When mode = 0, the spectral peak over the defined signal range'
  print,'    (tags jsnr1-jsnr2) is used, unless the /full keyword is set,'
  print,'    in which case shiftspec1 uses the peak for the full spectrum'
  print,' '
  print,"Interpolation (IDL's 'interpolate' procedure with cubic convolution)"
  print,'    is used; for shifts by an integer number of bins, this is'
  print,"    equivalent to a simple shift (IDL's 'shift' procedure), except"
  print,'    that instead of a circular shift the spectra are zero-filled'
  print,' '
  print,'The signal range (jsnr1 and jsnr2 tags) is shifted as well,'
  print,'    but the 0-Hz bin (xjcen tag) is left unchanged'
  print,' '
  print,'/full is relevant when mode = 0: it tells shiftspec1 to search for'
  print,'    the peak signal over the full spectrum rather than just over the'
  print,'    defined signal range (tags jsnr1-jsnr2)'
  print,' '
  print,'/roundshift rounds off fshift so that spectra are shifted by an'
  print,'    integer number of bins; otherwise no rounding is carried out'
  print,' '
  print,'/stack shifts every spectrum in the single-channel stack'
  print,'    rather than the loaded single-channel spectrum'
  print,' '
  print,'/silent suppresses screen output'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Shift the loaded single-channel spectrum

  if (*loaded1).ndata le 2 then begin
    print,'ERROR in shiftspec1: No single-channel spectrum is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Shift every spectrum in the single-channel stack

  if nstack1 eq 0 then begin
    print,'ERROR in shiftspec1: The single-channel stack is empty'
    return
  endif
  nuse = nstack1
endelse

; Shift spectra using the method specified via the mode parameter

if not keyword_set(silent) then print,' '

for n=0L,nuse-1 do begin

  ; Get a few elements of this spectrum

  useStruc = (keyword_set(st)) ? *stack1[n] : *loaded1
  spec = useStruc.spec
  tags1 = useStruc.tags
  df = tags1.dfreq
  posfr = tags1.posfr
  xjcen = tags1.xjcen
  sdev = tags1.sdev
  jsnr1 = tags1.jsnr1
  jsnr2 = tags1.jsnr2
  ndata = useStruc.ndata

  ; Figure out how many (floating-point) bins to shift the spectra

  if mode eq 0 then begin
    if keyword_set(full) then begin
      firstbin = 0L
      lastbin  = ndata - 1
    endif else begin
      firstbin = jsnr1
      lastbin  = jsnr2
    endelse
    peaksignal = max(spec[firstbin:lastbin], maxbin)
    maxbin = maxbin + firstbin
    binshift = -1.0*(maxbin - xjcen)
  endif else if mode eq 1 then begin
    binshift = posfr*fshift/df
  endif else begin
    binshift = 1.0*fshift
  endelse

  ; Round to an integer number of bins if desired

  if keyword_set(roundshift) then binshift = round(binshift)

  ; Shift the spectra via interpolation; avoid using the "shift" function,
  ; even for shifts by integer numbers of bins, because it produces
  ; circular shifts and I prefer filling in the shifted region with zeroes.
  ; Then correct the spectra (and the "sdev" tags) for the fact that
  ; interpolation by a fractional number of bins reduces the r.m.s. noise.

  coords = findgen(ndata) - binshift

  spec = interpolate(spec, coords, cubic=-1.0, missing=0.0)
  sdevfactor = sdevfactor_cubicconv_cw(binshift, cubic=-1.0)
  sdev = sdev*sdevfactor  ; 0.0 < sdevfactor <= 1.0
  spec = spec/sdevfactor
  jsnr1 = jsnr1 + floor(binshift)
  jsnr2 = jsnr2 + ceil(binshift)
  if (jsnr1 lt 0 and jsnr2 ge 0) $
         or (jsnr1 lt ndata and jsnr2 ge ndata) then begin
    print,'WARNING for spectrum #',n+1, $
          ': signal range partly extends beyond the spectrum', $
          format='(a,i0,a)'
    print,'        -- signal range truncated'
    jsnr1 = (jsnr1 > 0L) < (ndata - 1)
    jsnr2 = (jsnr2 > 0L) < (ndata - 1)
  endif else if (jsnr1 ge ndata or jsnr2 lt 0) then begin
    print,'WARNING for spectrum #',n+1, $
          ': signal range is entirely outside the spectrum', $
          format='(a,i0,a)'
    print,'        -- signal range reset to the full spectrum'
    jsnr1 = 0L
    jsnr2 = ndata - 1
  endif
  if keyword_set(st) then begin
    (*stack1[n]).spec = spec
    (*stack1[n]).tags.sdev = sdev
    (*stack1[n]).tags.jsnr1 = jsnr1
    (*stack1[n]).tags.jsnr2 = jsnr2
  endif else begin
    (*loaded1).spec = spec
    (*loaded1).tags.sdev = sdev
    (*loaded1).tags.jsnr1 = jsnr1
    (*loaded1).tags.jsnr2 = jsnr2
  endelse

  ; Display the results

  if keyword_set(st) then begin
    startstring = 'Stack #' + string(n+1, format='(i0)') + ': Shifting '
  endif else begin
    startstring = 'Shifting '
  endelse
  if keyword_set(roundshift) then begin
    shiftstring = 'by ' + string(binshift,format='(i+0)') + ' bins  ('  $
                  + string(posfr*df*binshift,format='(f+0)') + ' Hz)'
  endif else begin
    shiftstring = 'by ' + string(binshift,format='(f+0)') + ' bins  ('  $
                  + string(posfr*df*binshift,format='(f+0)') + ' Hz)'
  endelse
  if not keyword_set(silent) then print,startstring,chanstring,shiftstring

endfor

if not keyword_set(silent) then print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro im2spec,threshold=threshold,stack=st,help=help

; Sum the loaded image in delay over the signal range to produce
; a Doppler-only spectrum, then push this spectrum onto the
; single-channel stack
;
; If keyword st is set, sum every image in the image stack instead

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() ne 0 then begin
  print,' '
  print,'im2spec[,threshold=threshold][,/stack][,/help]'
  print,' '
  print,'Sum the loaded image in delay over the signal range to produce'
  print,'a Doppler-only spectrum, then push this spectrum onto the'
  print,'single-channel stack'
  print,' '
  print,'threshold is the minimum strength (noise standard deviations)'
  print,'    of pixels that are included in the sum; if this keyword'
  print,'    is omitted, all pixels within the signal delay range'
  print,"    (image rows 'delbin1' through 'delbin2') are summed"
  print,' '
  print,'/stack sums every stacki image instead of the loaded image'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Sum the loaded image

  if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
    print,'ERROR in im2spec: No image is loaded'
    return
  endif
  nuse = 1L
endif else begin

  ; Sum every image in the image stack

  if nstacki eq 0 then begin
    print,'ERROR in im2spec: The image stack is empty'
    return
  endif
  nuse = nstacki
endelse

; Store the loaded single-channel spectrum so we can use that "space" and
; reload the spectrum later on

storeloaded1 = ptr_new(*loaded1)

; Process the images

print,' '

for n=0L,nuse-1 do begin

  ; Get the relevant tags and other elements for this image

  if keyword_set(st) then begin
    useStruc = *stacki[n]
    image = useStruc.image
    width = useStruc.width
    height = useStruc.height
    pol = useStruc.pol
    dfreq = getextrai('fres',stack=(n+1))
    eph_col = getextrai('eph_col',stack=(n+1))
    eph_row = getextrai('eph_row',stack=(n+1))
    delayunit = getextrai('delayunit',stack=(n+1))
    rangeunit = getextrai('rangeunit',stack=(n+1))
    spb = getextrai('samples_per_baud',stack=(n+1))
    rowsPerBaud = getextrai('rows_per_baud',stack=(n+1))
    baud = getextrai('baudlen',stack=(n+1))
    codemethod = getextrai('codemethod',stack=(n+1))
    dopbin1 = getextrai('dopbin1',stack=(n+1))
    dopbin2 = getextrai('dopbin2',stack=(n+1))
    delbin1 = getextrai('delbin1',stack=(n+1))
    delbin2 = getextrai('delbin2',stack=(n+1))
    sdev = getextrai('sdev',stack=(n+1))
    jdstart = getextrai('jdstart',stack=(n+1))
    jdmean = getextrai('jdmean',stack=(n+1))
    jdend = getextrai('jdend',stack=(n+1))
    year = getextrai('year',stack=(n+1))
    month = getextrai('month',stack=(n+1))
    day = getextrai('day',stack=(n+1))
    hour = getextrai('hour',stack=(n+1))
    minute = getextrai('minute',stack=(n+1))
    second = getextrai('second',stack=(n+1))
    txoffset = getextrai('txoffset',stack=(n+1))
    int_time = getextrai('int_time',stack=(n+1))
    tx_power = getextrai('tx_power',stack=(n+1))
    nlooks = getextrai('nlooks',stack=(n+1))
    freqs = getextrai('freqs',stack=(n+1))
    elevation = getextrai('elevation',stack=(n+1))
    azimuth = getextrai('azimuth',stack=(n+1))
    rtt = getextrai('rtt',stack=(n+1))
    tsys = getextrai('tsys',stack=(n+1))
    gain_rxtx = getextrai('gain_rxtx',stack=(n+1))
    phase = getextrai('phase',stack=(n+1))
  endif else begin
    useStruc = *loadedi
    image = useStruc.image
    width = useStruc.width
    height = useStruc.height
    pol = useStruc.pol
    dfreq = getextrai('fres')
    eph_col = getextrai('eph_col')
    eph_row = getextrai('eph_row')
    delayunit = getextrai('delayunit')
    rangeunit = getextrai('rangeunit')
    spb = getextrai('samples_per_baud')
    rowsPerBaud = getextrai('rows_per_baud')
    baud = getextrai('baudlen')
    codemethod = getextrai('codemethod')
    dopbin1 = getextrai('dopbin1')
    dopbin2 = getextrai('dopbin2')
    delbin1 = getextrai('delbin1')
    delbin2 = getextrai('delbin2')
    sdev = getextrai('sdev')
    jdstart = getextrai('jdstart')
    jdmean = getextrai('jdmean')
    jdend = getextrai('jdend')
    year = getextrai('year')
    month = getextrai('month')
    day = getextrai('day')
    hour = getextrai('hour')
    minute = getextrai('minute')
    second = getextrai('second')
    txoffset = getextrai('txoffset')
    int_time = getextrai('int_time')
    tx_power = getextrai('tx_power')
    nlooks = getextrai('nlooks')
    freqs = getextrai('freqs')
    elevation = getextrai('elevation')
    azimuth = getextrai('azimuth')
    rtt = getextrai('rtt')
    tsys = getextrai('tsys')
    gain_rxtx = getextrai('gain_rxtx')
    phase = getextrai('phase')
  endelse
  if isnull(spb) then spb = 1L
  if isnull(rowsPerBaud) then rowsPerBaud = spb
  stride = spb/rowsPerBaud
  if isnull(delayunit) and notnull(baud) then delayunit = baud/rowsPerBaud
  if isnull(dfreq) or isnull(delayunit) or isnull(eph_col) or isnull(eph_row) then begin
    print,"ERROR in im2spec: Need the 'fres' and 'eph_col' and 'eph_row' tags plus either"
    print,"                      the 'delayunit' or 'baudlen' tag for image #",n+1,format='(a,i0)'
    return
  endif
  if notnull(sdev) and isnull(codemethod) then begin
    print,"WARNING in im2spec: Assuming codemethod = 'short' for image #",n+1,format='(a,i0)'
    codemethod = 'short'
  endif
  if isnull(delbin1) then begin
    print,"WARNING in im2spec: Setting missing 'delbin1' tag to 0"
    delbin1 = 0L
  endif
  if isnull(delbin2) then begin
    print,"WARNING in im2spec: Setting missing 'delbin2' tag to (nrows-1)"
    delbin2 = height - 1
  endif
  need_jdmean = isnull(jdmean)
  need_ymdhms = (isnull(year) or isnull(month) or isnull(day)         $
                 or isnull(hour) or isnull(minute) or isnull(second))
  if need_jdmean then begin
    if notnull(jdstart) and notnull(jdend) then begin
      jdmean = (jdstart + jdend)/2
    endif else if not need_ymdhms then begin
      print,"WARNING: image #",n+1," is missing 'jdstart' and/or 'jdend'",format='(a,i0,a)'
      print,"         --> will get 'jdmean' from date/time tags ('year' 'month' etc.)"
      jdmean = julday(month, day, year, hour, minute, second)
    endif else begin
      print,"WARNING: image #",n+1," is missing 'jdstart' and/or 'jdend'",format='(a,i0,a)'
      print,"                and is also missing the date/time tags ('year' 'month' etc.)"
      print,"         --> can't compute any date/time tags"
    endelse
  endif

  ; Use the specified signal threshold (if any) to determine how many image rows
  ; contain at least one pixel strong enough to be included in the sum over delay

  image = image[*, delbin1:delbin2]
  if keyword_set(threshold) then begin
    mask = (image ge threshold)
    image = image*mask
    pixels_per_row = round(total(mask, 1))
    rowmax = max(where(pixels_per_row gt 0), min=rowmin)
    if rowmax gt -1 then begin
      nrows = rowmax - rowmin + 1
    endif else begin
      nrows = 0L
      print,'WARNING: image #',n+1,' has no pixels above the specified threshold', $
            format='(a,i0,a)'
      print,'         --> no spectrum will be output'
    endelse
  endif else begin
    nrows = delbin2 - delbin1 + 1
  endelse

  ; Sum over delay, then reduce sdev to indicate reduced r.m.s. noise

  if nrows gt 0 then begin

    ; Must account for correlations between adjacent delay rows:
    ; For an explanation of the following steps, see comments
    ; in routine xseci in operations.pro

    n_covar = (codemethod eq 'long_orig') ? spb : 2*spb - 1
    j = indgen(n_covar)
    if codemethod eq 'long_orig' then begin
      covar = ( 1 - j/(1.0*spb) )^2
    endif else begin
      indexmask = intarr(n_covar)
      indexmask[0:spb-1] = 1
      covar = ( ( (2*spb-j-1)*(2*spb-j)*(2*spb-j+1)              $
                   - indexmask*4*(spb-j-1)*(spb-j)*(spb-j+1) )   $
                / (6.0*spb^3) )^2
    endelse
    corr = covar/covar[0]
    w = where(j mod stride eq 0, rowcount)
    corr_row = corr[w]
    delfactor_row = (nrows - indgen(rowcount)) > 0L
    corrsum = 2*total(delfactor_row*corr_row) - nrows   ; includes k < 0
    sdevfactor = sqrt(corrsum)

    ; Do the delay sum and the sdev reduction

    spec = total(image, 2) / sdevfactor
    if notnull(sdev) then sdevnew = sdev*sdevfactor

    ; Use extra tags to assign values to the (important) CW tags

    tags1 = blanktags1()
    ntags = n_tags(tags1)
    tags1.jcp = (pol eq 'OC') ? 1L : 2L
    tags1.dfreq = dfreq
    tags1.igw = round(delayunit)
    tags1.xjcen = round(eph_col)
    freq = tags1.dfreq*(findgen(width) - tags1.xjcen)
    tags1.jsnr1 = (notnull(dopbin1)) ? round(dopbin1) : 0L
    tags1.jsnr2 = (notnull(dopbin2)) ? round(dopbin2) : (width - 1)
    if notnull(sdev) then tags1.sdev = sdevnew
    if notnull(jdmean) then begin
      caldat_roundsec,jdmean,month,day,year,hour,minute,second
      tags1.iyy = year
      tags1.imm = month
      tags1.idd = day
      tags1.rchour = hour
      tags1.rcmin = minute
      tags1.rcsec = second
    endif
    if notnull(jdstart) then begin
      caldat_roundsec,jdstart,month,day,year,hour,minute,second
      tags1.rcsta = hour*3600.0 + minute*60.0 + second
    endif
    if notnull(jdend) then begin
      caldat_roundsec,jdend,month,day,year,hour,minute,second
      tags1.rcend = hour*3600.0 + minute*60.0 + second
    endif
    if notnull(txoffset) then tags1.doppl = txoffset
    if notnull(int_time) then tags1.tau = int_time
    if notnull(tx_power) then tags1.trpwr = tx_power
    if notnull(nlooks) then tags1.nffts = nlooks
    if notnull(freqs) then tags1.lfft = freqs
    if notnull(elevation) then tags1.elev = elevation
    if notnull(azimuth) then tags1.azim = azimuth
    if notnull(rtt) then tags1.rttim = rtt
    if notnull(tsys) then tags1.tsys = tsys
    if notnull(gain_rxtx) then tags1.gain = gain_rxtx
    if notnull(phase) then tags1.phase = phase

    ; Create a structure for the new single-channel spectrum,
    ; load it, add a comment saying how the spectrum was created,
    ; and push it onto the single-channel stack

    stack1Struc = {freq:freq, spec:spec, tags:tags1, extratags:useStruc.extratags, $
                   ndata:width, ntags:ntags, nextra:useStruc.nextra,               $
                   tname:useStruc.tname, pol:pol}
    *loaded1 = stack1Struc
    addcomment1,'Spectrum created by summing an image over delay'
    push1,/silent
    print,'Image #',n+1,' converted to stack1 spectrum #',nstack1,format='(a,i0,a,i0,a)'

  endif  ; if nrows gt 0

  ; Loop back and process the next image (if any)

endfor

print,' '

; Reload the single-channel spectrum that was there at the start

*loaded1 = *storeloaded1
ptr_free,storeloaded1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
