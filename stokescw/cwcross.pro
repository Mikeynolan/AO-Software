pro cwcross, out, range, stack=st, zero=zero, smooth=smooth, chan=chan, xrange=xrange, maxf=maxf, siglim=siglim, error=error, id=id, diff=diff, file=file, noset=noset, help=help,_extra=_extra
;
;
;  cwcross lives in the cmagri idl processing world and uses some
;  low-level interfaces.  It requires integ and f from ~nolan, which have
;  not been integrated into the package.
;  TODO: rename those to get them out of a trivial namespace and integrate
;  into this file.
;
;


; cwcross loads the requested pair, then plots the integral of the (by default OC) spectrum. It then asks you to click two end points.
; It integrates the signal between the two and computes the total and the (purely statstical) uncertainty. It includes both thermal
; and self-noise. If you specify a channel with chan=, it will plot that spectrum and use it for the "mouse" values.

; It also looks for where the first zero-crossing is, starting at the middle of the selected range.

; It is important that this routine be given raw spectra: If you do any smoothing outside, it will not compute the errors correctly. You can 
; give it sums made in the Magri cw reduction software, as that keeps track of looks.
; unless you include the keyword "noset" It also sets siglim, the cross-section and errors in the requested pair

; Common blocks to be able to set the cross-sections
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common loadedBlock,loadedi,loaded1,loaded
common channelBlock, chanstrings, maxchan

self = !dpi / (4 - !dpi) ; fading noise SNR factor from Elachi book.

if (keyword_set(help)) then begin
    print, 'Usage: cwcross, out, range, stack=n, zero=zero, smooth=smooth, file=file, error=error, /sc, /oc, /siglim, /noset'
    print, 'updates siglim, cross, and crerr for the stack entry or loaded pair'
    print, 'self-noise is only correct if no smoothing is done outside of this routine'
    print, 'If provided, out will contain [bw, zcbw] on return.'
    print, 'If provided, range will contain a list of columns between picked points'
    print, 'zero gives the y value to check for the "zero crossing bandwidth"'
    print, 'smooth is a frequency resolution to smooth the display by (not recommended)'
    print, 'file= will append to the given filename. /file will append to cwcross.out'
    print, 'error= will add the requested systematic relative error in computing the pol ratio (only)'
    print, 'chan=n will use that spectrum (oc is default) for plotting and computing the "mouse"'
    print, '/diff will plot the actual spectrum instad of the integral'
    print, 'cross sections. Mouse cross sections are computed using clicked points instead of'
    print, 'computed totals.'
    print, '/siglim  will use the limits from siglim instead of asking to click points'
    print, '/noset will *not* set the siglim, cross-section, and error in the stack'
    return
endif

; We'll be using the loaded pair, save it

saveloaded = *loaded

if keyword_set(st) then begin
  if st ge 1 and st le nstack then load, st else begin
    print, "ERROR in cwcross: stack entry "+strtrim(string(st),2)+" doesn't exist"
    return
  endelse
endif

pol=0
extract, s=s,freq=freq, tags=t, extratags=e, tname=tname
ds = double(s) ; do sums in double precision
npol = n_elements(s[*,0])

if keyword_set (chan) then begin
  if chan ge 1 and chan le npol then pol = chan-1 else begin
    print, "ERROR in cwcross: chan must be 1 .. npol"
    return
  endelse
endif

OCarray = double(reform(ds[0,*]))
SCarray = double(reform(ds[1,*]))
parray=double(reform(s[pol,*]))
looks = t[0].nffts
sdev = t.sdev


; Open output file if requested
if (keyword_set(file)) then begin
 if(size(file, /type) eq 7) then openw,lun, file, /append, /get_lun else openw, lun, 'cwcross.out', /append, /get_lun
 loops=2
endif else loops=1

if (keyword_set(maxf)) then xrange = [-maxf, maxf]
ss = size(xrange)
if (ss[ss[0]+1] eq 0) then xrange = [freq[0], freq[n_elements(freq)-1]]

if (keyword_set(smooth)) then begin
  olds = s
  oldt = t
  smoothf, smooth
  extract, s=s,tags=t
  newsdev = t.sdev
  parray = double(reform(s[pol,*]))
  changespec, s=olds, t=oldt
endif else newsdev = sdev

if (not keyword_set(error)) then error = 0.d0

dims = size(array, /dimensions)

out = dblarr(2)
if (not keyword_set(zero)) then zero = 0.0d0

ipa = integ(parray)
if (keyword_set(diff)) then begin
plot, freq, parray, xrange=xrange, _extra=_extra
endif else begin
plot, freq, ipa, xrange=xrange, _extra=_extra
endelse

if (keyword_set(siglim)) then begin
  l = t[0].jsnr1
  r = t[0].jsnr2
  x1 = freq[l]
  x2 = freq[r]
oplot, [x1,x1],[min(ipa), max(ipa)], linestyle=1
oplot, [x2,x2],[min(ipa), max(ipa)], linestyle=1
print, "Using siglim:", x1, x2
endif else begin

  print, 'click on left point'
  f, x1, y1

  wait, 0.3
  print, 'click on right point'
  f, x2, y2

  if x2 lt x1 then begin
	temp = x1
	x1 = x2
	x2 = temp
  endif

  l = max(where(freq le x1))
  r = max(where(freq le x2))

endelse

fmid = (x1 + x2 ) * 0.5d0

pmid = min (where (freq ge fmid))

OCleft = max(where(OCarray[0:pmid] le zero))
OCright = min(where(OCarray[pmid:*] le zero)) + pmid
SCleft = max(where(SCarray[0:pmid] le zero))
SCright = min(where(SCarray[pmid:*] le zero)) + pmid

sums = total(ds[*,l:r],2)
sumsq = total(ds[*,l:r]^2,2)
xsec = sums * sdev

noiseg = sqrt(r - l + 1.d0) / sums
noises = sqrt(sumsq / self^2 / looks) / sums
noise = sqrt(sumsq / self^2 / looks + (r-l+1.d0)) / sums

e = noise
e = sqrt(e*e+error*error)*xsec
ratioerr, xsec[1], e[1], xsec[0], e[0], 0., 1, ratio, ratlo, rathi

if (not keyword_set(siglim)) then begin
; smooth scales, put it back
msum = (y2-y1) * (newsdev[pol] / sdev[pol])
msumsq = sumsq[pol]
mnoiseg = sqrt(r-l+1.d0) / msum
mnoises = sqrt(msumsq / self^2 / looks) / msum
mnoise = sqrt(msumsq / self^2 /looks + (r-l+1.d0)) / msum
end

; Derived values need to be done by totalling S1..S4

if npol gt 2 then begin
  S1 = sums[0] * sdev[0] + sums[1] * sdev[1]
  S2 = 2 * sums[2] * sdev[2]
  S3 = -2 * sums[3] * sdev[3]
  S4 = sums[0]* sdev[0] - sums[1]*sdev[1]
  mu = sums[1] / sums[0] * sdev[1]/ sdev[0]
  rms2 = sqrt(S2^2 + S3^2)
  rms3 = sqrt(S2^2 + S3^2 + S4^2)
  DP = rms3 / S1
  DLP = rms2 / S1
  chi = 0.5 * asin(S4 / rms3)
  red = sqrt(0.5 * (rms3 + S4))
  green = sqrt(S1 - rms3)
  blue = sqrt(0.5 * (rms3 - S4))
endif
; Loop to print to file and screen
for i = 1, loops do begin
  if (i eq 1) then luni = -1 else luni = lun
  printf, luni, ' '
  printf, luni, format='("Stack:", I3, " Target: ", a16, " Time:", i5, i3.2, i3.2, I3.2, ":", I2.2, ":", I2.2)', keyword_set(st)?st:0, tname, t[0].iyy, t[0].imm, t[0].idd, t[0].rchour, t[0].rcmin, t[0].rcsec
  printf, luni, "Left: ", string(x1), " Right: ", string(x2), " Bandwidth = " + string(x2 - x1)
  printf, luni, strcompress("OC Zero-crossings at " + string(freq[OCleft]) + ' (y=' + string(OCarray[OCleft]) + ") and " + string(freq[OCright]) + " (y=" + string(OCarray[OCright]) + ')')
  printf, luni, strcompress("SC Zero-crossings at " + string(freq[SCleft]) + ' (y=' + string(SCarray[SCleft]) + ") and " + string(freq[SCright]) + " (y=" + string(SCarray[SCright]) + ')')
  printf, luni, "OC Zero-crossing bandwidth = " + string(freq[OCright] - freq[OCleft])
  printf, luni, "SC Zero-crossing bandwidth = " + string(freq[SCright] - freq[SCleft])
  
  printf, luni, "Signal and Noise:"
  
  printf, luni, "          SNR           xsec      thermal       self        total       total"
  printf, luni, "         sigmas         km^2      (frac)       (frac)      (frac)      (km^2)"
  for j = 0, (npol < 4)-1 do begin
      printf, luni, format='(a2,":    ",g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4)', chanstrings[j], sums[j], xsec[j], noiseg[j], noises[j], noise[j], noise[j]*xsec[j]
  endfor
  if npol gt 2 then begin
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "S1", S1
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "S2", S2
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "S3", S3
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "S4", S4
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "MU", mu
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "DP", DP
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "DL", DLP
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "CH", chi
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "RD", red
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "GN", green
     printf, luni, format='(a2,":    ",10x,2x,g10.4)', "BL", blue
  endif

  if ((not keyword_set(siglim)) and not keyword_set(diff)) then begin
  printf, luni, format='("mouse",A2,g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4)', (pol eq 1) ? "SC" : "OC", msum, msum*sdev[pol], mnoiseg, mnoises, mnoise, mnoise*msum*sdev[pol]
  endif
  printf, luni, format='("Looks:", i0, " Chans:", i0, " Pol ratio: ", f6.3, " (", f6.3, "-", f6.3,")")', looks, r-l+1, ratio, ratlo, rathi
  
endfor ; loops

; We did all kinds of stuff, reload
*loaded = saveloaded
if not (keyword_set(noset)) then begin
  if keyword_set(st) then begin
    p = st - 1
    (*stack[p]).tags.cross = xsec
    (*stack[p]).tags.crerr = noise*xsec
    (*stack[p]).tags.jsnr1 = l
    (*stack[p]).tags.jsnr2 = r
  endif else begin
    (*loaded).tags.cross = xsec
    (*loaded).tags.crerr = noise*xsec
    (*loaded).tags.jsnr1 = l
    (*loaded).tags.jsnr2 = r
  endelse
endif
if (n_params() gt 1) then begin
    out[0] = x2 - x1
    out[1] = freq[OCright] - freq[OCleft]
endif

if (n_params() gt 2) then range = lindgen(r-l+1) + l

if keyword_set(lun) then free_lun, lun

end
