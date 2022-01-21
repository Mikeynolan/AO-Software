pro cwcross, pair, out, range, zero=zero, smooth=smooth, sc=sc, oc=oc, xrange=xrange, maxf=maxf, siglim=siglim, error=error, id=id, diff=diff, file=file, noset=noset, _extra=_extra

; cwcross loads the requested pair, then plots the integral of the (by default OC) spectrum. It then asks you to click two end points.
; It integrates the signal between the two and computes the total and the (purely statstical) uncertainty. It includes both thermal
; and self-noise. If you specify /sc, it will plot that spectrum and use it for the "mouse" values.

; It also looks for where the first zero-crossing is, starting at the middle of the selected range.

; It is important that this routine be given raw spectra: If you do any smoothing outside, it will not compute the errors correctly. You can 
; give it sums made in the Magri cw reduction software, as that keeps track of looks.
; unless you include the keyword "noset" It also sets siglim, the cross-section and errors in the requested pair

; Common blocks to be able to set the cross-sections
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack
common loadedBlock,loadedi,loaded1,loaded

pol=0

if keyword_set (sc) then begin
	pol=1
endif

self = !dpi / (4 - !dpi) ; fading noise SNR factor from Elachi book.

if (n_params() lt 1) then begin
    print, 'Usage: cwcross, pair, out, range, zero=zero, smooth=smooth, file=file, error=error, /sc, /oc, /siglim, /noset'
    print, 'updates siglim, cross, and crerr for the stack and loaded pair'
    print, 'Note: the loaded pair is overwritten by the unsmoothed spectrum'
    print, 'self-noise is only correct if no smoothing is done outside of this routine'
    print, 'pair is the entry on the stack to load.'
    print, 'If provided, out will contain [bw, zcbw] on return.'
    print, 'If provided, range will contain a list of columns between picked points'
    print, 'zero gives the y value to check for the "zero crossing bandwidth"'
    print, 'smooth is a frequency resolution to smooth the display by (not recommended)'
    print, 'file= will append to the given filename. /file will append to cwcross.out'
    print, 'error= will add the requested systematic relative error in computing the pol ratio (only)'
    print, '/oc or /sc will use that spectrum (oc is default) for plotting and computing the "mouse"'
    print, '/diff will plot the actual spectrum instad of the integral'
    print, 'cross sections. Mouse cross sections are computed using clicked points instead of'
    print, 'computed totals.'
    print, '/siglim  will use the limits from siglim instead of asking to click points'
    print, '/noset will *not* set the siglim, cross-section, and error in the stack'
    return
endif

load, pair

extract, s=s,freq=freq, tags=t, extratags=e, tname=tname
OCarray = double(reform(s[0,*]))
SCarray = double(reform(s[1,*]))
parray=double(reform(s[pol,*]))
looks = t[pol].nffts
sdev = t[*].sdev

; Open output file if requested
if (keyword_set(file)) then begin
 if(size(file, /type) eq 7) then openw,lun, file, /append, /get_lun else openw, lun, 'cwcross.out', /append, /get_lun
 loops=2
endif else loops=1

if (keyword_set(maxf)) then xrange = [-maxf, maxf]
ss = size(xrange)
if (ss[ss[0]+1] eq 0) then xrange = [freq[0], freq[n_elements(freq)-1]]

if (keyword_set(smooth)) then begin
  smoothf, smooth
  extract, s=s,tags=t
  newsdev = t[*].sdev
  parray = double(reform(s[pol,*]))
  load, pair
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

OCsum = total((OCarray[l:r]))
OCsumsq = total((OCarray[l:r])^2)
OCxsec= OCsum * sdev[0]
SCsum = total((SCarray[l:r]))
SCsumsq = total((SCarray[l:r])^2)
SCxsec= SCsum * sdev[1]

OCnoiseg = sqrt(r-l+1.d0) / OCsum
SCnoiseg = sqrt(r-l+1.d0) / SCsum
OCnoises = sqrt(OCsumsq / self^2 / looks) / OCsum
SCnoises = sqrt(SCsumsq / self^2 / looks) / SCsum
OCnoise = sqrt(OCsumsq / self^2 /looks + (r-l+1.d0)) / OCsum
SCnoise = sqrt(SCsumsq / self^2 /looks + (r-l+1.d0)) / SCsum

OCe = OCnoise
OCe = sqrt(OCe*OCe+error*error)*OCxsec
SCe = SCnoise
SCe = sqrt(SCe*SCe+error*error)*SCxsec
ratioerr, SCxsec, SCe, OCxsec, OCe, 0., 1, ratio, ratlo, rathi

if (not keyword_set(siglim)) then begin
; smooth scales, put it back
msum = (y2-y1) * (newsdev[pol] / sdev[pol])
if (pol eq 0) then msumsq = OCsumsq else msumsq = SCsumsq
mnoiseg = sqrt(r-l+1.d0) / msum
mnoises = sqrt(msumsq / self^2 / looks) / msum
mnoise = sqrt(msumsq / self^2 /looks + (r-l+1.d0)) / msum
end

; Loop to print to file and screen
for i = 1, loops do begin
if (i eq 1) then luni = -1 else luni = lun
printf, luni, ' '
printf, luni, format='("Stack:", I3, " Target: ", a16, " Time:", i5, i3.2, i3.2, I3.2, ":", I2.2, ":", I2.2)', pair, tname, t[0].iyy, t[0].imm, t[0].idd, t[0].rchour, t[0].rcmin, t[0].rcsec
printf, luni, "Left: ", string(x1), " Right: ", string(x2), " Bandwidth = " + string(x2 - x1)
printf, luni, strcompress("OC Zero-crossings at " + string(freq[OCleft]) + ' (y=' + string(OCarray[OCleft]) + ") and " + string(freq[OCright]) + " (y=" + string(OCarray[OCright]) + ')')
printf, luni, strcompress("SC Zero-crossings at " + string(freq[SCleft]) + ' (y=' + string(SCarray[SCleft]) + ") and " + string(freq[SCright]) + " (y=" + string(SCarray[SCright]) + ')')
printf, luni, "OC Zero-crossing bandwidth = " + string(freq[OCright] - freq[OCleft])
printf, luni, "SC Zero-crossing bandwidth = " + string(freq[SCright] - freq[SCleft])

printf, luni, "Signal and Noise:"

printf, luni, "          SNR           xsec      thermal       self        total       total"
printf, luni, "         sigmas         km^2      (frac)       (frac)      (frac)      (km^2)"
printf, luni, format='("OC:    ",g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4)', OCsum, OCxsec, OCnoiseg, OCnoises, OCnoise, OCnoise*OCxsec
printf, luni, format='("SC:    ",g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4)', SCsum, SCxsec, SCnoiseg, SCnoises, SCnoise, SCnoise*SCxsec
if ((not keyword_set(siglim)) and not keyword_set(diff)) then begin
printf, luni, format='("mouse",A2,g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4, 2x, g10.4)', (pol eq 1) ? "SC" : "OC", msum, msum*sdev[pol], mnoiseg, mnoises, mnoise, mnoise*msum*sdev[pol]
endif
printf, luni, format='("Looks:", i0, " Chans:", i0, " Pol ratio: ", f6.3, " (", f6.3, "-", f6.3,")")', looks, r-l+1, ratio, ratlo, rathi

endfor

if not (keyword_set(noset)) then begin
  p = pair - 1
  (*loaded).tags[0].cross = OCxsec
  (*loaded).tags[0].crerr = OCnoise*OCxsec
  (*loaded).tags[1].cross = SCxsec
  (*loaded).tags[1].crerr = SCnoise*SCxsec
  (*loaded).tags[0].jsnr1 = l
  (*loaded).tags[1].jsnr1 = l
  (*loaded).tags[0].jsnr2 = r
  (*loaded).tags[1].jsnr2 = r
  (*stack[p]).tags[0].cross = OCxsec
  (*stack[p]).tags[0].crerr = OCnoise*OCxsec
  (*stack[p]).tags[1].cross = SCxsec
  (*stack[p]).tags[1].crerr = SCnoise*SCxsec
  (*stack[p]).tags[0].jsnr1 = l
  (*stack[p]).tags[1].jsnr1 = l
  (*stack[p]).tags[0].jsnr2 = r
  (*stack[p]).tags[1].jsnr2 = r
endif
if (n_params() gt 1) then begin
    out[0] = x2 - x1
    out[1] = freq[OCright] - freq[OCleft]
endif

if (n_params() gt 2) then range = lindgen(r-l+1) + l

if keyword_set(lun) then free_lun, lun

end
