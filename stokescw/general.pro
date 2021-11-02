;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function date2doy,month,day,year,help=help

; Turn a calendar date into the day of the year (1-366)

;                Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec
normal_year = [    0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334 ]
leap_year   = [    0,  31,  60,  91, 121, 152, 182, 213, 244, 274, 305, 335 ]

if n_params() ne 3 or keyword_set(help) then begin
  print,' '
  print,'doy = date2doy(month,day,year[,/help])'
  print,' '
  print,'Given the month (1-12), day, and year, return the day of the year (1-366)'
  print,' '
  return, 0
endif else if (year lt 1600 or year gt 2100) then begin
  print,' '
  print,'doy = date2doy(month,day,year[,/help])'
  print,' '
  print,'ERROR: date2doy only handles years 1600 through 2100'
  print,' '
  return, 0
endif

; Modify result depending on whether or not it's a leap year

subCentury = year mod 100
if ( (year mod 400 eq 0) or ( (subCentury ne 0) and (subCentury mod 4 eq 0) ) ) then begin
  doy = leap_year[month-1] + day
endif else begin
  doy = normal_year[month-1] + day
endelse

return, doy

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function is_little_endian,help=help

; Return true if running idl on a little-endian machine,
; false on a big-endian machine

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'Usage:  returnvalue = is_little_endian([/help])'
  print,' '
  print,'returnvalue is true if you are running idl on a little-endian machine,'
  print,'               false if you are running idl on a big-endian machine'
  print,' '
  return, 0
endif

return, byte(1,0)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function isnull,arg,help=help

; Return true if the argument is the null string.
; This function is needed because the simple check "arg eq ''"
; returns true if arg = 0.

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'Usage:  returnvalue = isnull(arg[,/help])'
  print,' '
  print,'returnvalue is true if arg is the null string'
  print,' '
  return, 0
endif

vartype = size(arg, /type)
if vartype eq 7 then begin
  return, (arg eq '')
endif else begin
  return, 0
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function notnull,arg,help=help

; Return true if the argument is a defined variable other than
; the null string.  This function is needed because the simple
; check "arg ne ''" returns false if arg = 0.

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'Usage:  returnvalue = notnull(arg[,/help])'
  print,' '
  print,'returnvalue is true if arg is a defined variable ', $
        'other than the null string',format='(2a)'
  print,' '
  return, 0
endif

vartype = size(arg, /type)
if vartype eq 0 then begin
  return, 0
endif else if vartype eq 7 then begin
  return, (arg ne '')
endif else begin
  return, 1
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function blanktags1,help=help

; Return a blanked CW tag structure for a single-channel spectrum:
;
;     rcsta    0   f   receive start (secs after midnight) (for a sum, refers to earliest run)
;     rcend    1   f   receive stop  (secs after midnight) (for a sum, refers to latest   run)
;    rchour    2   i   hour (0-23) for mean receive time
;     nffts    3   i   number of looks (ffts) incoherently summed in this spectrum
;      elev    4   f   mean telescope receive elevation (deg)
;      azim    5   f   mean telescope receive azimuth   (deg)
;     rttim    6   f   mean round-trip time (s)
;     doppl    7   f   transmit Doppler offset (Hz)
;       idd    8   i   day of month (1-31) for mean receive time
;       imm    9   i   month (1-12) for mean receive time
;       iyy   10   i   four-digit year for mean receive time
;     rcmin   11   i   minute (0-59) for mean receive time
;     rcsec   12   i   second (0-69) for mean receive time
;    rcnsec   13   i   nanosecond (0-999999999) for mean receive time
;     zepch   14   f   obsolete
;     phase   15   f   rotational phase (deg)
;       obs   16   i   observing station code (1 = Arecibo, 2 = Goldstone)
;      itar   17   i   target number
;      irun   18   i   run sequence number
;    jgroup   19   i   obsolete
;      lfft   20   i   fft length
;       igw   21   i   integer-rounded gate width (us)
;     dfreq   22   f   frequency resolution (Hz)
;       tau   23   f   integration time (s)
;      rmsc   24   f   calculated rms calibration noise / (k T_sys B)
;      rmsm   25   f   measured   rms calibration noise / (k T_sys B)
;     xjcen   26   i   0-Hz bin, counting from 0
;     jsnr1   27   i   leftmost  bin containing nonzero echo power, counting from 0
;     jsnr2   28   i   rightmost bin containing nonzero echo power, counting from 0
;      kpts   29   i   obsolete
;       jcp   30   i   polarization (1 = OC, 2 = SC)
;     posfr   31   i   frequency sense (1 = increases rightward, -1 = increases leftward)
;     trpwr   32   f   mean transmit power (kW)
;      tsys   33   f   mean system temperature (K)
;      gain   34   f   mean value of (transmit gain)*(receive gain)/1e12  (dimensionless)
;      sdev   35   f   radar cross-section equivalent of one noise standard deviation (km^2)
;     cross   36   f   calculated radar cross section (km^2)
;     crerr   37   f   error on calculated radar cross section (km^2)
;     nfreq   38   i   number of frequencies per hop cycle
;    frstep   39   f   size of frequency step in hop cycle (Hz)
;     color   40   i   hop color (0 to n-1, or else -1 for a dehopped spectrum)
;     freq1   41   i   frequency for hop color 0 (Hz)
;      util   42   f   utility tag
;
; There are a few differences between these tags and the 43 identically named tags
; used by the "tkplay" package at JPL:
;
; -- xjcen, jsnr1, and jsnr2 are integer here, floating-point in tkplay
; -- xjcen, jsnr1, and jsnr2 are counted from 0 here, from 1 in tkplay
;        (this difference was an accident)
; -- iyy, imm, idd, rchour, rcmin, rcsec, and rcnsec refer to the mean receive epoch here.
;        In tkplay, iyy, imm, and idd refer to the receive-start epoch (along with rcsta);
;        rchour, rcmin, and rcsec are obsolete; and rcnsec refers to the receive-end epoch.
; -- zepch is the epoch (Julian date) of zero rotation phase in tkplay; here, this tag is
;        obsolete (since single-precision isn't sufficient), and the "jd0" extra tag is
;        used to store this zero-phase epoch as a double-precision Julian date
;

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 0 or keyword_set(help) then print,'tags1 = blanktags1([/help])'

tags1 = {rcsta:0.0, rcend:0.0, rchour:0L, nffts:0L,  elev:0.0,     $
         azim:0.0,  rttim:0.0, doppl:0.0, idd:0L,    imm:0L,       $
         iyy:0L,    rcmin:0L,  rcsec:0L,  rcnsec:0L, zepch:2000.0, $
         phase:0.0, obs:1L,    itar:0L,   irun:0L,   jgroup:0L,    $
         lfft:0L,   igw:0L,    dfreq:0.0, tau:0.0,   rmsc:0.0,     $
         rmsm:0.0,  xjcen:0L,  jsnr1:0L,  jsnr2:0L,  kpts:0L,      $
         jcp:0L,    posfr:1L,  trpwr:0.0, tsys:0.0,  gain:0.0,     $
         sdev:1.0,  cross:0.0, crerr:0.0, nfreq:0L,  frstep:0.0,   $
         color:0L,  freq1:0.0, util:0.0 }

return, tags1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function blanktags,npol=npol,help=help

; Return a blanked npol-element vector of tag structures for an OC/SC pair
  
common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 0 or keyword_set(help) then print,'tags1 = blanktags([npol=npol][/help])'
if not keyword_set(npol) then npol=2

return, replicate(blanktags1(), npol)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function stackgroups,help=help

; Return a vector containing all group numbers in the pair stack,
; ordered and listed once each.
;
; If the stack is empty, return the null string.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,'Usage:  vectorvar = stackgroups([/help])'
  return, ''
endif else if nstack eq 0 then begin
  print,'stackgroups: The pair stack is empty'
  return, ''
endif

rawgroups = lindgen(nstack)
for n=0L,nstack-1 do rawgroups[n] = (*stack[n]).group
return, rawgroups[ uniq(rawgroups, sort(rawgroups)) ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function stackgroups1,help=help

; Return a vector containing all group numbers in the
; single-channel stack, ordered and listed once each.
;
; If stack1 is empty, return a null string.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,'Usage:  vectorvar = stackgroups1([/help])'
  return, ''
endif else if nstack1 eq 0 then begin
  print,'stackgroups1: The single-channel stack is empty'
  return, ''
endif

rawgroups = lindgen(nstack1)
for n=0L,nstack1-1 do rawgroups[n] = (*stack1[n]).group
return, rawgroups[ uniq(rawgroups, sort(rawgroups)) ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function stackgroupsi,help=help

; Return a vector containing all group numbers in the
; image stack, ordered and listed once each.
;
; If stacki is empty, return a null string.

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,'Usage:  vectorvar = stackgroupsi([/help])'
  return, ''
endif else if nstacki eq 0 then begin
  print,'stackgroupsi: The image stack is empty'
  return, ''
endif

rawgroups = lindgen(nstacki)
for n=0L,nstacki-1 do rawgroups[n] = (*stacki[n]).group
return, rawgroups[ uniq(rawgroups, sort(rawgroups)) ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getnstack,help=help

; Return the number of pairs in the pair stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,'Usage:  var = getnstack([/help])'
  return, ''
endif

return, nstack
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getnstack1,help=help

; Return the number of spectra in the single-channel stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,'Usage:  var = getnstack1([/help])'
  return, ''
endif

return, nstack1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function getnstacki,help=help

; Return the number of images in the image stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,'Usage:  var = getnstacki([/help])'
  return, ''
endif

return, nstacki
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro caldat_roundsec,juldate,mon,day,year,hr,min,sec,help=help

; Given a Julian date, return the calendar date/time rounded to the nearest second

if n_params() ne 7 or keyword_set(help) then begin
  print,' '
  print,'caldat_roundsec,juldate,mon,day,year,hr,min,sec[,/help]'
  print,' '
  print,'Given a Julian date, return the calendar date/time rounded to the nearest second'
  print,' '
  return
endif

juldate_long = floor(juldate + 0.5D0)
juldate_secs = round(86400.0*(juldate + 0.5D0 - juldate_long))
if juldate_secs eq 86400 then begin
  juldate_long = juldate_long + 1L
  juldate_secs = 0L
endif
caldat,juldate_long,mon,day,year
hr = juldate_secs/3600
juldate_secs = juldate_secs - hr*3600
min = juldate_secs/60
sec = juldate_secs - min*60

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro printcaldat,jd,help=help

; Round Julian date jd to the nearest second,
; then print the calendar date/time equivalent.

monthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'printcaldat,jd[,/help]'
  print,' '
  print,'Round Julian date jd to the nearest second and then print the calendar date/time equivalent'
  print,' '
  return
endif else if size(jd, /type) ne 5 then begin
  print,' '
  print,'printcaldat,jd[,/help]'
  print,' '
  print,'Make sure that Julian date jd is double-precision!'
  print,' '
  return
endif

caldat_roundsec,jd,mon,day,year,hr,min,sec
print,string(year,format='(i4)') + ' ' + monthnames[mon-1] + ' '              $
      + string(day,format='(i2.2)') + ' ' + string(hr,format='(i2.2)') + ':'  $
      + string(min,format='(i2.2)') + ':' + string(sec,format='(i2.2)')

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro pause,prompt=userprompt,reply=userreply

; Give a prompt and then wait for user input (e.g., pressing Enter)
;
; If the userreply keyword is present, return the user input

dummy = ''
if n_elements(userprompt) eq 0 then userprompt = 'Press Enter to continue: '
read,dummy,prompt=userprompt
if arg_present(userreply) then userreply = dummy

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro resettagIO,silent=silent,help=help

; (Re)set parameters for reading and displaying tags

common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline

if keyword_set(help) or n_params() gt 0L then begin
  print,'resettagIO[,/silent],[,/help]'
  return
endif

if not keyword_set(silent) then silent = 0

tags1 = blanktags1()
ntags = n_tags(tags1)

; Get type of each tag:
;   0=undef, 1=byte, 2=int, 3=long, 4=float, 5=double, 6=complex, 7=string,
;   8=structure, 9=complex double, 10=pointer, 11=object reference,
;   12=unsigned int, 13=unsigned long, 14=64-bit int, 15=unsigned 64-bit int
;
; Also define an I/O format corresponding to this type

tagtypes = intarr(ntags)
for n=0L,ntags-1 do begin
  tagtypes[n] = size(tags1.(n), /type)
endfor
tagformats = strarr(ntags)
setformat,/reset,silent=silent

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setformat,width=width,places=places,maxline=maxline,reset=reset, $
              silent=silent,help=help

; Define or list the format for tag output

common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline

default_iowidth = 11
default_iodecplaces = 3
default_iomaxline = 10

if n_params() ne 0 or keyword_set(help) or $
       (keyword_set(reset) and $
        (n_elements(width) gt 0 or n_elements(places) gt 0 or n_elements(maxline) gt 0)) $
       then begin
  print,' '
  print,'setformat[,width=width][,places=places][,maxline=maxline][,/reset]', $
            '[,/silent][,/help]',format='(2a)'
  print,' '
  print,'          width = width of each field, including spaces  (default = ', $
            default_iowidth,')',format='(a,i0,a)'
  print,'          places = # of decimal places                   (default = ', $
            default_iodecplaces,')',format='(a,i0,a)'
  print,'          maxline = maximum # of fields per output line  (default = ', $
            default_iomaxline,')',format='(a,i0,a)'
  print,' '
  print,'          /reset returns to default parameters'
  print,"          Use 'setformat' with no keywords set (or 'showformat') ", $
             "to list the current parameters",format='(2a)'
  print,' '
  return
endif

; Just list the current format parameters if requested

if not (n_elements(width) gt 0 or n_elements(places) gt 0 or n_elements(maxline) gt 0 $
                           or keyword_set(reset) or keyword_set(silent)) then begin
  print,' '
  print,'Tag format: width = ',iowidth,', places = ',iodecplaces,', maxline = ',iomaxline, $
        '  (defaults = ',default_iowidth,', ',default_iodecplaces,', ',default_iomaxline, $
        ')',format='(6(a,i0),a)'
  print,' '
  return
endif

; Get the tag names

tags1 = blanktags1()
ntags = n_tags(tags1)
tagnames = strlowcase(tag_names(tags1))

; (Re)set the formatting parameters and define some useful (pieces of) format strings

if keyword_set(reset) then begin
  iowidth = default_iowidth
  iodecplaces = default_iodecplaces
  iomaxline = default_iomaxline
endif else begin
  if n_elements(width) gt 0 then iowidth = width
  if n_elements(places) gt 0 then iodecplaces = places
  if n_elements(maxline) gt 0 then iomaxline = maxline
endelse

if not keyword_set(silent) then begin
  print,' '
  print,'Tag format: width = ',iowidth,', places = ',iodecplaces,', maxline = ',iomaxline, $
        '  (defaults = ',default_iowidth,', ',default_iodecplaces,', ',default_iomaxline, $
        ')',format='(6(a,i0),a)'
  print,' '
endif

endstring1 = string(iowidth,format='(i0)') + ')'
endstring2 = string(iowidth,format='(i0)') + '.' + string(iodecplaces,format='(i0)') + ')'

; Set the format for the tag names and the individual tags

nameformat = '(a' + endstring1

for n=0L,ntags-1 do begin

  ; Set I/O format based on the variable type of this tag:
  ;   0=undef, 1=byte, 2=int, 3=long, 4=float, 5=double, 6=complex, 7=string,
  ;   8=structure, 9=complex double, 10=pointer, 11=object reference,
  ;   12=unsigned int, 13=unsigned long, 14=64-bit int, 15=unsigned 64-bit int

  case tagtypes[n] of
    0:    tagformats[n] = ''
    1:    tagformats[n] = '(i' + endstring1
    2:    tagformats[n] = '(i' + endstring1
    3:    tagformats[n] = '(i' + endstring1
    4:    tagformats[n] = '(f' + endstring2
    5:    tagformats[n] = '(d' + endstring2
    6:    tagformats[n] = ''
    7:    tagformats[n] = '(a' + endstring1
    8:    tagformats[n] = ''
    9:    tagformats[n] = ''
    10:   tagformats[n] = ''
    11:   tagformats[n] = ''
    12:   tagformats[n] = '(i' + endstring1
    13:   tagformats[n] = '(i' + endstring1
    14:   tagformats[n] = '(i' + endstring1
    15:   tagformats[n] = '(i' + endstring1
    else: tagformats[n] = ''
  endcase
  if isnull(tagformats[n]) then begin
    print,' '
    print,'ERROR in setformat: No format defined for tag ',tagnames[n], $
          ' of type ',tagtypes[n],format='(3a,i0)'
    print,' '
    return
  endif
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showformat,help=help

; List the current format for tag output

common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline

if n_params() ne 0 or keyword_set(help) then begin
  print,'showformat[,/help]'
  return
endif

setformat

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro two2one,chan,stack=n,help=help

; Take the loaded pair and load one channel as the current single spectrum;
; if the stack keyword is set, use the nth stack pair as input rather than
; the loaded pair

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,'two2one,chan[,stack=n][,/help]'
  return
endif else if chan lt 1 or chan gt 4 then begin
  print,'two2one,chan[,stack=n][,/help]'
  print,'ERROR in two2one: Must have chan = 1 (OC) or 2 (SC) or 3 or 4'
  return
endif else if n_elements(n) gt 0 then begin
  if nstack eq 0 then begin
    print,'ERROR in two2one: The pair stack is empty'
    return
  endif else if n le 0 then begin
    print,'two2one,chan[,stack=n][,/help]'
    print,'ERROR in two2one: Must have n >= 1'
    return
  endif else if n gt nstack then begin
    print,'ERROR in two2one: There are only ',nstack, $
          ' spectra in the pair stack',format='(a,i0,a)'
    return
  endif
endif else if (*loaded).ndata le 2 then begin
  print,'ERROR in two2one: No pair has been loaded'
  return
endif

polstring = ['OC','SC']

; Get some elements of the input pair

useStruc = (n_elements(n) gt 0) ? *stack[n-1] : *loaded
tags1 = useStruc.tags[chan-1]
df = tags1.dfreq
posfr = tags1.posfr
xjcen = tags1.xjcen
ndata = useStruc.ndata

; Get the frequency vector

f = (n_elements(n) gt 0) ? posfr*df*(findgen(ndata) - xjcen) : useStruc.freq

; Get the extra tags relevant to this channel

extratags = useStruc.extratags
nextra = useStruc.nextra
splitExtra,chan,extratags,nextra,extratags1,nextra1

; Can't just assign element by element, because array dimensions may differ,
; and IDL cares deeply about this for structure fields.
;
; Must instead create a complete single-channel structure, then load it

chanStruc = {freq:f, spec:reform(useStruc.spec[chan-1,*]), $
             tags:tags1, extratags:extratags1, $
             ndata:ndata, ntags:useStruc.ntags, $
             nextra:nextra1, tname:useStruc.tname, $
             pol:polstring[chan-1]}

*loaded1 = chanStruc

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro one2two,help=help

; Take the loaded single-channel spectrum and write it into
; one channel of the loaded pair.  The polarization tag is used to
; determine which channel is involved.

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 0 or keyword_set(help) then begin
  print,'one2two,[/help]'
  print,'        -- all pair elements are overwritten except the spectrum ', $
        'and tags for the other channel',format='(2a)'
  return
endif

; Check which polarization channel we're dealing with

polstring = strupcase((*loaded1).pol)
if polstring eq 'OC' then begin
  chan = 1
endif else if polstring eq 'SC' then begin
  chan = 2
endif else if polstring eq 'RE' then begin
  chan = 3
endif else if polstring eq 'IM' then begin
  chan = 4
endif else begin
  print,"ERROR in one2two: Loaded single-channel spectrum has illegal ", $
        "polarization value = '",polstring,"', so the channel can't be determined", $
        format='(4a)'
  return
endelse

; Can't check them all. If > 1, assume OC is set
if chan eq 1 then otherchan = 2 else otherchan = 1

; Since we're changing only one of the two channels, make sure that the
; new channel has the same target name, the same number of spectral points,
; the same frequency vector, and the same number of tags as the old one.
;
; EXCEPTION: If the loaded spectrum has been reset (<= 2 points) then we'll
;            just make the other channel a vector of zeros as long as the new
;            channel, with blanked-out tags to go with it.

tname1 = (*loaded1).tname
ndata1 = (*loaded1).ndata
xjcen1 = (*loaded1).tags.xjcen
posfr1 = (*loaded1).tags.posfr
dfreq1 = (*loaded1).tags.dfreq
ntags1 = (*loaded1).ntags
ndata = (*loaded).ndata
if ndata gt 2 then begin

  ; Check target name

  tname = (*loaded).tname
  if tname1 ne tname then begin
    print,"ERROR in one2two: The new channel has a different target name ('", $
          tname1,"') than the other channel ('",tname,"')",format='(5a)'
    return
  endif

  ; Check spectrum size

  if ndata1 ne ndata then begin
    print,'ERROR in one2two: The new channel has ',ndata1, $
          ' spectral points but the other channel will still have ', $
          ndata,' points',format='(a,i0,a,i0,a)'
    return
  endif

  ; Check frequency vector:
  ; Rather than directly checking the frequencies and dealing with lots of tiny
  ; floating-point differences, check the tags from which the vectors were created
  ; and deal with just one roundoff uncertainty (dfreq).

  xjcen = (*loaded).tags[otherchan-1].xjcen
  posfr = (*loaded).tags[otherchan-1].posfr
  dfreq = (*loaded).tags[otherchan-1].dfreq
  if (xjcen1 ne xjcen) or (posfr1 ne posfr) or (dfreq1 ne dfreq) then begin
    print,'ERROR in one2two: The new channel has different frequencies (xjcen=', $
          xjcen1,', posfr=',posfr1,', dfreq=',dfreq1,') than the other channel (', $
          xjcen,', ',posfr,', ',dfreq,')',format='(6(a,i0),a)'
  endif

  ; Check number of tags

  ntags = (*loaded).ntags
  if ntags1 ne ntags then begin
    print,'ERROR in one2two: The new channel has a different number of tags (', $
          ntags1,') than the other channel (',ntags,')',format='(a,i0,a,i0,a)'
    return
  endif

endif

; Create the revised pair and its revised tags

pair = fltarr(npol,ndata1)
tags = blanktags(npol=npol)

; overwrite them all with the old one if needed then copy new one in

if ndata gt 2 then begin
  pair[otherchan-1,*] = (*loaded).spec[otherchan-1,*]
  for k=0L,ntags1-1 do tags[otherchan-1].(k) = (*loaded).tags[otherchan-1].(k)
endif
pair[chan-1,*] = (*loaded1).spec
for k=0L,ntags1-1 do tags[chan-1].(k) = (*loaded1).tags.(k)

; Inspect the extra tags of the loaded pair, delete any which are specific to
; the channel to be replaced, then merge the remainder with the extra tags of
; the replacement spectrum.
;
; When inspecting and deleting, go through the list BACKWARDS so that deleting
; one extra tag doesn't change the numbers of those not yet checked.

if ndata gt 2 then begin
  oldextratags = (*loaded).extratags
  oldnextra = (*loaded).nextra
  if chan eq 1 then begin
    chanstring1 = '_OC'
    chanstring2 = '# OC: '
  endif else if chan eq 2 then begin
    chanstring1 = '_SC'
    chanstring2 = '# SC: '
  endif else if chan eq 3 then begin
    chanstring1 = '_RE'
    chanstring2 = '# RE: '
  endif else if chan eq 4 then begin
    chanstring1 = '_IM'
    chanstring2 = '# IM: '
  endif else begin
    chanstring1 = 'NONONO'
    chanstring2 = 'NONONO'
  endelse
  for k=oldnextra-1L,0,-1 do begin
    if (strpos(oldextratags[k].name, chanstring1) ne -1) or $
       (isnull(oldextratags[k].name) and $
        strpos(oldextratags[k].comment, chanstring2) eq 0)  then $
      deleteextra,(k+1),/silent
  endfor

; No good way to merge pols 3 and 4. Just don't.

  if chan lt 3 then mergeExtra,(*loaded).extratags,(*loaded1).extratags, $
             (*loaded).nextra,(*loaded1).nextra,extratags,nextra
endif else begin
  extratags = (*loaded1).extratags
  nextra = (*loaded1).nextra
endelse

; Can't just assign element by element, because array dimensions may differ,
; and IDL cares deeply about this for structure fields.
;
; Must instead create a complete pair structure, then load it

pairStruc = {freq:(*loaded1).freq, spec:pair, tags:tags, extratags:extratags, $
             ndata:ndata1, ntags:ntags1, nextra:nextra, tname:tname1}

*loaded = pairStruc

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro splice,n_OC,n_SC,n_RE,n_IM,keep=keep,silent=silent,help=help

; Combine an OC and an SC spectrum from the single-channel stack to form a loaded pair,
; then delete the two single-channel spectra from the single-channel stack.

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() gt 4 or keyword_set(help) then begin
  print,' '
  print,'splice,n_OC,n_SC[,n_RE,n_IM][,/keep][,/silent][,/help]'
  print,' '
  print,'Combine an OC and an SC spectrum from the single-channel stack (stack1)'
  print,'       to form a dual-polarization spectral pair, then load this pair'
  print,'       and delete the two single-channel spectra from stack1'
  print,' '
  print,'       n_OC and n_SC are the numbers (counting from 1) within stack1
  print,'           of the OC and SC spectra that will be spliced together,
  print,"           as displayed by the 'showstack1' procedure"
  print,' '
  print,'       /keep prevents the two single-channel spectra from being deleted'
  print,'           from stack1 after being spliced together'
  print,' '
  return
endif

; Check that the two spectra slated for splicing actually exist 
; and have the correct polarizations

npol = n_params()

nOC = long(n_OC)
nSC = long(n_SC)
if npol gt 2 then nRE = long(n_RE)
if npol gt 3 then nIM = long(n_IM)

ar = [nOC,nSC]
if npol gt 2 then ar = [ar, nRE]
if npol gt 3 then ar = [ar, nIM]

ar = ar[reverse(sort(ar))]

nmin = ar[n_elements(ar) - 1]
nmax = ar[0]

if nstack1 eq 0 then begin
  print,'ERROR in splice: There are no spectra in the single-channel stack'
  return
endif else if nmin lt 1 then begin
  print,'ERROR in splice: Spectrum #',nmin,' is out of the valid stack1 range (1 - ', $
        nstack1,')',format='(a,i0,a,i0,a)'
  return
endif else if nmax gt nstack1 then begin
  print,'ERROR in splice: Spectrum #',nmax,' is out of the valid stack1 range (1 - ', $
        nstack1,')',format='(a,i0,a,i0,a)'
  return
endif else if (*stack1[nOC-1]).pol ne 'OC' then begin
  print,'ERROR in splice: Spectrum #',nOC,' is not an OC spectrum',format='(a,i0,a)'
  return
endif else if (*stack1[nSC-1]).pol ne 'SC' then begin
  print,'ERROR in splice: Spectrum #',nSC,' is not an SC spectrum',format='(a,i0,a)'
  return
endif else if (*stack1[nRE-1]).pol ne 'RE' then begin
  print,'ERROR in splice: Spectrum #',nRE,' is not an RE spectrum',format='(a,i0,a)'
  return
endif else if (*stack1[nIM-1]).pol ne 'IM' then begin
  print,'ERROR in splice: Spectrum #',nIM,' is not an IM spectrum',format='(a,i0,a)'
  return
endif

; Store the loaded single-channel spectrum so we can use that "space" and
; reload the spectrum later on

storeloaded1 = ptr_new(*loaded1)

; Splice the channels together to form a new loaded pair

unload, npol=npol
load1,nOC
one2two
load1,nSC
one2two
if npol gt 2 then begin
  load1,nRE
  one2two
endif
if npol gt 3 then begin
  load1,nIM
  one2two
endif

; Reload the single-channel spectrum that was there at the start

*loaded1 = *storeloaded1
ptr_free,storeloaded1

; Report the results

if not keyword_set(silent) then begin
  print,'Stack1 spectra #',nOC,' and #',nSC, $
        ' have been spliced and the resulting pair has been loaded',format='(2(a,i0),a)'
endif

; Delete the two single-channel spectra from the single-channel stack
;
; Delete the higher number first so that it doesn't change the number
; within stack1 of the remaining spectrum

if not keyword_set(keep) then begin
  for i = 0, npol - 1 do begin
    deletestack1, n=ar[i],silent=silent
  endfor
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro split,chan=chan,group=gr,stack=st,silent=silent,help=help

; Push each channel of the loaded pair onto the single-channel stack,
; unless the stack keyword is set, in which case each channel of
; each pair in the pair stack is pushed on the single-channel stack

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'split[,chan=1 or 2][,group=n][,/stack][,/silent][,/help]'
  print,' '
  print,'Push each channel of the loaded pair onto the single-channel stack'
  print,'      /stack does this instead for all pairs in the pair stack'
  print,'      group=n overrides the default group number for the new stack1 spectra'
  print,'      (next available # for loaded pair, existing group number(s) for /stack)'
  print,' '
  return
endif

if not keyword_set(st) then begin

  ; Split the loaded pair

  if (*loaded).ndata le 2 then begin
    print,'ERROR in split: No pair has been loaded'
    return
  endif
endif else begin

  ; Split every pair in the stack

  if nstack eq 0 then begin
    print,'ERROR in split: The pair stack is empty'
    return
  endif
endelse

; Check which channel(s) to send to the single-channel stack

if n_elements(chan) eq 0 then begin
  splitpol = [1,1]
  reportstr1 = 'Each channel'
endif else if chan eq 1 or chan eq 2 then begin
  splitpol = (chan eq 1) ? [1,0] : [0,1]
  reportstr1 = 'Channel ' + string(chan, format='(i0)')
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse

; Store the loaded single-channel spectrum so we can use that "space" and
; reload the spectrum later on

storeloaded1 = ptr_new(*loaded1)

; Do the splitting

if keyword_set(st) then begin
  for n=0L,nstack-1 do begin
    groupnumber = (n_elements(gr) gt 0) ? long(gr) : (*stack[n]).group
    for ch=1,2 do begin
      if splitpol[ch-1] then begin
        two2one,ch,stack=(n+1)
        push1,group=groupnumber,/silent
      endif
    endfor
  endfor
  reportstr2 = ' of each stack pair'
endif else begin
  groupnumber = (n_elements(gr) gt 0) ? long(gr) : 0L
  for ch=1,2 do begin
    if splitpol[ch-1] then begin
      two2one,ch
      push1,group=groupnumber,/silent
    endif
  endfor
  reportstr2 = ' of the loaded pair'
endelse

; Reload the single-channel spectrum that was there at the start

*loaded1 = *storeloaded1
ptr_free,storeloaded1

; Report the results

if not keyword_set(silent) then $
  print,reportstr1,reportstr2,' has been added to the single-channel stack'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mergeExtra,extratags1,extratags2,nextra1,nextra2,pair_extratags,pair_nextra

; Take the "extra" (i.e., rdf footer) tags from two channels of an OC/SC pair
; and merge them into a single set of tags with no OC / SC duplication
; note: There's no good way to do this for pols 3 and 4. Don't try.

etags1 = extratags1
etags2 = extratags2
matched2 = intarr(nextra2) ; starts with all zeros ("false"s)
pair_nextra = 0L

; Go through the extra pair1 tags looking for pair2 matches

for n1=0L,nextra1-1 do begin
  nmatch = 0L
  n2_arr = where(etags1[n1].name eq etags2.name, count)
  for n=0L,count-1 do begin
    n2 = n2_arr[n]
    if (etags1[n1].format eq etags2[n2].format) and $
       (etags1[n1].value eq etags2[n2].value) and $
       (etags1[n1].comment eq etags2[n2].comment) and $
       (not matched2[n2]) then begin
      matched2[n2] = 1
      nmatch = nmatch + 1L
    endif
  endfor

  ; If no SC match, adjust the tag name to show that
  ; the value only applies to the OC spectrum.

  if nmatch eq 0 then begin
    if notnull(etags1[n1].name) then begin
      etags1[n1].name = etags1[n1].name + '_OC'
    endif else begin
      if strpos(etags1[n1].comment, '# OC: ') ne 0 then $
        etags1[n1].comment = '# OC: ' + strtrim(strmid(etags1[n1].comment, 2), 2)
    endelse
  endif

  ; Add this extra OC tag to the list of extra pair tags

  if pair_nextra eq 0L then $
    pair_extratags = [etags1[n1]] $
  else $
    pair_extratags = [pair_extratags, etags1[n1]]
  pair_nextra = pair_nextra + 1L

endfor

; Add any unmatched extra SC tags to the list of extra pair tags,
; adjusting the tag names to show that the values only apply to the
; SC spectra.

for n2=0L,nextra2-1 do begin
  if not matched2[n2] then begin
    if notnull(etags2[n2].name) then begin
      etags2[n2].name = etags2[n2].name + '_SC'
    endif else begin
      if strpos(etags2[n2].comment, '# SC: ') ne 0 then $
        etags2[n2].comment = '# SC: ' + strtrim(strmid(etags2[n2].comment, 2), 2)
    endelse
    pair_extratags = [pair_extratags, pair_extratags[0]]
    for k=0L,3 do pair_extratags[pair_nextra].(k) = etags2[n2].(k)
    pair_nextra = pair_nextra + 1L
  endif
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro splitExtra,chan,extratags,nextra,extratags1,nextra1

; Take the "extra" (i.e., rdf footer) tags from two channels of an OC/SC pair
; and split off those relevant to a particular channel

if chan eq 1 then begin
  chanstring1 = '_OC'
  otherchanstring1 = '_SC'
  chanstring2 = '# OC: '
  otherchanstring2 = '# SC: '
endif else if chan eq 2 then begin
  chanstring1 = '_SC'
  otherchanstring1 = '_OC'
  chanstring2 = '# SC: '
  otherchanstring2 = '# OC: '
endif else begin
  print,'ERROR in splitExtra: Must have chan = 1 (OC) or 2 (SC)'
  return
endelse

extratags1 = extratags    ; We'll cut it down to size at the end
nextra1 = 0L

for n=0L,nextra-1 do begin
  chanpos1 = strpos(extratags[n].name, chanstring1)
  chanpos2 = strpos(extratags[n].comment, chanstring2)
  otherchanpos1 = strpos(extratags[n].name, otherchanstring1)
  otherchanpos2 = strpos(extratags[n].comment, otherchanstring2)
  if notnull(extratags[n].name) then begin
    if chanpos1 ne -1 then begin
      extratags1[nextra1].format = extratags[n].format
      extratags1[nextra1].name = strmid(extratags[n].name, 0, chanpos1)
      extratags1[nextra1].value = extratags[n].value
      extratags1[nextra1].comment = extratags[n].comment
      nextra1 = nextra1 + 1
    endif else if otherchanpos1 eq -1 then begin
      extratags1[nextra1].format = extratags[n].format
      extratags1[nextra1].name = extratags[n].name
      extratags1[nextra1].value = extratags[n].value
      extratags1[nextra1].comment = extratags[n].comment
      nextra1 = nextra1 + 1
    endif
  endif else begin

    ; This extra tag is a comment line

    if chanpos2 eq 0 then begin
      extratags1[nextra1].format = extratags[n].format
      extratags1[nextra1].name = extratags[n].name
      extratags1[nextra1].value = extratags[n].value
      extratags1[nextra1].comment = '# ' + strtrim(strmid(extratags[n].comment, 6), 2)
      nextra1 = nextra1 + 1
    endif else if otherchanpos2 ne 0 then begin
      extratags1[nextra1].format = extratags[n].format
      extratags1[nextra1].name = extratags[n].name
      extratags1[nextra1].value = extratags[n].value
      extratags1[nextra1].comment = extratags[n].comment
      nextra1 = nextra1 + 1
    endif
  endelse
endfor

extratags1 = extratags1[0:nextra1-1]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mergeExtraStack1,chan,extratags1,nextra1,groupvec=gvec

; Take the "extra" (i.e., rdf footer) tags from all spectra
; of a particular polarization channel in the single-channel
; stack, and merge them into a single set of tags which are
; common to those spectra.
;
; gvec is a vector listing the groups within stack1 to be merged;
; if this keyword isn't used, extra tags from the entire
; single-channel stack are merged.
;
; Written to facilitate writing spectra to disk in a
; tkplay-friendly format

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

extratags1 = {format:'', name:'', value:'', comment:''}
nextra1 = 0L

if chan eq 1 then begin
  chanstring = 'OC'
endif else if chan eq 2 then begin
  chanstring = 'SC'
endif else begin
  print,'ERROR in mergeExtraStack1: Must have chan = 1 (OC) or 2 (SC)'
  return
endelse

if nstack1 eq 0 then return

if n_elements(gvec) eq 0 then gvec = stackgroups1()

; See which stack1 spectra are in the specified polarization channel
; and collect the extratag formats, names, and values for those spectra

nInChan = 0L
for n=0L,nstack1-1 do begin
  group_el = where(gvec eq (*stack1[n]).group, count)
  if count gt 0 and (*stack1[n]).pol eq chanstring then begin
    if nInChan eq 0 then begin
      indx = [n]
    endif else begin
      indx = [indx, n]
    endelse
    nInChan = nInChan + 1
  endif
endfor

if nInChan eq 0 then return

; Initialize the extratags array: set it equal to the that of
; the first stack1 spectrum in the desired polarization channel

extratags1 = (*stack1[indx[0]]).extratags
nextra1 = (*stack1[indx[0]]).nextra

; Go through the other stack1 spectra in this polarization channel
; and find those extra tags which are present and identical for all
; spectra; keep those and discard the others

for n=1L,nInChan-1 do begin
  compareformat = strlowcase(extratags1.format)
  comparename = strlowcase(extratags1.name)
  comparevalue = strlowcase(extratags1.value)
  comparecomment = strlowcase(extratags1.comment)
  eformat = strlowcase((*stack1[indx[n]]).extratags.format)
  ename = strlowcase((*stack1[indx[n]]).extratags.name)
  evalue = strlowcase((*stack1[indx[n]]).extratags.value)
  ecomment = strlowcase((*stack1[indx[n]]).extratags.comment)
  kmax = nextra1 - 1L

  ; Go through the remaining extra tags in reverse order, so that deleting
  ; one tag won't change the numbers of the ones not yet checked

  for k=kmax,0,-1L do begin
    extranum = where(ename eq comparename[k], count)
    keepextra = 0
    keepcomment = 0
    for j=0L,count-1 do begin
      match = (notnull(comparename[k]) and $
               (eformat[extranum[j]] eq compareformat[k]) and $
               (evalue[extranum[j]]  eq comparevalue[k] )       ) or $
              (isnull(comparename[k]) and (ecomment[extranum[j]] eq comparecomment[k]))
      keepextra = keepextra or match
      keepcomment = keepcomment or (match and (ecomment[extranum[j]] eq comparecomment[k]))
    endfor

    if not keepextra then begin

      ; Don't create an error by trying to remove the last array element

      if nextra1 gt 1 then begin
        if k eq nextra1-1 then begin
          extratags1 = extratags1[0:nextra1-2]
        endif else if k gt 0 then begin
          extratags1 = [extratags1[0:k-1], extratags1[k+1:nextra1-1]]
        endif else begin
          extratags1 = extratags1[1:nextra1-1]
        endelse
      endif
      nextra1 = nextra1 - 1

    endif else if not keepcomment then begin
      extratags1[k].comment = ''
    endif
  endfor

endfor

; If no extra tags are left, create a new array consisting of
; the names and numbers of the cw tags

if nextra1 eq 0 then begin
  nextra1 = (*stack1[indx[0]]).ntags
  extratag = {format:'t', name:'', value:'', comment:''}
  extratags1 = replicate(extratag, nextra1)
  extratags1.name = strlowcase(tag_names((*stack1[indx[0]]).tags))
  extratags1.value = string(indgen(nextra1),format='(i0)')
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro unload,chan=chan,npol=npol,help=help

; Unload the currently loaded pair (i.e., zero things out)

common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 0 or keyword_set(help) then begin
  print,'unload[,chan=1 .. 4][,npol=npol][,/help]'
  print,'       -- if chan is specified, only the spectrum and tags are', $
        ' reset for that channel',format='(2a)'
  return
endif
if not keyword_set(npol) then npol = 2

; Check which channel(s) to unload

if (n_elements(chan) eq 0) or (n_elements(loaded) eq 0) then begin
  unloadboth = 1
endif else if (chan eq 1 or chan eq 2 or chan eq 3 or chan eq 4) then begin
  unloadboth = 0
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or omit (both)'
  return
endelse

; Do the unloading

if unloadboth then begin

  ; Define/reset the entire pair

  ; Initialize the frequency vector and the spectral value array
  ;
  ; Need to use >= 2 frequency points, because if you use a one-element vector
  ; in a structure definition, IDL ignores you and calls it a scalar; this might
  ; cause variable type mismatch problems later on.

  ndata = 2L
  freq = fltarr(ndata)
  pair = fltarr(npol,ndata)

  ; Initialize the other pair elements

  tags = blanktags(npol=npol)
  ntags = n_tags(tags)
  nextra = ntags
  extratag = {format:'t', name:'', value:'', comment:''}
  extratags = replicate(extratag, nextra)
  extratags.name = strlowcase(tag_names(tags))
  extratags.value = string(indgen(ntags), format='(i0)')

  pairStruc = {freq:freq, spec:pair, tags:tags, extratags:extratags, $
               ndata:ndata, ntags:ntags, nextra:nextra, tname:''}

  if n_elements(loaded) eq 0 then $
    loaded = ptr_new(pairStruc) $
  else $
    *loaded = pairStruc

endif else begin

  ; Only reset the spectrum and tags for the specified channel

  (*loaded).spec[chan-1,*] = (*loaded).spec[chan-1,*]*0.0
  (*loaded).tags[chan-1] = blanktags1()

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro unload1,help=help

; Unload the currently loaded single-channel spectrum (i.e., zero things out)
  
common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 0 or keyword_set(help) then begin
  print,'unload1[,/help]'
  return
endif

; Initialize the frequency and spectral value vectors
;
; Need to use >= 2 points, because if you use a one-element vector in a
; structure definition, IDL ignores you and calls it a scalar; this might
; cause variable type mismatch problems later on.

ndata = 2L
freq = fltarr(ndata)
spec = fltarr(ndata)

; Initialize the other spectrum elements

tags1 = blanktags1()
ntags = n_tags(tags1)
nextra = ntags
extratag = {format:'t', name:'', value:'', comment:''}
extratags = replicate(extratag, nextra)
extratags.name = strlowcase(tag_names(tags1))
extratags.value = string(indgen(ntags), format='(i0)')

specStruc = {freq:freq, spec:spec, tags:tags1, extratags:extratags, $
             ndata:ndata, ntags:ntags, nextra:ntags, tname:'', pol:''}

if n_elements(loaded1) eq 0 then $
  loaded1 = ptr_new(specStruc) $
else $
  *loaded1 = specStruc

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro unloadi,help=help

; Unload the currently loaded image (i.e., zero things out)
  
common loadedBlock,loadedi,loaded1,loaded

if n_params() ne 0 or keyword_set(help) then begin
  print,'unloadi[,/help]'
  return
endif

; Initialize the height, width, and image elements
;
; Use >= 2 points per dimension to make sure that IDL doesn't ignore you by
; calling a 1x1 array a scalar; this will avoid possible variable type mismatch
; problems later on.

width = 2L
height = 2L
image = fltarr(width, height)

; Initialize the extra tag array to a dummy comment line

extratags = [{format:'', name:'', value:'', comment:'# This is a dummy comment'}]

; Create the structure and load it

imageStruc = {image:image, extratags:extratags, width:width, height:height, $
              nextra:1L, tname:'', pol:''}

if n_elements(loadedi) eq 0 then $
  loadedi = ptr_new(imageStruc) $
else $
  *loadedi = imageStruc

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro load,n,chan=chan,help=help

; Load a pair from the stack

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'load,n[,chan=1 or 2][,/help]'
  print,'     n is the number (counting from 1) within the stack ', $
        'of the pair or single channel to load',format='(2a)'
  print,' '
  return
endif else if n le 0 then begin
  print,' '
  print,'ERROR in load: n = ',n,' is illegal (must be at least 1)',format='(a,i0,a)'
  print,' '
  return
endif else if n gt nstack then begin
  print,' '
  print,'ERROR in load: n = ',n,' but there are only ',nstack, $
        ' pairs in the stack',format='(a,i0,a,i0,a)'
  print,' '
  return
endif

; Check which channel(s) to load

if n_elements(chan) eq 0 then begin
  loadboth = 1
endif else if (chan eq 1 and chan le 4) then begin
  loadboth = 0
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or 3 or 4 or omit (all)'
  return
endelse

; Create the vector of frequencies

ndata = (*stack[n-1]).ndata
xjcen = (*stack[n-1]).tags[0].xjcen
dfreq = (*stack[n-1]).tags[0].dfreq
posfr = (*stack[n-1]).tags[0].posfr
freq = posfr*dfreq*(findgen(ndata) - xjcen)

; Load the frequencies along with the spectrum and tags from the stack

if loadboth then begin
  *loaded = {freq:freq, spec:(*stack[n-1]).spec, tags:(*stack[n-1]).tags, $
             extratags:(*stack[n-1]).extratags, ndata:ndata, $
             ntags:(*stack[n-1]).ntags, nextra:(*stack[n-1]).nextra, $
             tname:(*stack[n-1]).tname}
endif else begin
  polstring = (chan eq 1) ? 'OC' : 'SC'
  *loaded1 = {freq:freq, spec:reform((*stack[n-1]).spec[chan-1,*]), $
             tags:(*stack[n-1]).tags[chan-1], $
             extratags:(*stack[n-1]).extratags, ndata:ndata, $
             ntags:(*stack[n-1]).ntags, nextra:(*stack[n-1]).nextra, $
             tname:(*stack[n-1]).tname, pol:polstring}
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro load1,n,help=help

; Load a single-channel spectrum from the single-channel stack

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'load1,n[,/help]'
  print,'     n is the number (counting from 1) within stack1 ', $
        'of the spectrum to load',format='(2a)'
  print,' '
  return
endif else if n le 0 then begin
  print,' '
  print,'ERROR in load1: n = ',n,' is illegal (must be at least 1)',format='(a,i0,a)'
  print,' '
  return
endif else if n gt nstack1 then begin
  print,' '
  print,'ERROR in load1: n = ',n,' but there are only ',nstack1, $
        ' single-channel spectra in stack1',format='(a,i0,a,i0,a)'
  print,' '
  return
endif

; Create the vector of frequencies

ndata = (*stack1[n-1]).ndata
xjcen = (*stack1[n-1]).tags.xjcen
dfreq = (*stack1[n-1]).tags.dfreq
posfr = (*stack1[n-1]).tags.posfr
freq = posfr*dfreq*(findgen(ndata) - xjcen)

; Load the frequencies along with the spectrum and tags from the single-channel stack

*loaded1 = {freq:freq, spec:(*stack1[n-1]).spec, $
            tags:(*stack1[n-1]).tags, $
            extratags:(*stack1[n-1]).extratags, ndata:ndata, $
            ntags:(*stack1[n-1]).ntags, nextra:(*stack1[n-1]).nextra, $
            tname:(*stack1[n-1]).tname, pol:(*stack1[n-1]).pol}

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro loadi,n,help=help

; Load an image from the image stack

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'loadi,n[,/help]'
  print,'     n is the number (counting from 1) within stacki ', $
        'of the image to load',format='(2a)'
  print,' '
  return
endif else if n le 0 then begin
  print,' '
  print,'ERROR in loadi: n = ',n,' is illegal (must be at least 1)',format='(a,i0,a)'
  print,' '
  return
endif else if n gt nstacki then begin
  print,' '
  print,'ERROR in loadi: n = ',n,' but there are only ',nstacki, $
        ' images in stacki',format='(a,i0,a,i0,a)'
  print,' '
  return
endif

; Take everything except the group number from the image stack

*loadedi = {image:(*stacki[n-1]).image, $
            extratags:(*stacki[n-1]).extratags, $
            width:(*stacki[n-1]).width, height:(*stacki[n-1]).height, $
            nextra:(*stacki[n-1]).nextra, $
            tname:(*stacki[n-1]).tname, pol:(*stacki[n-1]).pol}

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro extract,freq=f,spec=s,ndata=nd,tags=t,ntags=nt,tname=tn, $
            extratags=ex,nextra=nex,pol=p,chan=chan,help=help

; Extract elements of the loaded pair so that they can be manipulated in ways
; not provided by these routines

common loadedBlock,loadedi,loaded1,loaded

channames = ['OC','SC','RE','IM']

extractsome = arg_present(f) or arg_present(s) or arg_present(nd) or arg_present(t) $
                             or arg_present(nt) or arg_present(tn) or arg_present(ex) $
                             or arg_present(nex) or arg_present(p)

if keyword_set(help) or n_params() ne 0 or (not extractsome) $
                     or (arg_present(pol) and n_elements(chan) eq 0) then begin
  print,'extract,freq=freq,spec=spec,ndata=ndata,tags=tags,ntags=ntags,extratags=extratags, $'
  print,'        nextra=nextra,tname=tname[[,pol=pol],chan=1 or 2][,/help]'
  print,' '
  print,'        -- can only use the pol keyword if extracting a single channel'
  return
endif

if (*loaded).ndata le 2 then begin
  print,'ERROR in extract: No pair is loaded, so nothing can be extracted'
  return
endif

; Check which channel(s) to load

if n_elements(chan) eq 0 then begin
  extractboth = 1
endif else if (chan ge 1 and chan le 4) then begin
  extractboth = 0
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or 3 or 4 or omit (both)'
  return
endelse

; Send the requested pair elements to the outside world

if arg_present(f) then f = (*loaded).freq
if arg_present(s) then s = (extractboth) ? (*loaded).spec : reform((*loaded).spec[chan-1,*])
if arg_present(nd) then nd = (*loaded).ndata
if arg_present(t) then t = (extractboth) ? (*loaded).tags : (*loaded).tags[chan-1]
if arg_present(nt) then nt = (*loaded).ntags
if arg_present(tn) then tn = (*loaded).tname
if arg_present(ex) then ex = (*loaded).extratags
if arg_present(nex) then nex = (*loaded).nextra
if arg_present(p) then p = channames[chan]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro extract1,freq=f,spec=s,ndata=nd,tags=t,ntags=nt,tname=tn, $
             extratags=ex,nextra=nex,pol=p,help=help

; Extract elements of the loaded single-channel spectrum so that they can be manipulated
; in ways not provided by these routines

common loadedBlock,loadedi,loaded1,loaded

extractsome = arg_present(f) or arg_present(s) or arg_present(nd) or arg_present(t) $
                             or arg_present(nt) or arg_present(tn) or arg_present(ex) $
                             or arg_present(nex) or arg_present(p)

if keyword_set(help) or n_params() ne 0 or (not extractsome) then begin
  print,'extract1,freq=freq,spec=spec,ndata=ndata,tags=tags,ntags=ntags,extratags=extratags, $'
  print,'         nextra=nextra,tname=tname,pol=pol[,/help]'
  return
endif

if (*loaded1).ndata le 2 then begin
  print,'ERROR in extract1: No spectrum is loaded, so nothing can be extracted'
  return
endif

; Send the requested spectrum elements to the outside world

if arg_present(f) then f = (*loaded1).freq
if arg_present(s) then s = (*loaded1).spec
if arg_present(nd) then nd = (*loaded1).ndata
if arg_present(t) then t = (*loaded1).tags
if arg_present(nt) then nt = (*loaded1).ntags
if arg_present(tn) then tn = (*loaded1).tname
if arg_present(ex) then ex = (*loaded1).extratags
if arg_present(nex) then nex = (*loaded1).nextra
if arg_present(p) then p = (*loaded1).pol

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro extracti,image=i,width=w,height=h,tname=tn, $
             extratags=ex,nextra=nex,pol=p,help=help

; Extract elements of the loaded image so that they can be manipulated
; in ways not provided by these routines

common loadedBlock,loadedi,loaded1,loaded

extractsome = arg_present(i) or arg_present(w) or arg_present(h) or arg_present(tn)  $
                             or arg_present(ex) or arg_present(nex) or arg_present(p)

if keyword_set(help) or n_params() ne 0 or (not extractsome) then begin
  print,'extracti,image=image,width=width,height=height,extratags=extratags, $'
  print,'         nextra=nextra,tname=tname,pol=pol[,/help]'
  return
endif

if (*loadedi).width le 2 and (*loadedi).height le 2 then begin
  print,'ERROR in extracti: No image is loaded, so nothing can be extracted'
  return
endif

; Send the requested image elements to the outside world

if arg_present(i) then i = (*loadedi).image
if arg_present(w) then w = (*loadedi).width
if arg_present(h) then h = (*loadedi).height
if arg_present(tn) then tn = (*loadedi).tname
if arg_present(ex) then ex = (*loadedi).extratags
if arg_present(nex) then nex = (*loadedi).nextra
if arg_present(p) then p = (*loadedi).pol

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro changespec,freq=f,spec=s,tags=t,tname=tn,extratags=ex, $
               chan=chan,reset=reset,help=help

; Input one or more pair elements which will replace the existing elements of
; the loaded pair; can optionally reset (unload) the pair first.
;
; If the user has extracted the loaded pair (procedure 'extract') in order to
; carry out manipulations not provided by these routines, changespec allows the
; result(s) to be loaded prior to further processing via these routines.

common loadedBlock,loadedi,loaded1,loaded

somechanges = arg_present(f) or arg_present(s)  or arg_present(t) $
                             or arg_present(tn) or arg_present(ex)

if keyword_set(help) or n_params() ne 0 or not somechanges then begin
  print,'changespec,freq=freq,spec=spec,tags=tags,extratags=extratags,', $
        'tname=tname[,chan=1 or 2][,/reset,/help]',format='(2a)'
  print,'        -- /reset resets (unloads) the loaded pair or channel prior', $
        ' to changing it',format='(2a)'
  print,'        -- the chan keyword only affects spec, tags, and /reset'
  return
endif

; Check which channel(s) to change

if n_elements(chan) eq 0 then begin
  changeboth = 1
endif else if (chan ge 1 and chan le 2) then begin
  changeboth = 0
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or 3 or 4 or omit (both)'
  return
endelse

; Do whatever unloading is required/requested; first store the loaded pair
; in case we find a problem with the input and we have to reload it.

storeloaded = ptr_new(*loaded)

if keyword_set(reset) then begin
  if changeboth then unload else unload,chan=chan
endif

; Compute some parameters which will be needed for input validation

fsize = (arg_present(f)) ? n_elements(f) : n_elements((*loaded).freq)
specinfo = (arg_present(s)) ? size(s) : size((*loaded).spec)
specdims = specinfo[0]
specsize = (specdims le 1) ? specinfo[1] : specinfo[2]
lspecinfo = size((*loaded).spec)
lspecdims = lspecinfo[0]
lspecsize = (lspecdims le 1) ? specinfo[1] : specinfo[2]
tagsize = (arg_present(t)) ? n_elements(t) : 2L
ltagsize = n_elements((*loaded).tags[0,*])

; Go through each argument present, check that there are no mismatches between
; dimensions or array lengths.  If mismatches are found, reload the initial
; loaded pair and return; otherwise set the appropriate pair elements.
 
if arg_present(f) then begin
  if fsize ne specsize then begin
    print,'ERROR in changespec: Mismatch between ',fsize,' frequencies and ',specsize, $
          ' spectral points',format='(a,i0,a,i0,a)'
    *loaded = *storeloaded
    ptr_free,storeloaded
    return
  endif else if fsize lt 2 then begin
    print,'ERROR in changespec: Need at least two frequencies'
    *loaded = *storeloaded
    ptr_free,storeloaded
    return
  endif
  (*loaded).freq = f
  (*loaded).ndata = fsize
endif

if arg_present(s) then begin
  if specsize ne fsize then begin
    print,'ERROR in changespec: Mismatch between ',fsize,' frequencies and ',specsize, $
          ' spectral points',format='(a,i0,a,i0,a)'
    *loaded = *storeloaded
    ptr_free,storeloaded
    return
  endif else if specsize lt 2 then begin
    print,'ERROR in changespec: Need at least two spectral values per channel'
    *loaded = *storeloaded
    ptr_free,storeloaded
    return
  endif else if (not changeboth) and (specsize ne lspecsize) then begin
    print,'ERROR in changespec: Mismatch between ',lspecsize, $
          ' spectral points in loaded pair and ',specsize, $
          ' points in the new channel',format='(a,i0,a,i0,a)'
    *loaded = *storeloaded
    ptr_free,storeloaded
    return
  endif else if changeboth and (specdims ne lspecdims) then begin
    print,"ERROR in changespec: Can't replace a two-channel spectral array with a 1-D vector"
    *loaded = *storeloaded
    ptr_free,storeloaded
    return
  endif
  if changeboth then begin
    (*loaded).spec = s
    (*loaded).ndata = specsize
  endif else begin
    (*loaded).spec[chan-1,*] = (specdims le 1) ? s[*] : s[chan-1,*]
  endelse
endif

if arg_present(t) then begin
  if n_tags(t) ne (*loaded).ntags then begin
    print,'ERROR in changespec: These routines are only set up to handle ',(*loaded).ntags, $
          'tags, not ',n_tags(t),format='(a,i0,a,i0)'
    *loaded = *storeloaded
    ptr_free,storeloaded
    return
  endif else if changeboth and (tagsize ne ltagsize) then begin
    print,"ERROR in changespec: Can't replace a two-channel tag-structure vector with a ", $
          "one-channel structure",format='(2a)'
    *loaded = *storeloaded
    ptr_free,storeloaded
    return
  endif
  if changeboth then begin
    (*loaded).tags = t
  endif else begin
    (*loaded).tags[chan-1] = (tagsize le 1) ? t : t[chan-1]
  endelse
endif

if arg_present(tn) then (*loaded).tname = tn

if arg_present(ex) then begin
  (*loaded).extratags = ex
  (*loaded).nextra = n_elements(ex)
endif

ptr_free,storeloaded

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro changespec1,freq=f,spec=s,tags=t,tname=tn,extratags=ex, $
                reset=reset,help=help

; Input one or more spectrum elements which will replace the existing elements of
; the loaded single-channel spectrum; can optionally reset (unload) the spectrum first.
;
; If the user has extracted the loaded spectrum (procedure 'extract1') in order to
; carry out manipulations not provided by these routines, changespec1 allows the
; result(s) to be loaded prior to further processing via these routines.

common loadedBlock,loadedi,loaded1,loaded

somechanges = arg_present(f) or arg_present(s)  or arg_present(t) $
                             or arg_present(tn) or arg_present(ex)

if keyword_set(help) or n_params() ne 0 or not somechanges then begin
  print,'changespec1,freq=freq,spec=spec,tags=tags,extratags=extratags,', $
        'tname=tname[,/reset,/help]',format='(2a)'
  print,'        -- /reset resets (unloads) the loaded spectrum prior', $
        ' to changing it',format='(2a)'
  return
endif

; Do whatever unloading is required/requested; first store the loaded spectrum
; in case we find a problem with the input and we have to reload it.

storeloaded1 = ptr_new(*loaded1)

if keyword_set(reset) then unload1

; Compute some parameters which will be needed for input validation

fsize = (arg_present(f)) ? n_elements(f) : n_elements((*loaded1).freq)
specsize = (arg_present(s)) ? n_elements(s) : n_elements((*loaded1).spec)
lspecsize = n_elements((*loaded1).spec)

; Go through each argument present, check that there are no mismatches between
; dimensions or array lengths.  If mismatches are found, reload the initial
; loaded spectrum and return; otherwise set the appropriate spectrum elements.
 
if arg_present(f) then begin
  if fsize ne specsize then begin
    print,'ERROR in changespec1: Mismatch between ',fsize,' frequencies and ',specsize, $
          ' spectral points',format='(a,i0,a,i0,a)'
    *loaded1 = *storeloaded1
    ptr_free,storeloaded1
    return
  endif else if fsize lt 2 then begin
    print,'ERROR in changespec1: Need at least two frequencies'
    *loaded1 = *storeloaded1
    ptr_free,storeloaded1
    return
  endif
  (*loaded1).freq = f
  (*loaded1).ndata = fsize
endif

if arg_present(s) then begin
  if specsize ne fsize then begin
    print,'ERROR in changespec1: Mismatch between ',fsize,' frequencies and ',specsize, $
          ' spectral points',format='(a,i0,a,i0,a)'
    *loaded1 = *storeloaded1
    ptr_free,storeloaded1
    return
  endif else if specsize lt 2 then begin
    print,'ERROR in changespec1: Need at least two spectral values'
    *loaded1 = *storeloaded1
    ptr_free,storeloaded1
    return
  endif else if specsize ne lspecsize then begin
    print,'ERROR in changespec1: Mismatch between ',lspecsize, $
          ' spectral points in loaded spectrum and ',specsize, $
          ' points in the new spectrum',format='(a,i0,a,i0,a)'
    *loaded1 = *storeloaded1
    ptr_free,storeloaded1
    return
  endif
  (*loaded1).spec = s
  (*loaded1).ndata = specsize
endif

if arg_present(t) then begin
  if n_tags(t) ne (*loaded1).ntags then begin
    print,'ERROR in changespec1: These routines are only set up to handle ', $
          (*loaded1).ntags,'tags, not ',n_tags(t),format='(a,i0,a,i0)'
    *loaded1 = *storeloaded1
    ptr_free,storeloaded1
    return
  endif
  (*loaded1).tags = t
endif

if arg_present(tn) then (*loaded1).tname = tn

if arg_present(ex) then begin
  (*loaded1).extratags = ex
  (*loaded1).nextra = n_elements(ex)
endif

ptr_free,storeloaded1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro changeim,image=i,tname=tn,extratags=ex,reset=reset,help=help

; Input one or more image elements which will replace the existing elements of
; the loaded image; can optionally reset (unload) the image first.
;
; If the user has extracted the loaded image (procedure 'extracti') in order to
; carry out manipulations not provided by these routines, changeim allows the
; result(s) to be loaded prior to further processing via these routines.

common loadedBlock,loadedi,loaded1,loaded

somechanges = arg_present(i) or arg_present(tn) or arg_present(ex)

if keyword_set(help) or n_params() ne 0 or not somechanges then begin
  print,'changeim,image=image,extratags=extratags,tname=tname[,/reset,/help]'
  print,'     -- /reset resets (unloads) the loaded image prior', $
        ' to changing it',format='(2a)'
  return
endif

; Do whatever unloading is required/requested; first store the loaded image
; in case we find a problem with the input and we have to reload it.

storeloadedi = ptr_new(*loadedi)

if keyword_set(reset) then unloadi

; Check that that the input array (if specified) is 2-D and has at least two
; elements in each dimension.  If there's a problem here, reload the initial
; loaded image and return; otherwise set the appropriate image elements.
 
if arg_present(i) then begin
  imageinfo = size(i)
  imagedims = imageinfo[0]
  width = imageinfo[1]
  height = (imagedims ge 2) ? imageinfo[2] : -1L
  if imagedims ne 2 then begin
    print,'ERROR in changeim: Input image is ',imagedims, $
          '-D, but only 2-D images are allowed',format='(a,i0,a)'
    *loadedi = *storeloadedi
    ptr_free,storeloadedi
    return
  endif else if width lt 2 or height lt 2 then begin
    print,'ERROR in changeim: Input image must be >= 2 pixels wide', $
          ' in each dimension',format='(2a)'
    *loadedi = *storeloadedi
    ptr_free,storeloadedi
    return
  endif
  (*loadedi).image = i
  (*loadedi).width = width
  (*loadedi).height = height
endif

if arg_present(tn) then (*loadedi).tname = tn

if arg_present(ex) then begin
  (*loadedi).extratags = ex
  (*loadedi).nextra = n_elements(ex)
endif

ptr_free,storeloadedi

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro push,chan=chan,group=gr,silent=silent,help=help

; Push (add) the loaded pair onto the end of the stack

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() ne 0 then begin
  print,'push,[,chan=1 or 2][,group=n][,/silent][,/help]'
  return
endif

if (*loaded).ndata le 2 then begin
  print,'ERROR in push: No pair is loaded, so nothing can be added to the stack'
  return
endif

; Check which channel(s) to push

if n_elements(chan) eq 0 then begin
  pushboth = 1
endif else if (chan ge 1 and chan le 4) then begin
  pushboth = 0
endif else begin
  print,'Must use chan = 1 (OC) or 2 (SC) or 3 or 4 or omit (both)'
  return
endelse

; Add the desired pair element(s) to the appropriate stack

if pushboth then begin

  ; Get a group number

  if n_elements(gr) gt 0 then begin
    groupToUse = long(gr)
  endif else begin
    groupToUse = (nstack gt 0) ? max(stackgroups()) + 1L : 1L
  endelse

  ; Create the structure which will be added to the pair stack

  stackStruc = {group:groupToUse, spec:(*loaded).spec, $
                tags:(*loaded).tags, extratags:(*loaded).extratags, $
                ndata:(*loaded).ndata, ntags:(*loaded).ntags, $
                nextra:(*loaded).nextra, tname:(*loaded).tname}

  ; Add the loaded pair to the pair stack

  if nstack eq 0 then begin
    stack = ptrarr(1, /allocate_heap)
    *stack[0] = stackStruc
  endif else begin
    stack = [stack, ptr_new(stackStruc)]
  endelse
  nstack = nstack + 1L
  if not keyword_set(silent) then $
    print,'Added pair #',nstack,' to the pair stack',format='(a,i0,a)'

endif else begin

  ; Store the loaded single-channel spectrum so we can use that "space" and
  ; reload the spectrum later on

  storeloaded1 = ptr_new(*loaded1)

  ; Load the desired channel of the loaded pair, then push it onto the
  ; single-channel stack

  two2one,chan
  push1,silent=silent

  ; Reload the single-channel spectrum that was there at the start

  *loaded1 = *storeloaded1
  ptr_free,storeloaded1

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro push1,group=gr,silent=silent,help=help

; Push (add) the loaded single-channel spectrum onto the end of
; the single-channel stack

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() ne 0 then begin
  print,'push1[,group=n],/silent][,/help]'
  return
endif

if (*loaded1).ndata le 2 then begin
  print,'ERROR in push1: No single-channel spectrum is loaded, so nothing can be', $
        'added to stack1',format='(2a)'
  return
endif

; Add the loaded single-channel spectrum to the end of stack1;
; start by getting a group number

if n_elements(gr) gt 0 then begin
  groupToUse = long(gr)
endif else begin
  groupToUse = (nstack1 gt 0) ? max(stackgroups1()) + 1L : 1L
endelse

; Create the structure which will be added to the single-channel stack

stack1Struc = {group:groupToUse, spec:(*loaded1).spec, tags:(*loaded1).tags, $
               extratags:(*loaded1).extratags, ndata:(*loaded1).ndata, $
               ntags:(*loaded1).ntags, nextra:(*loaded1).nextra, $
               tname:(*loaded1).tname, pol:(*loaded1).pol}

; Add it to stack1

if nstack1 eq 0 then begin
  stack1 = ptrarr(1, /allocate_heap)
  *stack1[0] = stack1Struc
endif else begin
  stack1 = [stack1, ptr_new(stack1Struc)]
endelse
nstack1 = nstack1 + 1L
if not keyword_set(silent) then $
  print,'Added spectrum #',nstack1,' to the single-channel stack',format='(a,i0,a)'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro pushi,group=gr,silent=silent,help=help

; Push (add) the loaded image onto the end of the image stack

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() ne 0 then begin
  print,'pushi[,group=n],/silent][,/help]'
  return
endif

if (*loadedi).width le 2 or (*loadedi).height le 2 then begin
  print,'ERROR in pushi: No image is loaded, so nothing can be', $
        'added to stacki',format='(2a)'
  return
endif

; Add the loaded image to the end of stacki;
; start by getting a group number

if n_elements(gr) gt 0 then begin
  groupToUse = long(gr)
endif else begin
  groupToUse = (nstacki gt 0) ? max(stackgroupsi()) + 1L : 1L
endelse

; Create the structure which will be added to the image stack

stackiStruc = {group:groupToUse, image:(*loadedi).image, $
               extratags:(*loadedi).extratags, width:(*loadedi).width, $
               height:(*loadedi).height, nextra:(*loadedi).nextra, $
               tname:(*loadedi).tname, pol:(*loadedi).pol}

; Add it to stacki

if nstacki eq 0 then begin
  stacki = ptrarr(1, /allocate_heap)
  *stacki[0] = stackiStruc
endif else begin
  stacki = [stacki, ptr_new(stackiStruc)]
endelse
nstacki = nstacki + 1L
if not keyword_set(silent) then $
  print,'Added image #',nstacki,' to the image stack',format='(a,i0,a)'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro deletestack,n=n,first=first,last=last,arr=arr,silent=silent,help=help

; Delete one or more pairs from the stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that a valid combination of keywords has been used

if n_elements(n) gt 0 then begin
  keywordsOK = n_elements(first) eq 0 and n_elements(last) eq 0 and n_elements(arr) eq 0
endif else if n_elements(first) gt 0 then begin
  keywordsOK = n_elements(last) gt 0 and n_elements(n) eq 0 and n_elements(arr) eq 0
endif else if n_elements(last) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(arr) gt 0 then begin
  keywordsOK = 1
endif else begin
  keywordsOK = 0
endelse

if n_params() ne 0 or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'deletestack[,n=n][,first=first,last=last][,arr=arr][,/silent][,/help]'
  print,'   Specify one and only one of the following:
  print,'   -  n is the number (counting from 1) of the pair to delete from the stack'
  print,'   -  first and last together specify a range of pairs to delete'
  print,'      (last will be truncated if it points beyond the end of the stack)'
  print,'   -  arr is an integer vector containing the numbers of the pairs to delete'
  print,' '
  return
endif

; Create a sorted vector containing the pair numbers to be deleted from the stack

if n_elements(n) gt 0 then begin
  ndel = 1L
  deletevec = [n]
endif else if n_elements(first) gt 0 then begin
  uselast = last < nstack
  if first gt last then begin
    print,'ERROR in deletestack: Must specify range with last >= first'
    return
  endif else if first gt uselast then begin
    print,"ERROR in deletestack: Can't delete pairs ",first,'-',last,', since', $
          format='(a,i0,a,i0,a)'
    print,'                      there are only ',nstack, $
          ' pairs in the stack',format='(a,i0,a)'
    return
  endif else begin
    ndel = uselast - first + 1L
    deletevec = lindgen(ndel) + first
  endelse
endif else begin
  ndel = n_elements(arr)
  deletevec = arr[ uniq(arr, sort(arr)) ]
  if n_elements(deletevec) lt ndel then begin
    print,"ERROR in deletestack: Don't include duplicate pair numbers in array arr"
    return
  endif
endelse

; Check that all pairs slated for deletion actually exist

nmin = min(deletevec)
nmax = max(deletevec)
if nstack eq 0 then begin
  print,'ERROR in deletestack: There are no pairs in the stack, so none can be deleted'
  return
endif else if nmin lt 1 then begin
  print,'ERROR in deletestack: Pair #',nmin,' is out of the valid stack range (1 - ', $
        nstack,')',format='(a,i0,a,i0,a)'
  return
endif else if nmax gt nstack then begin
  print,'ERROR in deletestack: Pair #',nmax,' is out of the valid stack range (1 - ', $
        nstack,')',format='(a,i0,a,i0,a)'
  return
endif

; Delete the specified pairs:
; Go through the sorted list of pair numbers BACKWARDS so that deleting one pair doesn't
; change the number within the stack of the remaining pairs

for j=ndel,1,-1L do begin
  k = deletevec(j-1)
  ptr_free,stack[k-1]
  if nstack gt 1 then begin
    if k eq nstack then begin
      stack = stack[0:nstack-2]
    endif else if k eq 1 then begin
      stack = stack[1:nstack-1]
    endif else begin
      stack = [stack[0:k-2], stack[k:nstack-1]]
    endelse
  endif
  nstack = nstack - 1L
endfor

; Report the results

if not keyword_set(silent) then begin
  print,' '
  formatstring = '(a,i0,a,' + string(ndel,format='(i0)') + '(2x,i0))'
  print,'Deleted ',ndel,' pairs from the stack:',deletevec,format=formatstring
  print,'Stack now contains ',nstack,' pairs',format='(a,i0,a)'
  print,' '
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro deletestack1,n=n,first=first,last=last,arr=arr,silent=silent,help=help

; Delete one or more spectra from the single-channel stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that a valid combination of keywords has been used

if n_elements(n) gt 0 then begin
  keywordsOK = n_elements(first) eq 0 and n_elements(last) eq 0 and n_elements(arr) eq 0
endif else if n_elements(first) gt 0 then begin
  keywordsOK = n_elements(last) gt 0 and n_elements(n) eq 0 and n_elements(arr) eq 0
endif else if n_elements(last) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(arr) gt 0 then begin
  keywordsOK = 1
endif else begin
  keywordsOK = 0
endelse

if n_params() ne 0L or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'deletestack1[,n=n][,first=first,last=last][,arr=arr][,/silent][,/help]'
  print,'   Specify one and only one of the following:
  print,'   -  n is the number (counting from 1) of the single-channel spectrum ', $
        'to delete from stack1',format='(2a)'
  print,'   -  first and last together specify a range of spectra to delete'
  print,'      (last will be truncated if it points beyond the end of the stack)'
  print,'   -  arr is an integer vector containing the numbers of the spectra to delete'
  print,' '
  return
endif

; Create a sorted vector containing the spectrum numbers to be deleted from the stack

if n_elements(n) gt 0 then begin
  ndel = 1L
  deletevec = [n]
endif else if n_elements(first) gt 0 then begin
  uselast = last < nstack1
  if first gt last then begin
    print,'ERROR in deletestack1: Must specify range with last >= first'
    return
  endif else if first gt uselast then begin
    print,"ERROR in deletestack1: Can't delete spectra ",first,'-',last,', since', $
          format='(a,i0,a,i0,a)'
    print,'                      there are only ',nstack1, $
          ' spectra in stack1',format='(a,i0,a)'
    return
  endif else begin
    ndel = uselast - first + 1L
    deletevec = lindgen(ndel) + first
  endelse
endif else begin
  ndel = n_elements(arr)
  deletevec = arr[ uniq(arr, sort(arr)) ]
  if n_elements(deletevec) lt ndel then begin
    print,"ERROR in deletestack1: Don't include duplicate spectrum numbers in array arr"
    return
  endif
endelse

; Check that all pairs slated for deletion actually exist

nmin = min(deletevec)
nmax = max(deletevec)
if nstack1 eq 0 then begin
  print,'ERROR in deletestack1: There are no spectra in the single-channel stack, ', $
        'so none can be deleted',format='(2a)'
  return
endif else if nmin lt 1 then begin
  print,'ERROR in deletestack1: Spectrum #',nmin,' is out of the valid stack1 range (1 - ', $
        nstack1,')',format='(a,i0,a,i0,a)'
  return
endif else if nmax gt nstack1 then begin
  print,'ERROR in deletestack1: Spectrum #',nmax,' is out of the valid stack1 range (1 - ', $
        nstack1,')',format='(a,i0,a,i0,a)'
  return
endif

; Delete the specified single-channel spectra:
; Go through the sorted list of spectrum numbers BACKWARDS so that deleting one spectrum doesn't
; change the number within stack1 of the remaining spectra

for j=ndel,1,-1L do begin
  k = deletevec(j-1)
  ptr_free,stack1[k-1]
  if nstack1 gt 1 then begin
    if k eq nstack1 then begin
      stack1 = stack1[0:nstack1-2]
    endif else if k eq 1 then begin
      stack1 = stack1[1:nstack1-1]
    endif else begin
      stack1 = [stack1[0:k-2], stack1[k:nstack1-1]]
    endelse
  endif
  nstack1 = nstack1 - 1L
endfor

; Report the results

if not keyword_set(silent) then begin
  print,' '
  formatstring = '(a,i0,a,' + string(ndel,format='(i0)') + '(2x,i0))'
  print,'Deleted ',ndel,' spectra from the single-channel stack:',deletevec,format=formatstring
  print,'Stack1 now contains ',nstack1,' spectra',format='(a,i0,a)'
  print,' '
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro deletestacki,n=n,first=first,last=last,arr=arr,silent=silent,help=help

; Delete one or more images from the image stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

; Check that a valid combination of keywords has been used

if n_elements(n) gt 0 then begin
  keywordsOK = n_elements(first) eq 0 and n_elements(last) eq 0 and n_elements(arr) eq 0
endif else if n_elements(first) gt 0 then begin
  keywordsOK = n_elements(last) gt 0 and n_elements(n) eq 0 and n_elements(arr) eq 0
endif else if n_elements(last) gt 0 then begin
  keywordsOK = 0
endif else if n_elements(arr) gt 0 then begin
  keywordsOK = 1
endif else begin
  keywordsOK = 0
endelse

if n_params() ne 0L or (not keywordsOK) or keyword_set(help) then begin
  print,' '
  print,'deletestacki[,n=n][,first=first,last=last][,arr=arr][,/silent][,/help]'
  print,'   Specify one and only one of the following:
  print,'   -  n is the number (counting from 1) of the image ', $
        'to delete from stacki',format='(2a)'
  print,'   -  first and last together specify a range of images to delete'
  print,'      (last will be truncated if it points beyond the end of the stack)'
  print,'   -  arr is an integer vector containing the numbers of the images to delete'
  print,' '
  return
endif

; Create a sorted vector containing the image numbers to be deleted from the stack

if n_elements(n) gt 0 then begin
  ndel = 1L
  deletevec = [n]
endif else if n_elements(first) gt 0 then begin
  uselast = last < nstacki
  if first gt last then begin
    print,'ERROR in deletestacki: Must specify range with last >= first'
    return
  endif else if first gt uselast then begin
    print,"ERROR in deletestacki: Can't delete images ",first,'-',last,', since', $
          format='(a,i0,a,i0,a)'
    print,'                      there are only ',nstacki, $
          ' images in stacki',format='(a,i0,a)'
    return
  endif else begin
    ndel = uselast - first + 1L
    deletevec = lindgen(ndel) + first
  endelse
endif else begin
  ndel = n_elements(arr)
  deletevec = arr[ uniq(arr, sort(arr)) ]
  if n_elements(deletevec) lt ndel then begin
    print,"ERROR in deletestacki: Don't include duplicate image numbers in array arr"
    return
  endif
endelse

; Check that all images slated for deletion actually exist

nmin = min(deletevec)
nmax = max(deletevec)
if nstacki eq 0 then begin
  print,'ERROR in deletestacki: There are no images in the image stack, ', $
        'so none can be deleted',format='(2a)'
  return
endif else if nmin lt 1 then begin
  print,'ERROR in deletestacki: Image #',nmin,' is out of the valid stacki range (1 - ', $
        nstacki,')',format='(a,i0,a,i0,a)'
  return
endif else if nmax gt nstacki then begin
  print,'ERROR in deletestacki: Image #',nmax,' is out of the valid stacki range (1 - ', $
        nstacki,')',format='(a,i0,a,i0,a)'
  return
endif

; Delete the specified image(s):
; Go through the sorted list of image numbers BACKWARDS so that deleting one image doesn't
; change the number within stacki of the remaining images

for j=ndel,1,-1L do begin
  k = deletevec(j-1)
  ptr_free,stacki[k-1]
  if nstacki gt 1 then begin
    if k eq nstacki then begin
      stacki = stacki[0:nstacki-2]
    endif else if k eq 1 then begin
      stacki = stacki[1:nstacki-1]
    endif else begin
      stacki = [stacki[0:k-2], stacki[k:nstacki-1]]
    endelse
  endif
  nstacki = nstacki - 1L
endfor

; Report the results

if not keyword_set(silent) then begin
  print,' '
  formatstring = '(a,i0,a,' + string(ndel,format='(i0)') + '(2x,i0))'
  print,'Deleted ',ndel,' images from the image stack:',deletevec,format=formatstring
  print,'Stacki now contains ',nstacki,' images',format='(a,i0,a)'
  print,' '
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showstack,page=page,first=nfirst,help=help

; Display all pairs in the pair stack
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

maxnamelen = 16
month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

if n_params() gt 0 or keyword_set(help) then begin
  print,'showstack[,page=n][,first=nfirst][,/help]'
  print,' '
  print,'n is the number of lines to display per page (default = all lines at once)'
  print,' '
  print,'nfirst is the first stack pair to display, counting from 1 (default = 1)'
  print,' '
  return
endif

if n_elements(page) eq 0 then page = nstack
if n_elements(nfirst) eq 0 then nfirst = 1L

; Check that there are some pairs in the stack to display

if nstack eq 0 then begin
  print,'The pair stack is empty'
  return
endif else if nstack lt nfirst then begin
  print,'There are only ',nstack,' pairs in the stack',format='(a,i0,a)'
  return
endif

; Format the header line

if nstack lt 1000 then begin
  numheadformat = '(a3,3x,'
  numformat = '(i3,3x,'
endif else begin
  numheadformat = '(a4,2x,'
  numformat = '(i4,2x,'
endelse
headformatstring = numheadformat + 'a,' + string((maxnamelen-4),format='(i0)') + $
                   'x,8x,a,7x,a,5x,a,2x,a,2x,a,2x,a,2x,a)'

n_so_far = nfirst - 1
userreply = ''

while n_so_far lt nstack and userreply ne 'q' and userreply ne 'Q' do begin

  ; Print the header line

  print,' '
  print,'Num','Name','Date and time','Zone','npts','dfreq (Hz)','sdev (km^2)', $
        'phase (deg)','Group',format=headformatstring

  ; Display information on each pair

  n1 = n_so_far
  n2 = (n_so_far + page - 1) < (nstack - 1)
  for n=n1,n2 do begin
    tname = (*stack[n]).tname
    namelen = strlen(tname)
    if namelen lt maxnamelen then begin
      nameformat = 'a,' + string((maxnamelen-namelen),format='(i0)') + 'x'
    endif else begin
      nameformat = 'a' + string(maxnamelen,format='(i0)')
      print,"showstack: Target name '",tname,"' will be truncated at ",maxnamelen, $
            " characters",format='(3a,i0,a)'
    endelse
    group = (*stack[n]).group
    ndata = (*stack[n]).ndata
    tags = (*stack[n]).tags
    extratags = (*stack[n]).extratags
    tznum = where(extratags.name eq 'timezone', count)
    if count gt 0 then begin
      tznum = tznum[0]
      timezone = extratags[tznum].value
    endif else begin
      timezone = '   '
    endelse
    dfreq = tags[0].dfreq
    sdev = tags[0].sdev
    sdev_format = (sdev ge 1.0) ? 'f11.4,1x' : 'f12.10'
    phase = tags[0].phase
    if (tags[0].iyy gt 0) then begin
      jd = julday(tags[0].imm, tags[0].idd, tags[0].iyy,   $
                  tags[0].rchour, tags[0].rcmin,           $
                  (tags[0].rcsec + tags[0].rcnsec/1.0e9) )
      caldat_roundsec,jd,mon,day,year,hour,min,sec
      if year eq -4713 then year = 0  ;  in case Julian date was set to 0
    endif else begin
      year = 1
      mon = 1
      day = 1
      hour = 0
      min = 0
      sec = 0
    endelse
    formatstring = numformat + nameformat + $
                   ',3x,i4,1x,a3,1x,i2.2,4x,i2.2,a1,i2.2,a1,i2.2,3x,a3,3x' + $
                   ',i6,3x,f7.3,3x,' + sdev_format + ',3x,f6.1,5x,i4)'
    print,n+1,tname,year,month[mon-1],day,hour,':',min,':',sec,timezone, $
          ndata,dfreq,sdev,phase,group,format=formatstring
  endfor
  n_so_far = n2 + 1
  if n_so_far lt nstack then begin
    print,' '
    pause,prompt='Press Enter to continue or q to quit: ',reply=userreply
  endif
endwhile

print,' '

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showstack1,page=page,first=nfirst,help=help

; Display all spectra in the single-channel stack
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

maxnamelen = 16
month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

if n_params() gt 0 or keyword_set(help) then begin
  print,'showstack1[,page=n][,first=nfirst][,/help]'
  print,' '
  print,'n is the number of lines to display per page (default = all lines at once)'
  print,' '
  print,'nfirst is the first stack1 spectrum to display, counting from 1 (default = 1)'
  print,' '
  return
endif

if n_elements(page) eq 0 then page = nstack1
if n_elements(nfirst) eq 0 then nfirst = 1L

; Check that there are some spectra in the single-channel stack to display

if nstack1 eq 0 then begin
  print,'The single-channel stack is empty'
  return
endif else if nstack1 lt nfirst then begin
  print,'There are only ',nstack1,' spectra in the single-channel stack', $
        format='(a,i0,a)'
  return
endif

; Format the header line

if nstack1 lt 1000 then begin
  numheadformat = '(a3,3x,'
  numformat = '(i3,3x,'
endif else begin
  numheadformat = '(a4,2x,'
  numformat = '(i4,2x,'
endelse
headformatstring = numheadformat + 'a,' + string((maxnamelen-4),format='(i0)') + $
                   'x,3x,a,9x,a,7x,a,5x,a,2x,a,2x,a,2x,a,2x,a)'

n_so_far = nfirst - 1
userreply = ''

while n_so_far lt nstack1 and userreply ne 'q' and userreply ne 'Q' do begin

  ; Print the header line

  print,' '
  print,'Num','Name','Pol','Date and time','Zone','npts','dfreq (Hz)','sdev (km^2)', $
        'phase (deg)','Group',format=headformatstring

  ; Display information on each pair

  n1 = n_so_far
  n2 = (n_so_far + page - 1) < (nstack1 - 1)
  for n=n1,n2 do begin
    tname = (*stack1[n]).tname
    namelen = strlen(tname)
    if namelen lt maxnamelen then begin
      nameformat = 'a,' + string((maxnamelen-namelen),format='(i0)') + 'x'
    endif else begin
      nameformat = 'a' + string(maxnamelen,format='(i0)')
      print,"showstack1: Target name '",tname,"' will be truncated at ",maxnamelen, $
            " characters",format='(3a,i0,a)'
    endelse
    group = (*stack1[n]).group
    ndata = (*stack1[n]).ndata
    pol = (*stack1[n]).pol
    tags1 = (*stack1[n]).tags
    extratags = (*stack1[n]).extratags
    tznum = where(extratags.name eq 'timezone', count)
    if count gt 0 then begin
      tznum = tznum[0]
      timezone = extratags[tznum].value
    endif else begin
      timezone = '   '
    endelse
    dfreq = tags1.dfreq
    sdev = tags1.sdev
    sdev_format = (sdev ge 1.0) ? 'f11.4,1x' : 'f12.10'
    phase = tags1.phase
    if (tags1.iyy gt 0) then begin
      jd = julday(tags1.imm, tags1.idd, tags1.iyy,     $
                  tags1.rchour, tags1.rcmin,           $
                  (tags1.rcsec + tags1.rcnsec/1.0e9) )
      caldat_roundsec,jd,mon,day,year,hour,min,sec
      if year eq -4713 then year = 0  ;  in case Julian date was set to 0
    endif else begin
      year = 1
      mon = 1
      day = 1
      hour = 0
      min = 0
      sec = 0
    endelse
    formatstring = numformat + nameformat + $
                   ',3x,a3,4x,i4,1x,a3,1x,i2.2,4x,i2.2,a1,i2.2,a1,i2.2,3x,a3,3x' + $
                   ',i6,3x,f7.3,3x,' + sdev_format + ',3x,f6.1,5x,i4)'
    print,n+1,tname,pol,year,month[mon-1],day,hour,':',min,':',sec,timezone, $
          ndata,dfreq,sdev,phase,group,format=formatstring
  endfor
  n_so_far = n2 + 1
  if n_so_far lt nstack1 then begin
    print,' '
    pause,prompt='Press Enter to continue or q to quit: ',reply=userreply
  endif
endwhile

print,' '
   
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showstacki,page=page,first=nfirst,help=help

; Display all images in the image stack
 
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

maxnamelen = 16
month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

if n_params() gt 0 or keyword_set(help) then begin
  print,'showstacki[,page=n][,first=nfirst][,/help]'
  print,' '
  print,'n is the number of lines to display per page (default = all lines at once)'
  print,' '
  print,'nfirst is the first stacki image to display, counting from 1 (default = 1)'
  print,' '
  return
endif

if n_elements(page) eq 0 then page = nstacki
if n_elements(nfirst) eq 0 then nfirst = 1L

; Check that there are some images in the image stack to display

if nstacki eq 0 then begin
  print,'The image stack is empty'
  return
endif else if nstacki lt nfirst then begin
  print,'There are only ',nstacki,' images in the image stack',format='(a,i0,a)'
  return
endif

; Format the header line

if nstacki lt 1000 then begin
  numheadformat = '(a3,3x,'
  numformat = '(i3,3x,'
endif else begin
  numheadformat = '(a4,2x,'
  numformat = '(i4,2x,'
endelse
headformatstring = numheadformat + 'a,' + string((maxnamelen-4),format='(i0)') + $
                   'x,3x,a,7x,a,6x,a,2x,a,1x,a,4x,a,4x,a,3x,a,3x,a,2x,a)'

n_so_far = nfirst - 1
userreply = ''

while n_so_far lt nstacki and userreply ne 'q' and userreply ne 'Q' do begin

  ; Print the header line

  print,' '
  print,'Num','Name','Pol','Date and time','Zone','cols','rows','dfreq', $
        'baud','ddel','phase','Group',format=headformatstring

  ; Display information on each image

  n1 = n_so_far
  n2 = (n_so_far + page - 1) < (nstacki - 1)
  for n=n1,n2 do begin
    tname = (*stacki[n]).tname
    namelen = strlen(tname)
    if namelen lt maxnamelen then begin
      nameformat = 'a,' + string((maxnamelen-namelen),format='(i0)') + 'x'
    endif else begin
      nameformat = 'a' + string(maxnamelen,format='(i0)')
      print,"showstacki: Target name '",tname,"' will be truncated at ",maxnamelen, $
            " characters",format='(3a,i0,a)'
    endelse
    group = (*stacki[n]).group
    width = (*stacki[n]).width
    height = (*stacki[n]).height
    pol = (*stacki[n]).pol
    extratags = (*stacki[n]).extratags
    tznum = where(extratags.name eq 'timezone', count)
    if count gt 0 then begin
      tznum = tznum[0]
      timezone = extratags[tznum].value
    endif else begin
      timezone = '   '
    endelse
    jdmean_num = where(extratags.name eq 'jdmean', count)
    if count gt 0 then begin
      jdmean_num = jdmean_num[0]
      jdmean = double(extratags[jdmean_num].value)
    endif else begin
      jdstart_num = where(extratags.name eq 'jdstart', count1)
      jdend_num = where(extratags.name eq 'jdend', count2)
      if count1 gt 0 and count2 gt 0 then begin
        jdstart_num = jdstart_num[0]
        jdend_num = jdend_num[0]
        jdstart = double(extratags[jdstart_num].value)
        jdend = double(extratags[jdend_num].value)
        jdmean = (jdstart + jdend)/2D0
      endif else begin
        jdmean = 0D0
      endelse
    endelse
    caldat_roundsec,jdmean,mon,day,year,hour,min,sec
    if year eq -4713 then year = 0  ;  in case Julian date was set to 0
    dfreq_num = where(extratags.name eq 'fres', count)
    dfreq = (count gt 0) ? float(extratags[dfreq_num[0]].value) : 0.0
    phase_num = where(extratags.name eq 'phase', count)
    if count gt 0 then begin
      phase = string(float(extratags[phase_num[0]].value), format='(f5.1)')
    endif else begin
      phase = '     '
    endelse
    baud_num = where(extratags.name eq 'baudlen', count)
    baud = (count gt 0) ? float(extratags[baud_num[0]].value) : ''
    delayunit_num = where(extratags.name eq 'delayunit', count)
    if count gt 0 then begin
      delayunit = float(extratags[delayunit_num[0]].value)
    endif else if notnull(baud) then begin
      spb_num = where(extratags.name eq 'samples_per_baud', count)
      spb = (count gt 0) ? long(extratags[spb_num[0]].value) : 1L
      rows_per_baud_num = where(extratags.name eq 'rows_per_baud', count)
      rows_per_baud = (count gt 0) ? long(extratags[rows_per_baud_num[0]].value) : spb
      delayunit = baud/rows_per_baud
    endif else begin
      delayunit = -9.99
    endelse
    if isnull(baud) then baud = -9.99
    formatstring = numformat + nameformat + $
                   ',3x,a3,3x,i4,1x,a3,1x,i2.2,2x,i2.2,a1,i2.2,a1,i2.2,3x,a3,2x' + $
                   ',i4,1x,i4,2x,f7.3,2x,f6.2,1x,f6.2,3x,a5,1x,i4)'
    print,n+1,tname,pol,year,month[mon-1],day,hour,':',min,':',sec,  $
          timezone,width,height,dfreq,baud,delayunit,phase,group,    $
          format=formatstring
  endfor
  n_so_far = n2 + 1
  if n_so_far lt nstacki then begin
    print,' '
    pause,prompt='Press Enter to continue or q to quit: ',reply=userreply
  endif
endwhile

print,' '
   
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro resetstack,help=help

; Define/reset the pair stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'resetstack[,/help]'
  return
endif

if n_elements(nstack) eq 0 then begin
  nstack = 0L
endif else begin
  if nstack gt 0 then begin
    ptr_free,stack
    nstack = 0L
  endif
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro resetstack1,help=help

; Define/reset the single-channel stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'resetstack1[,/help]'
  return
endif

if n_elements(nstack1) eq 0 then begin
  nstack1 = 0L
endif else begin
  if nstack1 gt 0 then begin
    ptr_free,stack1
    nstack1 = 0L
  endif
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro resetstacki,help=help

; Define/reset the image stack

common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'resetstacki[,/help]'
  return
endif

if n_elements(nstacki) eq 0 then begin
  nstacki = 0L
endif else begin
  if nstacki gt 0 then begin
    ptr_free,stacki
    nstacki = 0L
  endif
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; MAIN PROGRAM which will be run when this file is first read in
; (i.e., at the command ".r general")
;
; Define the loaded image, the loaded pair, the loaded single-channel
; spectrum, the image stack, the pair stack, and the single-channel stack;
; set the tag variable types and initial display format; and set some
; plotting parameters

unload
unload1
unloadi
resetstack
resetstack1
resetstacki
resettagIO,/silent
!x.style = 1
!y.style = 3

idl_version = float(!VERSION.release)
minimum_idl_version = 5.6
if idl_version ge minimum_idl_version then begin
  finish_compilation = 1
endif else begin
  print,' '
  print,'*************************** WARNING *****************************'
  print,' '
  print,'Some procedures in this package require IDL version ', $
        minimum_idl_version,' or higher',format='(a,f3.1,a)'
  print,' '
  print,'You are using version ',idl_version,format='(a,f3.1)'
  print,' '
  userchoice = ''   ; initialize as a string or else read will assume float
  userprompt = 'Continue anyway? (y/n) '
  read,userchoice,prompt=userprompt
  userchoice = strlowcase(strmid(strtrim(userchoice), 0, 1))
  finish_compilation = (userchoice ne 'n')
  print,' '
  print,'*****************************************************************'
  print,' '
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
