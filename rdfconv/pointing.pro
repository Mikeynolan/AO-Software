pro pointing, datafile=datafile, yymmdd=yymmdd, buffer=buffer, perr=perr, herr=herr, average=average, silent=silent, noplot=noplot, outfile=outfile, plotfile=plotfile, multi=multi, stack=st, _extra=ext, help=help

; For chris stack stuff. Note that this defines 'stack', so keyword needs
; to be something else - used st
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) then begin
   print, " "
   print, "Usage:  pointing[, datafile='datafile', yymmdd=yymmdd, buffer=#, perr=#, herr=#, /average, /silent, /noplot, /multi, /outfile, /plotfile, plotKeywords]"
   print, " "
   print, "Plots (1) pointing errors due to platform tilt and height, (2) platform rotation about two axes,"
   print, " and (3) platform height versus time"
   print, " "
   print, "No inputs plots all data for the current date, yymmdd = ###### plots all data for requested date,"
   print, " and datafile = 'filename' plots a subset of time bracketing the scans in the datafile.  Note"
   print, " that using a datafile will override use of the yymmdd keyword for a specific date.  Using"
   print, " /datafile, rather than datafile='blah', will default to /share/olda/datafile"
   print, " "
   print, "If given a datafile, prints interpolated pointing errors at start of receive as a rough midpoint"
   print, " for the experiment TX/RX cycle.  If also given /average, instead average over the pointing errors"
   print, " between the TX start and RX stop times; this is only recommended for RTTs of several minutes,"
   print, " e.g., Mercury and MBAs, interpolation is better for RTTs of seconds, i.e., less than the 2 minute"
   print, " cadence of the tiedown measurements"
   print, " "
   print, "/silent will suppress text output to screen except error messages"
   print, "/noplot will suppress graphical output to screen"
   print, "/multi will plot all three plots at once, ignoring it will plot them sequentially"
   print, " "
   print, "/outfile will output comma-separated text files of the plotted data named pointing_error.txt,"
   print, " pointing_rotation.txt, and pointing_height.txt; output on multiple days will be *#.txt"
   print, " "
   print, "/plotfile will output PostScript and PNG plots named pointing_error, pointing_rotation, and"
   print, "/stack will see if there are spectra on the stack with matching times and add the errors as extratqags if there are."
   print, " pointing_height as .ps and .png files; plots for multiple days will be *#.ps and *#.png"
   print, " "
   print, "Also accepts standard plotting keywords, e.g., xrange (do not recommend yrange due to different"
   print, " ranges on the vertical axis in each plot)"
   print, " "
   print, "buffer = # is the search buffer around receive times for interpolation (default = 10 min)"
   print, "perr   = # is the warning limit for pointing error (default = 40 arcsec)"
   print, "herr   = # is the warning limit for height error (default = lambda/4 = 0.0315 m)"
   print, " "
   return
endif

; extend user's IDL path to include Patrick's IDL files
; preferably this should already be in your idlstartup file

!path = expand_path('+/home/ptaylor/idl')+':'+!path

; Phil's wrapper program, pnterrrot, defaults to IDL v7.1 because it
;  is the default IDL version for bash, the shell the wrapper is written
;  for.  If you use a more recent version of IDL, you must set:
;    setenv IDL_DLM_PATH /pkg/pcidl/idl71/bin/bin.linux.x86_64/
;  (in tcsh; and a related command in bash) to use Phil's version of
;  pnterrrot or else IDL will barf on startup

; instead one can use this program and not have to worry about 
;  the IDL versions and IDL_DLM_PATH

; the following are necessary to load for the procedure to run
;  properly; however, putting them here means this cannot be a true
;  IDL procedure because procedures cannot compile other procedures...
;  these must all instead be in the user's idlstartup file or loaded
;  by hand at the IDL prompt!

!Quiet = 1
;@phil
;@rirawinit
;@geninit
;@lrinit
;@agcinit
;@tdinit
;!Quiet = 0

pmulti = !p.multi        ; store original multiplot setting
;xmargin = !x.margin
;ymargin = !y.margin
!x.margin = [10.0, 3.0]  ; use default margins (riwatch can do weird things)
!y.margin = [4.0, 2.0]

; set up default arrays 
; (care about the mid-cycle time because pointing matters for both TX and RX)

hours = [0, 24]   ; allow xrange to work even if datafile not used
datearr = []      ; dates of scans read in from datafile
timearr = []      ; times of all records from datafile
scanarr = []      ; scannumbers read in from datafile
scanrecarr = []   ; record numbers where scan number increments
interpRot = []    ; interpolated values for RX start times
interphght = []
interpExtr0 = []
interpExtr1 = []
interpavgh = []

; set up other default parameters

if not keyword_set(buffer) then buffer = 10.  ; look for data within buffer of receive start
if not keyword_set(perr) then perr = 40.      ; print * if pointing error exceeds 40 arcsec
if not keyword_set(herr) then herr = 0.0315   ; print * if height error exceeds lambda/4

; IF we're going to save values to the stack, preload the stack with NaNs so the spectra are consistent.
;
;if keyword_set(st) then begin
  ;for ii = 1, nstack do begin
    ;setextra, 'f', 'perrorrot', !values.d_nan, comment='Rotation pointing error due to tiedown settings [arcsec]', stack=ii
    ;setextra, 'f', 'perrorhgt', !values.d_nan, comment='Height pointing error due to tiedown settings [arcsec]', stack=ii
    ;setextra, 'f', 'focusoff', !values.d_nan, comment='Focus error [m]', stack=ii
    ;setextra, 'i', 'badcal', 0, comment='No pointing information', stack=ii
  ;endfor
;endif
; if a date is not given by user, grab current date

if not keyword_set(yymmdd) then begin 
   if not keyword_set(datafile) then print, "Usage:  pointing[, datafile='datafile', yymmdd=yymmdd, /help]"
   spawn, 'yymmdd', result   ; get current date
   yymmdd = long(result[0])  ; switch from string to numerical value
endif

; break down date into year, month, and day

yy = yymmdd/10000l               ; two-digit year
yr = yy + 2000l                  ; four-digit year
mon = (yymmdd/100) mod 100l      ; month number
day = yymmdd mod 100l            ; day of month
dayno = dmtodayno(day, mon, yr)  ; day of year

if keyword_set(datafile) then begin  ; if a datafile is given

   ; if /datafile (not a quoted string), default to /share/olda/datafile

   if (size(datafile, /type) ne 7) then datafile = '/share/olda/datafile'

   ; open and try to read the datafile

   openr, unit, datafile, /get_lun  ; open datafile
   n = riget(unit, b, /complex)     ; read in first record to struct b
   if (n eq 0) then begin
      print, " "
      print, " ERROR:  No readable RI data found in datafile!"
      print, " "
      return
   endif
   point_lun, unit, 0  ; go back to beginning of file

   yyyydoy = b.h.std.date        ; date given as four-digit year (yyyy) and day of year (doy)
   secofday = b.h.std.time       ; time given as seconds from AST midnight
   scannum = b.h.std.scannumber  ; scan number of record 

   yr = yyyydoy/1000l               ; four-digit year
   yy = yr-(yr/1000l)*1000l         ; two-digit year
   dayno = yyyydoy-yr*1000l         ; day of year
   daymon = daynotodm(dayno, yr)    ; day and month in an array
   mon = daymon[1]                  ; month number
   day = daymon[0]                  ; day of month
   yymmdd = yy*10000l+mon*100l+day  ; yymmdd format 

   datearr = [yymmdd]    ; keep track of dates, times, and scans
   timearr = [secofday]  ; time for each record
   scanarr = [scannum]   ; scan number
   scanrecarr = [1]      ; first record of each scan
 
endif

; set up date arrays in case datafile spans multiple AST dates

dates = [yymmdd]  ; yymmdd found either by default, command line, or from datafile
daynos = [dayno]  ; day number by similar means
yrs = [yr]        ; year

; with date established, start reading through datafile again

rec = 0  ; number of records read in

if keyword_set(datafile) then begin  ; if a datafile is given

   minhour = 25.  ; absurd minimum timestamp in hours
   maxhour = 0.   ; absurd maximum timestamp in hours

   repeat begin  ; loop through datafile record by record

      n = riget(unit, b, /complex)  ; read in a record
      rec = rec + 1                 ; increment record

      scannum = b.h.std.scannumber  ; scan number of record

      ; record time of record

      secofday = b.h.std.time        
      timearr = [timearr, secofday]

      ; check date for midnight rollover or concatenated dates

      yyyydoy = b.h.std.date         
      yr = yyyydoy/1000l             
      yy = yr-(yr/1000l)*1000l       
      dayno = yyyydoy-yr*1000l       
      daymon = daynotodm(dayno, yr)  
      mon = daymon[1]
      day = daymon[0]                
      yymmdd = yy*10000l+mon*100l+day
      hour = secofday/3600.           

      ; if the datafile has multiple dates, either from combining multiple
      ;  tracks or crossing AST midnight, split the dates up or else the pointing
      ;  files that increment at midnight could get confused

      if (yymmdd ne dates[n_elements(dates)-1]) then begin  ; new date
         dates = [dates, yymmdd]   ; add yymmdd to array for use later
         daynos = [daynos, dayno]  ; add dayno to array
         yrs = [yrs, yr]           ; add yr to array
         hours = [hours, minhour, maxhour]  ; array of min/max hours per date
         minhour = 25.             ; reset min value for new date
         maxhour = 0.              ; reset max value for new date
      endif

      if (hour lt minhour) then minhour = hour  ; new minimum time for date
      if (hour gt maxhour) then maxhour = hour  ; new maximum time for date

      ; check for new scan number

      if (scannum ne scanarr[n_elements(scanarr)-1]) then begin  ; record begins new scan
         scanarr = [scanarr, scannum]    ; add new scan number to array
         scanrecarr = [scanrecarr, rec]  ; add first record of new scan to array
         datearr = [datearr, yymmdd]     ; add date of new scan to array
      endif

   endrep until eof(unit)  ; keep looping until end of file reached

   hours = [hours, minhour, maxhour]  ; must put something in the hours array

   ; since this is a datafile with specific hours, delete the placeholders from
   ;  the first two positions

   hours = hours[2:n_elements(hours)-1] 

   ; account for last scan (EOF prevents new scan number from triggering save)

   scanrecarr = [scanrecarr, rec]

   close, unit     ; close datafile
   free_lun, unit  ; free file unit

endif

; with dates collected, collect the relevant data

for i=0, n_elements(dates)-1 do begin  ; loop through dates in case of date rollover

   if keyword_set(multi) then !p.multi = [0, 1, 3] else !p.multi = 0  ; multiple plots
   print, " "

   ; get distomat data from laser rangefinders for requested date

   n = lrpcinp(dates[i], lr, /ext, /quiet)

   if (n le 0) then begin  ; bail out if no distomat data found
      print, " "
      print," ERROR:  No distomat data found for date (yymmdd):  ", strtrim(string(dates[i],format='(I)'),1)
      print, " "
      return
   endif

   ; check that data in the distomat file is for the requested date

   ii = where(long(lr.date) eq daynos[i], cnt)  ; grab indices of matching data
   if cnt eq 0 then begin                       ; no data found for requested date
      print, " "
      print, " ERROR:  Distomat file has no data for date (yymmdd):  ", strtrim(string(dates[i],format='(I)'),1)
      print, " "
      return
   endif
      
   if cnt ne n then begin  ; if some of the data is from the wrong date
      lr = lr[ii]          ;  only use the data from the correct date
      n = cnt 
   endif

   ; everything below is necessary to make the plots
   ; Phil's comment is that "if current month, need to get agc, td separately"

   if (lr[0].az eq -1) then begin
      n = agcinpday(dates[i], b)  ; get azimuth encoder data
      az = b.cb.pos[0]            ; azimuth data
      za = b.cb.pos[1]            ; zenith angle data
      hragc = b.cb.time
      hragc[0] = hragc[1]
      hragc /= 3600.              ; change seconds to hours
      lrdate = lr.date            ; array of fractional day of year
      lrdate0 = long(lrdate[0])   ; first timestamp
      ii = where((long(lrdate) ge (lrdate0 + 1)), cnt)  ; check for date rollover
      hrlr = (lrdate mod 1d)*24.           ; change days to fraction of days then to hours
      if cnt gt 0 then hrlr[ii] += 24.     ; add a day to those after rollover (due to mod)
      lr.az = interpol(az, hragc, hrlr)    ; interpolate the azimuth data to distomat timestamps
      lr.zagr = interpol(za, hragc, hrlr)  ; interpolate the zenith angle data
   endif

   ; call platform pointing model with year and distomat data
   ; returns results in the rotI structure

   n = platrotmodel(yrs[i], lr, rotI)
   ii = where(rotI.ok eq 1, cnt)  ; indices with valid data
   if cnt eq 0 then begin 
      print, " "
      print, " ERROR:  Cannot find laser rangefinder data for all 6 distomats.  Plots cannot be made!"
      print, " "
      return
   endif

   ; grab min, max, spread of pointing errors due to platform tilt and height
   ; determine min and max only between xrange limits

   ; have to be careful to re-define the arrays after each where
   ;  statement to ensure proper mapping of the times and values

   plothr = rotI[ii].hr            ; only keep indices of valid data (time)
   plotRot = rotI[ii].pntErrRot    ; valid errors from platform rotation in radians
   plothght = rotI[ii].pntErrhght  ; valid errors from platform height in radians

   iii = where(plothr ge (hours[2*i]-1))  ; only keep times after minimum request - 1 
   plothr = plothr[iii]                   ; valid times after time request begins
   plotRot = plotRot[iii]                 ; remaining errors from platform rotation
   plothght = plothght[iii]               ; remaining errors from platform height

   iiii = where(plothr le (hours[2*(i+1)-1]+1))  ; only keep times before maximum request + 1
   plothr = plothr[iiii]                  ; valid times before time request ends
   plotRot = plotRot[iiii]                ; remaining errors between time requests
   plothght = plothght[iiii]              ; remaining errors between time requests

   plotExtr0 = rotI[ii].rotExtr[0]
   plotExtr1 = rotI[ii].rotExtr[1]
   plotExtr0 = plotExtr0[iii]
   plotExtr1 = plotExtr1[iii]
   plotExtr0 = plotExtr0[iiii]
   plotExtr1 = plotExtr1[iiii]

   plotavgh = rotI[ii].avgh  ; average platform height on date and time range requested
   plotavgh = plotavgh[iii]
   plotavgh = plotavgh[iiii]
   
   ftom = 12.*2.54/100.  ; convert feet to meters
   fh = 1256.22*ftom     ; convert nominal platform height from feet to meters

   iiiii = where(abs(plotavgh-fh) lt 1., cnt)  ; note where the height difference is less than 1 m
                                               ;  (to avoid illogical values/zeros)
   plotavgh = plotavgh[iiiii]
   plothrs = plothr[iiiii]

   ; do interpolation or averaging for each scan

   first = 1  ; first time through, print a header comment

   for j=0, n_elements(datearr)-1 do begin  ; loop through dates of scans
      warning = ' '  ; reset warning symbol for questionable pointing

      ; for interpolation, use RX start time from timearr[scanrecarr[j]]
      ; for /average, want estimate tiedown measurement between TX start and RX stop
      ;  RX start = timearr[scanrecarr[j]]
      ;  RX stop  = timearr[scanrecarr[j+1]-1]
      ;  TX start = RX start - (RX stop - RX start + 8)
      ;  any calculated time must be scaled by 3600 from sec of day to decimal hours

      rxstart = timearr[scanrecarr[j]]
      rxstop = timearr[scanrecarr[j+1]-1]
      if (rxstop lt rxstart) then rxstop = 86399
      txstart = rxstart - (rxstop-rxstart + 8)
      if (txstart lt 0) then txstart = 1

      if keyword_set(average) then begin  ; use TX start/RX stop in hours
         mintime = txstart/3600.
         maxtime = rxstop/3600.
      endif else begin                    ; use RX start in hours
         mintime = rxstart/3600.
         maxtime = rxstart/3600.
      endelse

      if (datearr[j] eq dates[i]) then begin  ; if on current date
         hr1 = 0       ; reset time limits and positions bracketing mid-cycle time
         hr2 = 24
         kmin = 1e30
         kmax = 1e-30
         for k=0, n_elements(plothr)-1 do begin   ; loop through valid times
            if (plothr[k] le mintime) then begin  ; bracket from below
               hr1 = plothr[k]
               kmin = k
            endif
            if (plothr[k] ge maxtime) then begin  ; bracket from above
               hr2 = plothr[k]
               kmax = k
               break  ; once bracketed, break out of the loop
            endif
         endfor ; should now have bracketed scan

         ; if kmin/kmax never changed or the closest times are more than 'buffer'
         ;  from the requested times, warn user because there likely is not enough
         ;  useful data from tiedowns

         if ((kmin eq 1e30) or (kmax eq 1e-30) or (abs(hr1-mintime) gt buffer/60.) or (abs(hr2-maxtime) gt buffer/60.)) then begin
            interpRot = [interpRot, -999]      ; set to bogus values as placeholders if printed
            interphght = [interphght, -999]
            interpExtr0 = [interpExtr0, -999]
            interpExtr1 = [interpExtr1, -999]
            interpavgh = [interpavgh, -999]
            if not keyword_set(silent) then begin
               print, " WARNING:  No valid tiedown data within ", strtrim(string(buffer,format='(I)'),1), " minutes of bracketing time for scan:  ", strtrim(string(scanarr[j]),1)
            endif
         endif else begin
            if ((not keyword_set(silent)) and (first eq 1)) then begin
               print, "# Scan Number      ErrTilt       ErrHeight      RotWest       RotNorth    Height    * = ", strtrim(string(perr,format='(I)'),1), "+ arcsec or ", strtrim(string(herr,format='(F6.4)'),1)," m"
               first = 0
            endif

;           once reasonable bracketing found, take average of tiedown measurements

;            print, kmin, kmax, hr1, hr2

            if keyword_set(average) then begin
               interpRot = [interpRot, total(plotRot[kmin:kmax])/(kmax-kmin+1)*3600.]
               interphght = [interphght, total(plothght[kmin:kmax])/(kmax-kmin+1)*3600.]
               interpExtr0 = [interpExtr0, total(plotExtr0[kmin:kmax])/(kmax-kmin+1)*3600.]
               interpExtr1 = [interpExtr1, total(plotExtr1[kmin:kmax])/(kmax-kmin+1)*3600.]
               interpavgh = [interpavgh, total(plotavgh[kmin:kmax])/(kmax-kmin+1)-fh]
            endif else begin  ; interpolate
               interpRot = [interpRot, (plotRot[kmin] + (plotRot[kmax]-plotRot[kmin])/(hr2-hr1)*(rxstart/3600.-hr1))*3600.]
               interphght = [interphght, (plothght[kmin] + (plothght[kmax]-plothght[kmin])/(hr2-hr1)*(rxstart/3600.-hr1))*3600.]
               interpExtr0 = [interpExtr0, (plotExtr0[kmin] + (plotExtr0[kmax]-plotExtr0[kmin])/(hr2-hr1)*(rxstart/3600.-hr1))*3600.]
               interpExtr1 = [interpExtr1, (plotExtr1[kmin] + (plotExtr1[kmax]-plotExtr1[kmin])/(hr2-hr1)*(rxstart/3600.-hr1))*3600.]
               interpavgh = [interpavgh, plotavgh[kmin] + (plotavgh[kmax]-plotavgh[kmin])/(hr2-hr1)*(rxstart/3600.-hr1)-fh]
            endelse

;           add a warning symbol for large excursions in pointing and platform height

            if ((interpRot[j] ge perr) or (interphght[j] ge perr) or (abs(interpavgh[j]) gt herr)) then warning = '*'

            if not keyword_set(silent) then begin
               print, "Scan ", strtrim(string(scanarr[j],format='(I)'),1), ": ", string(interpRot[j],format='(F6.1)'), " arcsec ", string(interphght[j],format='(F6.1)'), " arcsec ", string(interpExtr0[j],format='(F6.1)'), " arcsec ", string(interpExtr1[j],format='(F6.1)'), " arcsec ", string(interpavgh[j],format='(F7.3)'), " m   ", warning
            endif
         endelse
      endif
   endfor

   if not keyword_set(silent) then print, " "

   ; plot pointing errors due to platform tilt and height
   ; limit xrange to one hour before and after records in datafile

   if not keyword_set(noplot) then begin

      min = min([plotRot, plothght])*3600. ; minimum error due to platform rotation/height
      max = max([plotRot, plothght])*3600. ; maximum error due to platform rotation/height
      dif = max-min                        ; difference in min/max for plot ranges
  
      ; charsize is a function of !p.multi:  if three plots are made at
      ;  once, charsize is increased, otherwise it is 1.5

      plot, plothr, plotRot*3600., xrange=[(hours[2*i]-1)>0, (hours[2*(i+1)-1]+1)<24], /xstyle, yrange=[min-0.1*dif, max+0.1*dif], charsize=1.5*(!p.multi[2]/3.+1.), xtitle='Time in AST on ' + strtrim(string(dates[i]),1) + ' [h]', ytitle='Pointing Error [arcsec]', title='Pointing Error from Platform Tilt and Height', _extra=ext
      cgoplot, plothr, plotRot*3600., color='green', psym=-1 
      cgoplot, plothr, plothght*3600., color='red', psym=-1

      ; overplot nominal pointing and 1 arminute error

      oplot, [!X.CRANGE[0],!X.CRANGE[1]], [0,0], linestyle=2
      cgoplot, [!X.CRANGE[0],!X.CRANGE[1]], [60,60], linestyle=2, color='red'
      cgoplot, [!X.CRANGE[0],!X.CRANGE[1]], [-60,-60], linestyle=2, color='red'

      ; plot first and last timestamps of records in datafile

      if keyword_set(datafile) then begin
         oplot, [hours[2*i], hours[2*i]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2
         oplot, [hours[2*(i+1)-1], hours[2*(i+1)-1]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2
      endif

      ; include legend for pointing errors

      al_legend, ['Platform Tilt','Platform Height'], linestyle=[0,0], charsize=2, box=0, /top, /left, thick=2, linsize=0.2, colors=['green', 'red'], textcolors=['white', 'white']

      if not keyword_set(multi) then pause ; pause if plotting one at a time

      ; grab min, max, spread of platform rotation
      ; determine min and max only between xrange limits

      min = min([plotExtr0,plotExtr1])*3600.
      max = max([plotExtr0,plotExtr1])*3600.
      dif = max-min

      ; plot platform rotation about two axes

      plot, plothr, plotExtr0*3600, xrange=[(hours[2*i]-1)>0, (hours[2*(i+1)-1]+1)<24], /xstyle, yrange=[min-0.1*dif, max+0.1*dif], charsize=1.5*(!p.multi[2]/3.+1), xtitle='Time in AST on ' + strtrim(string(dates[i]),1) + ' [h]', ytitle='Platform Tilt [arcsec]', title='Platform Tilt (wrt model)', _extra=ext
      cgoplot, plothr, plotExtr0*3600, color='green', psym=-1
      cgoplot, plothr, plotExtr1*3600, color='red', psym=-1

      ; overplot nominal pointing

      oplot, [!X.CRANGE[0],!X.CRANGE[1]], [0,0], linestyle=2

      ; plot first and last timestamps of records in datafile
      
      if keyword_set(datafile) then begin
         oplot, [hours[2*i], hours[2*i]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2
         oplot, [hours[2*(i+1)-1], hours[2*(i+1)-1]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2
      endif

      ; include legend for platform rotations

      al_legend, ['X (west) Rotation', 'Y (north) Rotation'], linestyle=[0,0], charsize=2, box=0, /top, /left, thick=2, linsize=0.2, colors=['green', 'red'], textcolors=['white', 'white']

      if not keyword_set(multi) then pause ; pause if plotting one at a time

      ; grab min, max, and diff of platform height
      ; determine min and max only between xrange limits

      min = min([plotavgh-fh, 0]) ; get minimum deviation from nominal height
      max = max([plotavgh-fh, 0]) ; get maximum deviation from nominal height
      dif = max - min

      ; plot platform height

      plot, plothrs, plotavgh-fh, xrange=[(hours[2*i]-1)>0, (hours[2*(i+1)-1]+1)<24], /xstyle, yrange=[min-0.1*dif, max+0.1*dif], charsize=1.5*(!p.multi[2]/3.+1), xtitle='Time in AST on ' + strtrim(string(dates[i]),1) + ' [h]', ytitle='Offset from Nominal Platform Height [m]', title='Platform Height (Focus)', _extra=ext
      cgoplot, plothrs, plotavgh-fh, color='green', psym=-1

      ; overplot nominal and half-wavelength (S band) offsets of platform height

      oplot, [!X.CRANGE[0],!X.CRANGE[1]], [0,0], linestyle=2
      cgoplot, [!X.CRANGE[0],!X.CRANGE[1]], [-0.063,-0.063], linestyle=2, color='red'
      cgoplot, [!X.CRANGE[0],!X.CRANGE[1]], [0.063,0.063], linestyle=2, color='red'

      ; plot first and last timestamps of records in datafile

      if keyword_set(datafile) then begin
         oplot, [hours[2*i], hours[2*i]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2
         oplot, [hours[2*(i+1)-1], hours[2*(i+1)-1]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2
      endif

      ; include legend for platform height

      al_legend, ['Nominal Platform Height', 'Half-Wavelength Offset (0.063 m)'], linestyle=[2,2], charsize=2, box=0, /top, /left, thick=2, linsize=0.2, colors=['white', 'red'], textcolors=['white', 'white']

      if not keyword_set(multi) then pause ; pause if plotting one at a time

   endif

   ; output to PostScript files if requested
   
   if keyword_set(plotfile) then begin

      !p.multi = 0  ; always output plots to individual files
         
      ; plot pointing error to file

      thisdevice = !D.Name
      set_plot, 'PS'

      if (n_elements(dates) gt 1) then plotstem = 'pointing_error'+strtrim(string(i+1),1) else plotstem = 'pointing_error'
;      psobject = obj_new("fsc_psconfig", filename=plotstem, defaultsetup="Color (Landscape)")
;      pskeywords = psobject->getkeywords()
;      device, _extra=pskeywords
      device, file=plotstem+'.ps', /landscape, /color, bits=8

      min = min([plotRot,plothght])*3600.
      max = max([plotRot,plothght])*3600.
      dif = max-min
      plot, plothr, plotRot*3600., xrange=[(hours[2*i]-1)>0, (hours[2*(i+1)-1]+1)<24], /xstyle, yrange=[min-0.1*dif, max+0.1*dif], charsize=1.5, xtitle='Time in AST on ' + strtrim(string(dates[i]),1) + ' [h]', ytitle='Pointing Error [arcsec]', title='Pointing Error from Platform Tilt and Height', thick=6, xthick=6, ythick=6, charthick=6, font=0, _extra=ext
      cgoplot, plothr, plotRot*3600., color='green', psym=-1, thick=6
      cgoplot, plothr, plothght*3600., color='red', psym=-1, thick=6

      oplot, [!X.CRANGE[0],!X.CRANGE[1]], [0,0], linestyle=2, thick=6
      cgoplot, [!X.CRANGE[0],!X.CRANGE[1]], [60,60], linestyle=2, color='red', thick=6
      cgoplot, [!X.CRANGE[0],!X.CRANGE[1]], [-60,-60], linestyle=2, color='red', thick=6
      if keyword_set(datafile) then begin
         oplot, [hours[2*i], hours[2*i]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2, thick=6
         oplot, [hours[2*(i+1)-1], hours[2*(i+1)-1]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2, thick=6
      endif
      al_legend, ['Platform Tilt','Platform Height'], linestyle=[0,0], charsize=1.5, box=0, /top, /left, thick=6, charthick=4, linsize=0.2, colors=['green', 'red'], textcolors=['black', 'black'], font=0
      sharpcorners, thick=8
      device, /close_file

      convertstring = 'convert -rotate 270 -flatten '+plotstem+'.ps '+plotstem+'.png'
      spawn, convertstring

      ; plot platform rotation to file

      if (n_elements(dates) gt 1) then plotstem = 'pointing_rotation'+strtrim(string(i+1),1) else plotstem = 'pointing_rotation'
;      psobject = obj_new("fsc_psconfig", filename=plotstem, defaultsetup="Color (Landscape)")
;      pskeywords = psobject->getkeywords()
;      device, _extra=pskeywords
      device, file=plotstem+'.ps', /landscape, /color, bits=8

      min = min([plotExtr0,plotExtr1])*3600.
      max = max([plotExtr0,plotExtr1])*3600.
      dif = max-min
      plot, plothr, plotExtr0*3600, xrange=[(hours[2*i]-1)>0, (hours[2*(i+1)-1]+1)<24], /xstyle, yrange=[min-0.1*dif, max+0.1*dif], charsize=1.5, xtitle='Time in AST on ' + strtrim(string(dates[i]),1) + ' [h]', ytitle='Platform Tilt [arcsec]', title='Platform Tilt (wrt model)', thick=6, xthick=6, ythick=6, charthick=6, font=0, _extra=ext
      cgoplot, plothr, plotExtr0*3600, color='green', psym=-1, thick=6
      cgoplot, plothr, plotExtr1*3600, color='red', psym=-1, thick=6
      oplot, [!X.CRANGE[0],!X.CRANGE[1]], [0,0], linestyle=2, thick=6
      if keyword_set(datafile) then begin
         oplot, [hours[2*i], hours[2*i]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2, thick=6
         oplot, [hours[2*(i+1)-1], hours[2*(i+1)-1]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2, thick=6
      endif
      al_legend, ['X (west) Rotation', 'Y (north) Rotation'], linestyle=[0,0], charsize=1.5, box=0, /top, /left, thick=6, charthick=4, linsize=0.2, colors=['green', 'red'], textcolors=['black', 'black'], font=0
      sharpcorners, thick=8
      device, /close_file

      convertstring = 'convert -rotate 270 -flatten '+plotstem+'.ps '+plotstem+'.png'
      spawn, convertstring

      ; plot platform height to file

      if (n_elements(dates) gt 1) then plotstem = 'pointing_height'+strtrim(string(i+1),1) else plotstem = 'pointing_height'
;      psobject = obj_new("fsc_psconfig", filename=plotstem, defaultsetup="Color (Landscape)")
;      pskeywords = psobject->getkeywords()
;      device, _extra=pskeywords
      device, file=plotstem+'.ps', /landscape, /color, bits=8

      min = min([plotavgh-fh, 0])
      max = max([plotavgh-fh, 0])
      dif = max - min
      plot, plothrs, plotavgh-fh, xrange=[(hours[2*i]-1)>0, (hours[2*(i+1)-1]+1)<24], /xstyle, yrange=[min-0.1*dif, max+0.1*dif], charsize=1.5, xtitle='Time in AST on ' + strtrim(string(dates[i]),1) + ' [h]', ytitle='Offset from Nominal Platform Height [m]', title='Platform Height (Focus)', thick=6, xthick=6, ythick=6, charthick=6, font=0, _extra=ext
      cgoplot, plothrs, plotavgh-fh, color='green', psym=-1, thick=6
      oplot, [!X.CRANGE[0],!X.CRANGE[1]], [0,0], linestyle=2, thick=6
      cgoplot, [!X.CRANGE[0],!X.CRANGE[1]], [-0.063,-0.063], linestyle=2, color='red', thick=6
      cgoplot, [!X.CRANGE[0],!X.CRANGE[1]], [0.063,0.063], linestyle=2, color='red', thick=6
      if keyword_set(datafile) then begin
         oplot, [hours[2*i], hours[2*i]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2, thick=6
         oplot, [hours[2*(i+1)-1], hours[2*(i+1)-1]], [!Y.CRANGE[0],!Y.CRANGE[1]], linestyle=2, thick=6
      endif
      al_legend, ['Nominal Platform Height', 'Half-Wavelength Offset (0.063 m)'], linestyle=[2,2], charsize=1.5, box=0, /top, /left, thick=6, charthick=4, linsize=0.2, colors=['black', 'red'], textcolors=['black', 'black'], font=0
      sharpcorners, thick=8
      device, /close_file

      convertstring = 'convert -rotate 270 -flatten '+plotstem+'.ps '+plotstem+'.png'
      spawn, convertstring
      
      set_plot, thisdevice

   endif

   if keyword_set(outfile) then begin  ; CSV text output

      ; above was to create plots, now do similarly for CSV text output if asked
      ;  but this time do not use the +/- 1 hour buffer to the time request

      ; but sometimes the important values are not within the RX brackets
      ;  so print the same range as plotted; min/max can be limited

      printhr = rotI[ii].hr             ; only keep indices of valid data (time)
      printRot = rotI[ii].pntErrRot     ; valid errors from platform rotation in radians
      printhght = rotI[ii].pntErrhght   ; valid errors from platform height in radians
      printExtr0 = rotI[ii].rotExtr[0]  ; valid tilts from X rotation
      printExtr1 = rotI[ii].rotExtr[1]  ; valid tilts from Y rotation
      printavgh = rotI[ii].avgh         ; valid average platform heights

      iii = where(printhr ge hours[2*i], cnt1)  ; only keep data after start of request
      printhr = printhr[iii]                    ; if cnt1 = 0, there was no valid data
      printRot = printRot[iii]
      printhght = printhght[iii]
      printExtr0 = printExtr0[iii]
      printExtr1 = printExtr1[iii]
      printavgh = printavgh[iii]
 
      cnt2 = 0  ; set counter for data found prior to end of request

      if (cnt1 gt 0) then begin  ; there is valid data after start of request
         iiii = where(printhr le hours[2*(i+1)-1], cnt2)  ; only keep times before end of request
         printhr = printhr[iiii]
         printRot = printRot[iiii]
         printhght = printhght[iiii]
         printExtr0 = printExtr0[iiii]
         printExtr1 = printExtr1[iiii]
         printavgh = printavgh[iiii]

         iiiii = where(abs(printavgh-fh) lt 1.)  ; reject any height differences > 1 meter
         printavgh = printavgh[iiiii]
         printhrs = printhr[iiiii]               ; printhrs never used
     
         min1 = min([printRot, printhght])*3600.     ; min/max of each dataset
         max1 = max([printRot, printhght])*3600.
         min2 = min([printExtr0, printExtr1])*3600.
         max2 = max([printExtr0, printExtr1])*3600.
         min3 = min(printavgh-fh)
         max3 = max(printavgh-fh)
      endif
  
      ; create filenames for output text file(s)

      if (n_elements(dates) gt 1) then begin  ; multiple dates means multiple files
         printstem = 'pointing_error'+strtrim(string(i+1),1)
      endif else begin  ; single date and file
         printstem = 'pointing_error'
      endelse

      printname = printstem+'.txt'

      ; open text output file

      get_lun, unit
      openw, unit, printname

      ; print header for text output file
      
      printf, unit, "# "
      printf, unit, "# Pointing error due to platform tilt and height, platform rotation about X (west)"
      printf, unit, "#  and Y (north) axes, and relative platform height during time request based on"
      printf, unit, "#  tiedown data and telescope pointing model"
      printf, unit, "# "
      printf, unit, "# A value of 60 arcsec or more roughly exceeds the half-power half-beam width"
      printf, unit, "# A value of +/- 0.063 m is a half wavelength at S-band and is out of focus"

      if keyword_set(datafile) then begin
         printf, unit, "# "
         printf, unit, "# Datafile = ", datafile
         if keyword_set(average) then begin
            method = "# Average values for full-cycle (TX start to RX stop) times of scans with 9 columns:"
         endif else begin
            method = "# Interpolated values for mid-cycle (RX start) times of scans with 9 columns:"
         endelse
         printf, unit, "# "
         printf, unit, method
         printf, unit, "# "
         printf, unit, "# Scan number, Date (yyyy-mm-dd), Time (hh:mm:ss, AST), Time (h, AST), Tilt error (arcsec),"
         printf, unit, "#  Height error (arcsec), X rotation (arcsec), Y rotation (arcsec), Relative platform height (m)"
         printf, unit, "# "
      endif
         
      for j=0, n_elements(datearr)-1 do begin ; loop through dates of scans
         if (datearr[j] eq dates[i]) then begin ; if on current date

            ; gather date in yyyy-mm-dd format, time in hh:mm, and decimal hours in AST for output

            if (datearr[j] lt 210000) then century = 2000 else century = 1900
            printdate = strtrim(string(century+strmid(strtrim(datearr[j],1),0,2),format='(I)'),1)+'-'+strmid(strtrim(datearr[j],1),2,2)+'-'+strmid(strtrim(datearr[j],1),4,2)  ; yyyy-mm-dd

            dechour = timearr[scanrecarr[j]]/3600.  ; decimal hour
            printhour = strtrim(floor(dechour),1)
            if (printhour lt 10) then printhour = '0'+strtrim(string(printhour),1)  ; zero pad hours

            printmin = strtrim(floor((dechour-floor(dechour))*60.),1)  ; integer minutes
            if (printmin lt 10) then printmin = '0'+strtrim(string(printmin),1)  ; zero pad minutes

            printsec = strtrim(round(((dechour-floor(dechour))*60.-floor((dechour-floor(dechour))*60.))*60.),1)  ; integer seconds
            if (printsec lt 10) then printsec = '0'+strtrim(string(printsec),1)  ; zero pad seconds

            printtime = printhour+":"+printmin+":"+printsec  ; hh:mm:ss

            string = strtrim(string(scanarr[j],format='(I)'),1)+","+printdate+","+printtime+","+strtrim(string(timearr[scanrecarr[j]]/3600.,format='(D)'),1)+","+strtrim(string(interpRot[j],format='(F7.2)'),1)+","+strtrim(string(interphght[j],format='(F7.2)'),1)+","+strtrim(string(interpExtr0[j],format='(F7.2)'),1)+","+strtrim(string(interpExtr1[j],format='(F7.2)'),1)+","+strtrim(string(interpavgh[j],format='(F8.3)'),1)
            printf, unit, string
         endif
; MCN If /stack keyword is used and there are data loaded, add extratags with info.
        if keyword_set(st) then begin
; See if the date is one of the ones on the stack. If so, add extratag
; Get JD
          myear = long(datearr[j]/10000)
          mmon = long((datearr[j] - myear * 10000) / 100)
          mday = long(datearr[j] mod 100)
          myear = myear + century
          mhr = long(dechour)
          mmn = long((dechour - mhr) * 60)
          msec = long(((dechour - mhr) * 60 - mmn) *60 + .5)
          mjd = julday(mmon, mday, myear, mhr, mmn,msec)
          mjd = mjd + 4.d0 / 24

          if nstack eq 0 then begin
            print,'ERROR in pointing with stack: The pair stack is empty'
          endif

          for ii = 1, nstack do begin
            jds = getextra('jdstart', stack=ii)
            if (mjd - jds lt 1.d-5) then begin
              perror = sqrt(interpRot[j]*interpRot[j]+interphght[j]+interphght[j])
              setextra, 'f', 'perror', perror, comment='pointing error due to tiedown settings [arcsec]', stack=ii
              setextra, 'f', 'focusoff', interpavgh[j], comment='Focus error [m]', stack=ii
              calflg = perror gt perr or abs(interpavgh[j]) gt herr
              comstr = strcompress("true if pointing error exceeds " + string(perr) + " or focus error is more than " + string(herr) + " m")
              setextra, 'f', 'badcal', long(calflg), comment=comstr, stack=ii
            endif
          endfor
        endif ; if st

      endfor

      printf, unit, "# "
      printf, unit, "# "
      printf, unit, "# Available data from tiedowns (during time request):"

      ; if the counters returned zero valid data, print "n/a" to output files and quit
      ;  else print the valid min/max values to the headers

      if ((cnt1 eq 0) or (cnt2 eq 0)) then begin

         printf, unit, "# "
         printf, unit, "# Minimum value (platform tilt/height) = n/a, no valid data found for time request!"
         printf, unit, "# Maximum value (platform tilt/height) = n/a"
         printf, unit, "# "
         printf, unit, "# Minimum value (platform rotation) = n/a, no valid data found for time request!"
         printf, unit, "# Maximum value (platform rotation) = n/a"
         printf, unit, "# "
         printf, unit, "# Minimum value (platform height) = n/a, no valid data found for time request!"
         printf, unit, "# Maximum value (platform height) = n/a"
         printf, unit, "# "

      endif else begin

         printf, unit, "# "
         printf, unit, "# Minimum value (platform tilt/height) = "+strtrim(string(min1, format='(F7.2)'),1)+" arcsec"
         printf, unit, "# Maximum value (platform tilt/height) = "+strtrim(string(max1, format='(F7.2)'),1)+" arcsec"
         printf, unit, "# "
         printf, unit, "# Minimum value (platform rotation) = "+strtrim(string(min2, format='(F7.2)'),1)+" arcsec"
         printf, unit, "# Maximum value (platform rotation) = "+strtrim(string(max2, format='(F7.2)'),1)+" arcsec"
         printf, unit, "# "
         printf, unit, "# Minimum value (platform height) = "+strtrim(string(min3, format='(F6.3)'),1)+" m"
         printf, unit, "# Maximum value (platform height) = "+strtrim(string(max3, format='(F6.3)'),1)+" m"
         printf, unit, "# "

         ; print valid lines to the output text files

         printf, unit, "# "
         printf, unit, "# Available data from tiedowns from time request +/- 1 h (6 columns):"
         printf, unit, "# "
         printf, unit, "# Time (h, AST), Tilt error (arcsec), Height error (arcsec), X rotation (arcsec),"
         printf, unit, "#  Y rotation (arcsec), Relative platform height (m)"
         printf, unit, "# "

         for n=0L, n_elements(plothr)-1 do begin
            string = strtrim(string(plothr[n], format='(D)'),1)+','+strtrim(string(plotRot[n]*3600., format='(F7.2)'),1)+','+strtrim(string(plothght[n]*3600., format='(F7.2)'),1)+','+strtrim(string(plotExtr0[n]*3600., format='(F7.2)'),1)+','+strtrim(string(plotExtr1[n]*3600., format='(F7.2)'),1)+','+strtrim(string(plotavgh[n]-fh, format='(F8.3)'),1)
            printf, unit, string
         endfor
      endelse

      ; close text output file and free unit

      close, unit
      free_lun, unit

   endif

   if keyword_set(multi) then pause  ; ensures a pause between multiplots on multiple days

endfor

;!Quiet = 0

!p.multi = pmulti  ; restore original settings
;!x.margin = xmargin
;!y.margin = ymargin

end
