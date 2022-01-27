;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
forward_function t_prob

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function ratioexpectation,fracsigma_numer,fracsigma_denom,rms=rms

; Compute the expectation value of the ratio of two independent normally
; distributed variates, for the case where the denominator is positive-definite,
; the two means are equal, the standard deviation of the numerator is a fraction
; fracsigma_numer of the mean, and the standard deviation of the denominator
; is a fraction fracsigma_denom of the mean.
;
; The expectation value is greater than one, even though the means are equal,
; because positive fluctuations in the denominator produce negative fluctuations
; in the ratio but negative fluctuations in the denominator produce even larger
; positive fluctuations in the ratio.
;
; This routine is meant to be used for CW dehopping:

; Assuming that the ratio of noise-free spectrum to background spectrum is
; unity prior to subtraction leaves a small positive bias in the
; subtracted, normalized spectra, particularly noticeable (several tenths
; of a standard deviation) at fine frequency resolution.
;
; The expected rms deviation of this ratio about the mean ratio can be returned
; via the rms argument
;
; CM  2002 Oct 23

ratio_min = 0.d0       ; power can't be negative
ratio_max = 3.d0       ; high enough so that sum of pdf is ~ 1
ratio_spacing = 1.d-4

fracvar_numer = (1.d0*fracsigma_numer)^2
fracvar_denom = (1.d0*fracsigma_denom)^2
n_points = round((ratio_max - ratio_min)/ratio_spacing)
ratio = (ratio_max - ratio_min)*(dindgen(n_points) + 0.5)/n_points + ratio_min

pdf_factor1 = 1/sqrt(2*!DPI)
pdf_factor2 = (fracvar_denom*ratio + fracvar_numer)             $
              / (fracvar_denom*ratio^2 + fracvar_numer)^1.5
exponent_factor3 = -0.5*( (ratio - 1)^2 / (fracvar_denom*ratio^2 + fracvar_numer) ) $
                   > (-5.0d2)
pdf_factor3 = exp(exponent_factor3)
pdf = pdf_factor1*pdf_factor2*pdf_factor3*ratio_spacing

pdfsum = total(pdf)
if pdfsum lt 0.99 then begin
  print,' '
  print,'WARNING in ratioexpectation: Large fractional fluctuations (sum of pdf = ', $
        pdfsum,format='(a,f5.3)'
  print,'                             -- must increase ratio_max in routine'
  print,' '
endif

meanratio = total(pdf*ratio)
if arg_present(rms) then rms = float(sqrt(total(pdf*(ratio-meanratio)^2)))
return, float(meanratio)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro setup,parameter,value,reset=reset,stackreset=stackreset,help=help,list=list

; List or set the dehopping parameters,
; or reset the dehopping parameters and the tags (and, if desired, the stacks)

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common setup,srate,hopbw,hopf1,txoff
common target,tname,diameter,period,jd0,phase0

if keyword_set(help) or n_params() eq 1 or n_params() gt 2 then begin
  print,'setup,param,value[,/reset][,/list][,/help]'
  print,'where param can be: file,fstem,fsuf,npts,hopbw,hopf1,nhops,hop1,txoff,', $
             'dwell,srate,date,',format='(2a)'
  print,'                    zfile,zstem,target,nblock,lambda,tzcorr,diam,period'
  return
endif

; List the dehopping parameters if requested

if keyword_set(list) or (n_params() eq 0 and not keyword_set(reset)) then begin
  showparams
  return
endif

; Reset the dehopping parameters and the tags if requested 

if keyword_set(reset) then begin
  if keyword_set(stackreset) then $
    resetparams,/stackreset $
  else $
    resetparams
  return
endif

; Set the desired dehopping parameter

case parameter of
  'file'  : begin
              if size(value, /type) ne 7 then begin
                print,' '
                print,'ERROR: Bad value of file'
                print,'Make sure that file is a quoted string!'
                print,' '
              endif else begin
                infile = value
                filestem = ''
                filesuffix = ''
              endelse
            endcase
  'fstem' : begin
              if size(value, /type) ne 7 then begin
                print,' '
                print,'ERROR: Bad value of fstem'
                print,'Make sure that fstem is a quoted string!'
                print,' '
              endif else begin
                filestem = value
                infile = ''
              endelse
            endcase
  'zfile' : begin
              if size(value, /type) ne 7 then begin
                print,' '
                print,'ERROR: Bad value of zfile'
                print,'Make sure that zfile is a quoted string!'
                print,' '
              endif else begin
                zfile = value
                zstem = ''
                filesuffix = ''
                if notnull(zfile) then readazel
                if isnull(zfile) and n_zdata gt 0 then begin
                  n_zdata = 0L
                  radec = 0
                  print,'Azel info zeroed out'
                endif
              endelse
            endcase
  'zstem' : begin
              if size(value, /type) ne 7 then begin
                print,' '
                print,'ERROR: Bad value of zstem'
                print,'Make sure that zstem is a quoted string!'
                print,' '
              endif else begin
                zstem = value
                zfile = ''
                if notnull(zstem) then begin
                  print,'Azel file will be read when fsuf is (re)set'
                endif else if n_zdata gt 0 then begin
                  n_zdata = 0L
                  radec = 0
                  print,'Azel info zeroed out'
                endif
              endelse
            endcase
  'fsuf' :  begin
              if size(value, /type) ne 7 and  $
                    not (size(value, /type) ge 1 and size(value, /type) le 3) then begin
                print,' '
                print,'ERROR: Bad value of fsuf'
                print,'Make sure that fsuf is a quoted string or integer!'
                print,' '
              endif else begin
                if size(value, /type) eq 7 then begin
                  filesuffix = (strmid(value,0,1) eq '.') ? value : '.' + value
                endif else begin
                  filesuffix = '.' + string(value,format='(i0)')
                endelse
                infile = ''
                if notnull(zstem) then begin
                  zfile = zstem + filesuffix
                  readazel
                  if isnull(zfile) then begin
                    filesuffix = ''
                    if n_zdata gt 0 then begin
                      n_zdata = 0L
                      radec = 0
                      print,'Azel info zeroed out'
                    endif
                  endif
                endif
                zfile = ''
              endelse
            endcase
  'npts'  : begin
              npts = round(1L*value)
              if 2*(npts/2L) ne npts then begin
                print,' '
                print,'WARNING: odd number of fft points  (npts =',npts,')', $
                      format='(a,i0,a)'
                print,'         -- the fft and shift routines which got you here'
                print,'            may or may not properly handle such spectra'
                print,' '
              endif
              df = (npts ne 0) ? float(srate)/npts : 0.0
            endcase
  'srate' : begin
              srate = 1d0*value
              df = (npts ne 0) ? float(srate)/npts : 0.0
            endcase
  'nblock': begin
              if value ge 1 then begin
                nblock = round(1L*value)
              endif else begin
                print,'ERROR: nblock = ',value,' is too low, ', $
                      'so we instead set nblock = 1',format='(a,i0,a)'
                nblock = 1L
              endelse
            endcase
  'hopbw' : hopbw = 1.0*value
  'hopf1' : hopf1 = 1.0*value
  'nhops' : nhops = round(1L*value)
  'hop1'  : hop1 = round(1L*value)
  'txoff' : txoff = 1.0*value
  'dwell' : dwell = 1d0*value
  'lambda': lambda = 1.0*value
  'date'  : date = julday(value[0],value[1],value[2])
  'tzcorr': begin
              tzcorr = round(1L*value)
              if tzcorr eq 0 then begin
                outputTimeZone = 'AST'
              endif else if tzcorr eq 4 then begin
                outputTimeZone = 'UTC'
              endif else begin
                outputTimeZone = ''
              endelse
            endcase
  'target': begin
              if size(value, /type) ne 7 then begin
                print,' '
                print,'ERROR: Bad value of target'
                print,'Make sure that target is a quoted string!'
                print,' '
              endif else begin
                tname = value
                case strlowcase(tname) of
                  'europa'  : begin
                                diameter = 3138.0
                                period = 24*3.55
                              endcase
                  'ganymede': begin
                                diameter = 5262.0
                                period = 24*7.16
                              endcase
                  'callisto': begin
                                diameter = 4800.0
                                period = 24*16.69
                              endcase
                  'titan'   : begin
                                diameter = 5150.0
                                period = 24*15.95
                              endcase
                  else      :
                endcase
              endelse
            endcase
  'diam'  : diameter = 1.0*value
  'period': period = 1.0*value
  else    : begin
              print,'ERROR: Parameter ',parameter,' is unknown and will be ignored'
              print,'Legal parameters: file,fstem,fsuf,npts,hopbw,hopf1,', $
                         'nhops,hop1,txoff,dwell,srate,date,',format='(2a)'
              print,'                  zfile,zstem,target,nblock,lambda,tzcorr,diam,period'
            endcase
endcase

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro resetparams,help=help,stackreset=stackreset

; Reset the dehopping parameters, the tags, and (if requested) the stacks

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common setup,srate,hopbw,hopf1,txoff
common target,tname,diameter,period,jd0,phase0
common tagIO,tagtypes,tagformats,nameformat,iowidth,iodecplaces,iomaxline
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if keyword_set(help) or n_params() gt 0 then begin
  print,'resetparams[,/help,/stackreset]'
  return
endif

; Reset dehopping parameters

infile = ''             ; spectra file
filestem = ''           ; stem of spectra files
zfile = ''              ; azel file
zstem = ''              ; stem of azel files
filesuffix = ''         ; suffix of spectra and azel files (after '.p1' or '.p2' or 'azel')
tname = ''              ; target name
npts = 6250L            ; points per spectrum
hopbw = 10000.0         ; Hz
hopf1 = -15000.0        ; Hz
txoff = 0.0             ; Hz
nhops = 4L
hop1 = 0L               ; first hop number
dwell = 10d0            ; sec
srate = 62500d0         ; Hz
nblock = 6L             ; hop sequences to sum before processing
date = 0L               ; [mm,dd,yyyy]
lambda = 0.1259632      ; m
period = 1.0            ; hr
diameter = 1.0          ; km
tzcorr = 0L             ; hr
outputTimeZone = 'AST'

df = (npts ne 0) ? float(srate)/npts : 0.0  ; Hz
n_zdata = 0L                         ; # of azel data points read in
radec = 0                            ; no right ascension / declination data

; Unload (reset) the loaded pair and the loaded single-channel spectrum

unload
unload1

; Zero out the stacks if requested

if keyword_set(stackreset) then begin
  resetstack
  resetstack1
endif

; Reset some plotting parameters

!x.style = 1
!y.style = 3

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro showparams,help=help

; Display the dehopping parameters

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common setup,srate,hopbw,hopf1,txoff
common target,tname,diameter,period,jd0,phase0

if keyword_set(help) or n_params() gt 0 then begin
  print,'showparams[,/help]'
  return
endif

; List the parameters

print,'datafile:            ',infile
print,'datafile stem:       ',filestem
print,'za file:             ',zfile
print,'za file stem:        ',zstem
print,'data/za file suffix: ',filesuffix
print,'target:              ',tname
if date ne 0 then begin
  caldat,date,mon,day,year
  print,'date:                ',mon,'/',day,'/',year, $
        format = '(a,i2.2,a,i2.2,a,i4)'
endif else begin
  print,'date:                 0/ 0/   0'
endelse
print,'output time zone:    ',outputTimeZone,format='(a,a8)'
print,'spec points:         ',npts,format='(a,i8)'
print,'nblock:              ',nblock,format='(a,i8)
print,'nhops:               ',nhops,format='(a,i8)'
print,'1st hop (start@0):   ',hop1,format='(a,i8)'
print,'hop bandwidth:       ',hopbw,' Hz',format='(a,f8.1,a)'
print,'1st hop frequency:   ',hopf1,' Hz',format='(a,f8.1,a)'
print,'txoffset:            ',txoff,' Hz',format='(a,f8.1,a)'
print,'dwell time:          ',float(dwell),' s',format='(a,f8.1,a)'
print,'sample rate:         ',float(srate),' Hz',format='(a,f8.1,a)'
print,'freq. resolution:    ',df,' Hz',format='(a,f8.3,a)'
print,'wavelength:         ',lambda,' m',format='(a,f9.7,a)'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro readazel,help=help

; Open a file which was output from stripAzEl; read header information
; about the target and the observing setup; and read values of the
; azimuth, elevation, gain, system temperatures, etc., as a function of
; time during the run.

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common setup,srate,hopbw,hopf1,txoff
common target,tname,diameter,period,jd0,phase0

if keyword_set(help) or n_params() ne 0 then begin
  print,'readazel[,/help]'
  return
endif

; Open zfile

err = openInfile(lun,zfile,'azel',/get_lun)
if err ne 0 then return

; Read the azel header
; -- assume that each line starts with a 7-character label and then a space
; -- take extra care with the "TARGET" line in case the target name includes spaces

string = '' ; readf needs this defined as a string in advance, or else it assumes float
readf,lun,string,year,mon,day,format='(a8,3i0)'
date = julday(mon,day,year)
readf,lun,string,format='(a)'
pieces = stregex(strmid(string,8),'^(.+) ([^ ]+ [^ ]+ [^ ]+ [^ ]+)$', /subexpr, /extract)
tname_Z = pieces[1]

; jd0 (Julian date when rotation phase = phase0) has to be defined in advance as
; double-precision or else reads will assume float, with terrible roundoff consequences.

jd0 = 0.0D
reads,pieces[2],diameter_Z,period_Z,jd0,phase0,format='(4f0)'
readf,lun,string,nhops_Z,hop1_Z,hopbw_Z,hopf1_Z,dwell_Z,txoff_Z,format='(a8,2i0,4f0)'
readf,lun,string,srate_Z,format='(a8,f0)'
readf,lun,string,lambda_Z,format='(a8,f0)'

; Check whether or not delay and Doppler corrections are included

if eof(lun) then begin
  deldopcorr = 0
endif else begin
  readf,lun,string,format='(a8)'
  deldopcorr = (string eq 'EPHEM   ') ? 1 : 0
endelse
nheadlines = (deldopcorr) ? 6L : 5L

close,lun
free_lun,lun

; Let the user know about any dehopping parameters which are being changed
; as a result of new values being read in from the header of zfile.

if tname ne tname_Z then begin
  if notnull(tname) then print,'Switching from tname = ',tname,' to ',tname_Z
  tname = tname_Z
endif
if diameter ne diameter_Z then begin
  if diameter ne 1.0 then $
    print,'Switching from diameter = ',diameter,' km to ',diameter_Z,' km'
  diameter = 1.0*diameter_Z
endif
if period ne period_Z then begin
  if period ne 1.0 then $
    print,'Switching from period = ',period,' hr to ',period_Z,' hr'
  period = 1.0*period_Z
endif
if nhops ne nhops_Z then begin
  print,'Switching from nhops = ',nhops,' to ',nhops_Z
  nhops = 1L*nhops_Z
endif
if hop1 ne hop1_Z then begin
  print,'Switching from hop1 = ',hop1,' to ',hop1_Z
  hop1 = 1L*hop1_Z
endif
if hopbw ne hopbw_Z then begin
  print,'Switching from hopbw = ',hopbw,' Hz to ',hopbw_Z,' Hz'
  hopbw = hopbw_Z
endif
if hopf1 ne hopf1_Z then begin
  print,'Switching from hopf1 = ',hopf1,' Hz to ',hopf1_Z,' Hz'
  hopf1 = hopf1_Z
endif
if float(dwell) ne dwell_Z then begin
  print,'Switching from dwell = ',float(dwell),' s to ',dwell_Z,' s'
  dwell = double(dwell_Z)
endif
if txoff ne txoff_Z then begin
  print,'Switching from txoff = ',txoff,' Hz to ',txoff_Z,' Hz'
  txoff = txoff_Z
endif
if float(srate) ne srate_Z then begin
  print,'Switching from srate = ',float(srate),' Hz to ',srate_Z,' Hz'
  srate = double(srate_Z)
endif
if lambda ne lambda_Z then begin
  print,'Switching from lambda = ',lambda,' m to ',lambda_Z,' m',format='(a,f9.7,a,f9.7,a)'
  lambda = 1.0*lambda_Z
endif

; Read the azel data into a dummy structure:
; hop color, azimuth, elevation, gain, Tsys, etc., for each dwell of the run

template = $
    {version:1.0,datastart:nheadlines,delimiter:32b, $
     missingvalue:!VALUES.F_NAN,commentsymbol:'', $
     fieldcount:14L,fieldtypes:[3,4,3,4,4,4,4,4,4,4,4,4,5,5], $
     fieldnames:["rec","sec","color","az","za","t1","t2","grgt","pwr","rtt","ra","dec","delcorr","dopcorr"], $
     fieldlocations:[5L,7L,15L,20L,30L,37L,45L,52L,60L,69L,78L,87L,99L,111L], $
     fieldgroups:[0,1,2,3,4,5,6,7,8,9,10,11,12,13]}
zdummy = read_ascii(zfile, template=template, count=n_zdata)

if n_zdata eq 0 then begin
  print,'ERROR: No azel data found in file ',zfile
  return
endif

; After multiplication by 1e12,
; grgt is (dimensionless rx gain)*(dimensionless tx gain)

zdummy.grgt = zdummy.grgt * 1.0e12

; Construct two azimuth vectors, using the ranges [0,360) and [180,540).

az1 = zdummy.az - 360*floor(zdummy.az/360.0)
az2 = az1 - 360*floor( (az1 - 180)/360.0 )

; Check whether the azel file contains valid right ascension and declination data
; (it won't if it was created prior to April 2002, or if stripAzEl was run without
; the ephemeris option).  If it does, construct two RA vectors, using the ranges
; [0,24) and [12,36) hours.

if not finite(zdummy.ra[0]) then begin
  radec = 0
endif else begin
  radec = zdummy.ra[0] ge 0 and zdummy.ra[0] lt 24
endelse
if radec then begin
  ra1 = zdummy.ra - 24*floor(zdummy.ra/24.0)
  ra2 = ra1 - 24*floor( (ra1 - 12.0)/24.0 )
endif else begin
  ra1 = zdummy.ra   ;  Can't do math on NaN values
  ra2 = zdummy.ra
endelse

; Create the zdata structure which includes all of these vectors

zdata = {rec:zdummy.rec, sec:zdummy.sec, color:zdummy.color, $
         az1:az1, az2:az2, za:zdummy.za, t1:zdummy.t1, t2:zdummy.t2, $
         grgt:zdummy.grgt, pwr:zdummy.pwr, rtt:zdummy.rtt, $
         ra1:ra1, ra2:ra2, dec:zdummy.dec, delcorr:zdummy.delcorr, dopcorr:zdummy.dopcorr}

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro dehop,flipchan=flipchan,npol=npol,help=help,_extra=_ext

; Dehop an OC/SC spectral pair
;
; Modified 2008 Jul 6 by CM to add the degreemax and plotfit keywords

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'dehop[,/noweight][,/merge OR ,/blocks][,maxf=maxf]    $'
  print,'     [,/nobgfit][,/bias_allow][,degreemax=degreemax]  $'
  print,'     [,/roundephcorr][,flipchan=flipchan][,npol=npol][,/plotfit][,/help]'
  print,' '
  print,'      /merge produces one output spectrum for the entire file, ', $
        'rather than one per scan',format='(2a)'
  print,'      /blocks produces one output spectrum per block rather than one per scan'
  print,'      maxf only processes frequency range [-maxf,+maxf] for each hop'
  print,'      /nobgfit normalizes to background spectra rather than ', $
        'to polynomial fits to background spectra',format='(2a)'
  print,'      /bias_allow is relevant only if /nobgfit is set;'
  print,'            it prevents subtraction of the bias produced when dividing ', $
        'by noisy background spectra',format='(2a)'
  print,'      degreemax is the maximum polynomial degree allowed when fitting ', $
        'background spectra',format='(2a)'
  print,'            (default = 4)'
  print,'      /roundephcorr rounds Doppler ephemeris corrections so that the spectra'
  print,'            are shifted by an integer number of bins; otherwise no rounding'
  print,'            is carried out and interpolation (cubic convolution) is performed'
  print,'            rather than simple shifting.'
  print,'
  print,'            (Note that interpolation can artificially reduce the noise variance.)'
  print,'
  print,'      flipchan is used to correct a cabling error (switching I and Q):'
  print,'            flipping is done about the middle of the spectrum, and neither'
  print,'            xjcen, posfr, nor the frequency vector are changed.  Flipping is'
  print,'            done before any shifting is carried out due to /roundephcorr.'
  print,'            Permitted values are flipchan = 1, 2, or 3 (flip both channels).'
  print,'      npol sets the number of channels it reads. npol=2 (default) is normal'
  print,'            OC/SC processing. npol=4 adds the RE and Im Stokes data. Note:'
  print,'            flipchan must be 0 or 3 if npol=4: you can''t flop the crosses'
  print,'      /plotfit plots the polynomial fit to the background spectrum ', $
        'for each hop of each block',format='(2a)'
  print,' '
  return
endif

; Default to standard processing, mainly because file checking is done at a lower level.
if not keyword_set(npol) then npol=2
if npol lt 2 or npol gt 4 then begin
  print, 'ERROR: Dehop: npol mist be 2, 3, or 4'
  return
endif

if (notnull(zstem)) and (isnull(filesuffix)) then begin
  print,"ERROR in dehop: You forgot to set fsuf, so zfile isn't specified"
  return
endif

; Deal with the flipchan flag

if not keyword_set(flipchan) then begin
  iqerror_OC = 0
  iqerror_SC = 0
endif else begin
  if npol gt 2 and flipchan ne 3 then begin
    print, "ERROR in dehop: Don't know how to do Stokes with differently-flipped spectra"
    return
  endif
  iqerror_OC = (flipchan eq 1 or flipchan eq 3) ? 1 : 0
  iqerror_SC = (flipchan eq 2 or flipchan eq 3) ? 1 : 0
endelse

; Find out how many spectra are in the single-channel stack at the outset

nstack1_start = nstack1

; Store the loaded single-channel spectrum so we can use that "space" and
; reload the spectrum later on

storeloaded1 = ptr_new(*loaded1)

; Dehop each polarization in turn, each time setting the datafile name
; to the appropriate value ('.p1' or '.p2' suffix for OC, SC, respectively).
; If necessary set the azel filename as well.
;
; Start with the OC data; the output spectra will be added to the
; single-channel stack

infile = filestem + '.p1' + filesuffix
if notnull(filesuffix) then zfile = zstem + filesuffix
dehop1,1,iqerror=iqerror_OC,/silent,_extra=_ext
n_OC = nstack1 - nstack1_start

; Dehop the SC data

infile = filestem + '.p2' + filesuffix
dehop1,2,iqerror=iqerror_SC,/silent,_extra=_ext
n_SC = nstack1 - (nstack1_start + n_OC)

if npol gt 2 then begin

; Dehop the RE data

infile = filestem + '.p3' + filesuffix
dehop1,2,iqerror=iqerror_SC,/silent,_extra=_ext
n_RE = nstack1 - (nstack1_start + n_SC + n_OC)

endif

if npol gt 3 then begin

; Dehop the IM data

infile = filestem + '.p4' + filesuffix
dehop1,2,iqerror=iqerror_SC,/silent,_extra=_ext
n_IM = nstack1 - (nstack1_start + n_RE + n_SC + n_OC)

endif

; Decide what to do with the two channels, then do it

if npol eq 2 then processRaw_combineChans,n_OC,n_SC,nstack1_start
if npol eq 3 then processRaw_combineChans,n_OC,n_SC,nstack1_start,n_RE
if npol eq 4 then processRaw_combineChans,n_OC,n_SC,nstack1_start,n_RE,n_IM

; If everything worked correctly and all output channels were spliced into
; OC/SC pairs, reload the single-channel spectrum that was there at the start

if npol eq 2 then if n_OC eq n_SC then *loaded1 = *storeloaded1 else $
if npol eq 3 then if n_OC eq n_SC and n_OC eq n_RE then *loaded1 = *storeloaded1 else $
if npol eq 4 then if n_OC eq n_SC and n_OC eq n_RE and n_OC eq n_IM then *loaded1 = *storeloaded1

; Reset infile and zfile to their previous values and clean up pointers

infile = ''
if notnull(filesuffix) then zfile = ''
ptr_free,storeloaded1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro dehop1,chan,noweight=noweight,merge=merge,blocks=blocks,maxf=maxf,nobgfit=nobgfit, $
           bias_allow=bias_allow,degreemax=degreemax,roundephcorr=roundephcorr,        $
           iqerror=iqerror,plotfit=plotfit,silent=silent,help=help

; Dehop a single-channel spectrum and load the result
;
; Modified by CM 2003 Dec 13: Check whether nout > 0 before pushing the spectrum onto the
;          single-channel stack, so that we don't get the last complete block pushed twice
;          when /blocks is set and only unused data follow this block
;
; Modified 2008 Jul 6 by CM to add the degreemax and plotfit keywords

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common target,tname,diameter,period,jd0,phase0
common setup,srate,hopbw,hopf1,txoff
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common loadedBlock,loadedi,loaded1,loaded
common fileInit,doperr,nper,gconst,tsys,sdev_over_rmsc,yin,y,lun,nrec,endOfScan,endOfFile,nscan
common scanInit,r,nout,iblock,nrec_scan,rcsta_min,rcend_max,az1_min,az1_max,az2_min, $
                az2_max,rctime_block,za_block,az1_block,az2_block,rtt_block,pwr_block, $
                tsys_block,gain_block,rctime_mean,za_mean,az1_mean,az2_mean,rtt_mean, $
                pwr_mean,tsys_mean,gain_mean,nffts,tau,sdev_over_rmsc_block,hop, $
                ra1_min,ra1_max,ra2_min,ra2_max,ra1_block,ra2_block,dec_block, $
                ra1_mean,ra2_mean,dec_mean,dec_min,dec_max,rtt_min,rtt_max,blocknumWithinScan

if not keyword_set(merge) then merge = 0
if not keyword_set(blocks) then blocks = 0
if not keyword_set(roundephcorr) then roundephcorr = 0
if not keyword_set(iqerror) then iqerror = 0
if not keyword_set(plotfit) then plotfit = 0
if n_elements(degreemax) eq 0 then degreemax = -1

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'dehop1,chan[,/noweight][,/merge OR ,/blocks][,maxf=maxf][,/nobgfit]  $'
  print,'           [,/bias_allow][,degreemax=degreemax][,/roundephcorr]      $'
  print,'           [,/iqerror][,/plotfit][,/silent][,/help]'
  print,' '
  print,'       /merge produces one output spectrum for the entire file, ', $
        'rather than one per scan',format='(2a)'
  print,'       /blocks produces one output spectrum per block rather than one per scan'
  print,'       maxf only processes frequency range [-maxf,+maxf] for each hop'
  print,'       /nobgfit normalizes to background spectra rather than ', $
        'to polynomial fits to background spectra',format='(2a)'
  print,'       /bias_allow is relevant only if /nobgfit is set;'
  print,'             it prevents subtraction of the bias produced when dividing ', $
        'by noisy background spectra',format='(2a)'
  print,'       degreemax is the maximum polynomial degree allowed when fitting ', $
        'background spectra',format='(2a)'
  print,'             (default = 4)'
  print,'       /roundephcorr rounds Doppler ephemeris corrections so that the spectrum'
  print,'             is shifted by an integer number of bins; otherwise no rounding'
  print,'             is carried out and interpolation (cubic convolution) is performed'
  print,'             rather than simple shifting'
  print,'
  print,'             (Note that interpolation can artificially reduce the noise variance.)'
  print,'
  print,'       /iqerror is used to correct a cabling error (switching I and Q):'
  print,'             the spectrum is flipped about its middle, and neither xjcen,'
  print,'             posfr, nor the frequency vector are changed.  Flipping is'
  print,'             done before any shifting is carried out due to /roundephcorr.'
  print,'       /plotfit plots the polynomial fit to the background spectrum ', $
        'for each hop of each block',format='(2a)'
  print,'       /silent suppresses output when the processed spectrum is added to the ', $
        'single-channel stack',format='(2a)'
  print,' '
  return
endif else if chan ne 1 and chan ne 2 then begin
  print,' '
  print,'dehop1,chan[,/noweight][,/merge OR ,/blocks][,maxf=maxf][,/nobgfit]  $'
  print,'           [,/bias_allow][,degreemax=degreemax][,/roundephcorr]      $'
  print,'           [,/iqerror][,/plotfit][,/silent][,/help]'
  print,'ERROR in dehop1: Must have chan = 1 (OC) or 2 (SC)'
  print,' '
  return
endif else if merge and blocks then begin
  print,' '
  print,'dehop1,chan[,/noweight][,/merge OR ,/blocks][,maxf=maxf][,/nobgfit]  $'
  print,'           [,/bias_allow][,degreemax=degreemax][,/roundephcorr]      $'
  print,'           [,/iqerror][,/plotfit][,/silent][,/help]'
  print,"ERROR in dehop1: Can't set both /merge and /blocks"
  print,' '
  return
endif else if (notnull(zstem)) and (isnull(filesuffix)) then begin
  print,"ERROR in dehop1: You forgot to set fsuf, so zfile isn't specified"
  return
endif

if n_elements(maxf) gt 0 then begin
  maxf = 1.0*maxf
  if maxf le 0 or maxf gt hopbw/2 then begin
    print,'ERROR in dehop1: Must have 0 < maxf <= hopbw/2'
    print,' '
    return
  endif
endif else begin
  maxf = hopbw/2
endelse

; Print a warning about hop colors for data with no azel info

if n_zdata eq 0 then begin
  print,' '
  print,'WARNING in dehop1: infile should contain just ONE scan'
  print,"                   -- without azel info, there's no way to synch the hop cycle ", $
        "beyond scan 1",format='(2a)'
  print,' '
endif

; Initialize various quantities before starting work

processRaw_fileInit1,chan,noweight=noweight,maxf=maxf

; Start reading spectra

repeat begin

  ; Initialize the next output spectrum

  processRaw_scanInit1

  ; Process the blocks of raw spectra which go into that output spectrum

  repeat begin
    processRaw_doBlock1,chan,merge=merge,nobgfit=nobgfit,bias_allow=bias_allow,             $
                        degreemax=degreemax,roundephcorr=roundephcorr,iqerror=iqerror, $
                        plotfit=plotfit
  endrep until (endOfFile or blocks or (endOfScan and not merge))

  ; If any "good" blocks were processed, create the finished spectrum

  if nout gt 0 then begin

    ; Sum the blocks, set the tags, and load the output spectrum

    processRaw_sumBlocks1,chan,merge=merge

    ; Push the new spectrum onto the single-channel stack
    ; (especially important if there are more spectra to come)

    push1,silent=silent

  endif

endrep until (endOfFile)

; All raw spectra have been read in

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro processRaw_fileInit1,chan,hoptouse=hoptouse,noweight=noweight,maxf=maxf

; When processing hopped or unhopped data,
; initialize various quantities prior to reading any raw data
; (i.e., before starting work on the first processed output spectrum)

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common target,tname,diameter,period,jd0,phase0
common setup,srate,hopbw,hopf1,txoff
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common fileInit,doperr,nper,gconst,tsys,sdev_over_rmsc,yin,y,lun,nrec,endOfScan,endOfFile,nscan

if not keyword_set(noweight) then noweight = 0

; Set doperr (ephemeris error for some Titan runs)
; -- note that this parameter isn't yet used in the reduction

if strlowcase(tname) eq 'titan' then begin
  case date of
    julday(10,23,1999): doperr = -896.44
    julday(10,25,1999): doperr = 2144.32
    julday(10,27,1999): doperr = 3998.17
    julday(10,29,1999): doperr = 647.21
    julday(10,31,1999): doperr = -3456.35
    else              : doperr = 0.0
  endcase
endif

; Compute lbin and rbin, vectors containing the leftmost and rightmost bin numbers
; to be saved for each hop; if the hop bandwidth isn't an integer number of bins
; (e.g., due to being Doppler shifted) then roundoff must be dealt with so that
; all hop colors have the same number of spectral points (nper).
;
; If the data are unhopped, or if we're saving just one hop color, keep
; a single spectrum covering the full unaliased bandwidth

if n_elements(hoptouse) eq 0 then begin
  if n_elements(maxf) gt 0 then begin
    maxf = 1.0*maxf
  endif else begin
    maxf = hopbw/2
  endelse
  cbin = round( npts/2.0 + npts*(txoff + hopf1 + findgen(nhops)*hopbw)/srate )
  if 2*(npts/2L) ne npts then cbin = cbin - 1L
  nper = round(2*npts*maxf/srate)
  if 2*(nper/2L) ne nper then nper = nper - 1L
  lbin = cbin - nper/2L
  rbin = lbin + nper - 1L
endif else begin
  cbin = npts/2L
  nper = 2*cbin
  lbin = [0L]
  rbin = [nper - 1L]
endelse

; To get sdev, the cross section equivalent of the rms noise power,
; equate the received signal power (given by the radar equation)
; to the rms noise power fluctuation (k*Tsys*df*rmsc) and solve for the
; cross section.  Here df is the frequency resolution and rmsc is the
; calculated fractional rms noise power.
;
; For frequency-hopped data with nhops colors per hop cycle, only 
; complete cycles used, and dehopping done using raw background spectra
; for baselining and normalization (i.e., /nobgfit set), we have
;      rmsc = sqrt{ nhops / [(nhops-1)*df*tau] },
; with tau the integration time for a full block being summed.
; See Ostro et al. 1992, JGR, vol 97, no. E11, pp. 18227-18244.
;
; For unhopped data, rmsc = 1/sqrt(df*tau)
;
; For hopped data which are dehopped using low-order polynomial fits to
; background spectra for baselining and normalization (i.e., /nobgfit
; *not* set), rmsc should be only slightly larger than the value for
; unhopped spectra.

; gconst = (4 pi)^3 * c^4 * k / (16 * 1e3 * 1e6 * lambda^2)
; -- this is a bunch of constants and unit conversion factors
;    taken from the radar equation and from the expression for 
;    rms noise power

gconst = 871730.0 ; Arecibo (lambda = 0.1259632 m)

; calculate sdev/rmsc for each raw spectrum
; (rmsc involves tau, and we don't yet know how many dwells we'll use
; in the sum so we don't yet know the integration time)

if n_zdata gt 0 then begin
  if chan lt 3 then tsys = (chan eq 1) ? zdata.t1 : zdata.t2 $
  else tsys = sqrt(zdata.t1 * zdata.t2)
endif

if noweight or (n_zdata eq 0) then begin

  ; uniform weighting if no other info available
  ; (do this for 10000 dwells -- surely enough for any scan)

  sdev_over_rmsc = fltarr(10000) + 1.0

endif else begin

  ; In the next expression tsys is in K, rtt in sec, df in Hz, pwr in kW and
  ; grgt = (rx gain)*(tx gain) is dimensionless
  ;
  ; If any of these are set to zero or negative values for a given dwell,
  ; set the corresponding sdev_over_rmsc value to be a very large negative
  ; number; that way, if somehow this dwell accidentally gets used later on,
  ; the block to which it contributes will have a large negative sdev, and
  ; you'll see this if you're paying attention to the output display

  sdev_over_rmsc = fltarr(n_zdata)
  for n = 0L, n_zdata-1 do begin
    if tsys[n] gt 0 and zdata.rtt[n] gt 0  $
                    and zdata.pwr[n] gt 0 and zdata.grgt[n] gt 0 then begin
      sdev_over_rmsc[n] = gconst * tsys[n] * zdata.rtt[n]^4 * df  $
                          / (zdata.pwr[n] * zdata.grgt[n])
    endif else begin
      sdev_over_rmsc[n] = -9.99e20
    endelse
  endfor
endelse

; Define some spectral arrays for each hop cycle:
;     yin = raw spectrum
;       y = raw spectra sorted by hop color

yin = fltarr(npts)
if n_elements(hoptouse) eq 0 then begin
  y = fltarr(nhops,npts)                ; hopped data to be dehopped
endif else if hoptouse ne -1 then begin
  y = fltarr(nhops,npts)                ; hopped data, saving just one hop color
endif else begin
  y = fltarr(1,npts)                    ; unhopped data
endelse

; Open input file and initialize the number of raw spectra read in

err = openInfile(lun,infile,/get_lun)
if err ne 0 then retall

nrec = 0L   ; total individual raw spectra encountered (but not necessarily used)
nscan = 0L  ; initialize the number (counting from 1) of the scan currently being read in

endOfScan = 1L   ; time to start a new scan
endOfFile = 0L   ; we're actually at the start of the file

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro processRaw_scanInit1

; When processing hopped or unhopped data,
; initialize various quantities prior to processing each scan 
; (i.e., when starting to work on each processed output spectrum)

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common fileInit,doperr,nper,gconst,tsys,sdev_over_rmsc,yin,y,lun,nrec,endOfScan,endOfFile,nscan
common scanInit,r,nout,iblock,nrec_scan,rcsta_min,rcend_max,az1_min,az1_max,az2_min, $
                az2_max,rctime_block,za_block,az1_block,az2_block,rtt_block,pwr_block, $
                tsys_block,gain_block,rctime_mean,za_mean,az1_mean,az2_mean,rtt_mean, $
                pwr_mean,tsys_mean,gain_mean,nffts,tau,sdev_over_rmsc_block,hop, $
                ra1_min,ra1_max,ra2_min,ra2_max,ra1_block,ra2_block,dec_block, $
                ra1_mean,ra2_mean,dec_mean,dec_min,dec_max,rtt_min,rtt_max,blocknumWithinScan

; Print the header line of the block-parameter display

print,'  chan    block#    nseq      aveback         sdev         rmsc          tau'

; Initialize an array holding the sum of all raw spectra in a block for each hop color

r = fltarr(nhops,npts)

; Initialize counters

nout = 0L        ; # of input blocks contributing to the next output spectrum
iblock = 0L      ; # of raw spectra summed into block
nrec_scan = 0L   ; # of raw spectra encountered in this scan (but not necessarily used)
if endOfScan then begin
  blocknumWithinScan = 0L ; block number within this scan
endif

; Initialize parameters used to compute tag values for the processed spectrum

if n_zdata gt 0 then begin
  rcsta_min = 999999.9D
  rcend_max = -999999.9D
  az1_min = 999.9
  az1_max = -999.9
  az2_min = 999.9
  az2_max = -999.9
  rtt_min = 999999.9
  rtt_max = -999999.9
  rctime_block = 0.0D
  za_block = 0.0
  az1_block = 0.0
  az2_block = 0.0
  rtt_block = 0.0
  pwr_block = 0.0
  tsys_block = 0.0
  gain_block = 0.0
  rctime_mean = 0.0D
  za_mean = 0.0
  az1_mean = 0.0
  az2_mean = 0.0
  rtt_mean = 0.0
  pwr_mean = 0.0
  tsys_mean = 0.0
  gain_mean = 0.0
  if radec then begin
    ra1_min = 999.9
    ra1_max = -999.9
    ra2_min = 999.9
    ra2_max = -999.9
    dec_min = 999.9
    dec_max = -999.9
    ra1_block = 0.0
    ra2_block = 0.0
    dec_block = 0.0
    ra1_mean = 0.0
    ra2_mean = 0.0
    dec_mean = 0.0
  endif
endif

nffts = 0L
tau = 0.0
sdev_over_rmsc_block = 0.0

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro processRaw_doBlock1,chan,merge=merge,hoptouse=hoptouse,nobgfit=nobgfit, $
                        bias_allow=bias_allow,degreemax=degreemax,          $
                        roundephcorr=roundephcorr,iqerror=iqerror,plotfit=plotfit

; Construct, sum, and (if data were hopped) background-subtract
; and normalize a block of raw spectra
;
; The two output parameters, endOfFile and endOfScan (in common block fileInit), let other
; procedures know whether or not there are more blocks to process for a given scan or file.
;
; Spectra are normalized using a polynomial fit to the background,
; unless nobgfit is set, in which case the background spectrum itself is used.
; Division by noisy background spectra results in dehopped spectra which are biased
; to positive values (see header comments in function ratioexpectation).
;
; Modified 2003 Jul 18 by CM to go on to the next highest polynomial degree rather than choking
;          if a baseline fit is impossible (e.g., if input measurement errors are negative).
;
; Modified 2008 Jul 6 by CM to add the chan argument and the degreemax and plotfit keywords

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common setup,srate,hopbw,hopf1,txoff
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common fileInit,doperr,nper,gconst,tsys,sdev_over_rmsc,yin,y,lun,nrec,endOfScan,endOfFile,nscan
common scanInit,r,nout,iblock,nrec_scan,rcsta_min,rcend_max,az1_min,az1_max,az2_min, $
                az2_max,rctime_block,za_block,az1_block,az2_block,rtt_block,pwr_block, $
                tsys_block,gain_block,rctime_mean,za_mean,az1_mean,az2_mean,rtt_mean, $
                pwr_mean,tsys_mean,gain_mean,nffts,tau,sdev_over_rmsc_block,hop, $
                ra1_min,ra1_max,ra2_min,ra2_max,ra1_block,ra2_block,dec_block, $
                ra1_mean,ra2_mean,dec_mean,dec_min,dec_max,rtt_min,rtt_max,blocknumWithinScan
common blockSpectra,s,sdev_block,rmsc_block,bias_block
chanstrings = ['OC','SC','RE','IM']

chanstring = chanstrings[chan-1]

if not keyword_set(nobgfit) then begin
  nobgfit = 0
  xfit = 2*dindgen(nper)/(nper - 1) - 1  ; "frequency" vector for background fitting
endif
if not keyword_set(bias_allow) then bias_allow = 0
if not keyword_set(iqerror) then iqerror = 0
if not keyword_set(plotfit) then plotfit = 0
if n_elements(degreemax) eq 0 then begin
  degreemax = 4
endif else if degreemax lt 0 then begin
  degreemax = 4
endif

; Determine whether or not the data are hopped, and also how many spectra you'd
; like to read at a time if you don't hit EOF or the end of a scan

if n_elements(hoptouse) gt 0 then begin
  hopped = hoptouse ne -1
  nin_desired = 1L
endif else begin
  hopped = 1
  nin_desired = nhops
endelse

; Loop through all raw spectra in a block

repeat begin

  if endOfScan then nscan = nscan + 1  ; about to start a new scan

  ; Read in a whole hop cycle
  ; (unless end of scan or end of azel file is reached first).
  ;
  ; If dealing with unhopped data, or if saving just one color of hopped data,
  ; read in one spectrum at a time rather than a cycle
  ;
  ; Do NOT allow any hop cycle to include data from two different scans.

  if n_elements(hoptouse) gt 0 then begin
    nread = 1L
  endif else if n_zdata eq 0 then begin
    nread = nhops
  endif else begin
    nreadmax = nhops < (n_zdata - nrec)
    zeroRec = where(zdata.rec[nrec:nrec+nreadmax-1] eq 0, count)
    if zdata.rec[nrec] gt 0 then begin
      nread = (count eq 0) ? nreadmax : zeroRec[0]
    endif else begin
      nread = (count eq 1) ? nreadmax : zeroRec[1]
    endelse
  endelse

  nin = 0L
  for k=0L,nread-1 do begin
    if not eof(lun) then begin 
      readu,lun,yin
      yin = swap_endian(temporary(yin), /swap_if_little_endian)
      if not hopped then begin
        hop = 0L
      endif else if n_zdata gt 0 then begin
        hop = zdata.color[nrec]
      endif else if nrec gt 0 then begin
        hop = (hop + 1L) mod nhops
      endif else begin
        hop = hop1
      endelse

      ; Remove the DC spike

      dcbin = npts/2L
      yin[dcbin] = 0.5*(yin[dcbin-1] + yin[dcbin+1])

      ; Flip the spectrum (to account for a cabling error) if requested

      if iqerror then yin = reverse(yin)

      ; Shift or interpolate according to the Doppler correction

      if deldopcorr then begin
        if keyword_set(roundephcorr) then begin
          dopbincorr = round(zdata.dopcorr[nrec]/df)
          yin = shift(yin, -dopbincorr)
        endif else begin
          dopbincorr = zdata.dopcorr[nrec]/df
          coords = findgen(npts) + dopbincorr
          yin = interpolate(yin, coords, cubic=-1.0)
        endelse
      endif
      y[hop,*] = yin
      nin = nin + 1L
      nrec = nrec + 1L
    endif
  endfor

  nrec_scan = nrec_scan + nin

  ; Check whether or not to use this hop cycle (or single spectrum):
  ; -- Don't use an incomplete cycle 
  ; -- Use a complete cycle if there's no azel info to say otherwise
  ; -- Don't use a complete cycle if there's azel info which ends
  ;    before the end of the cycle
  ; -- Use a complete cycle with complete azel info only if none of
  ;    the dwells in the cycle have been flagged bad (tx power < 1.0 kW)
  ; -- If you're reducing only one hop color in the cycle, make sure
  ;    that this spectrum has the right color

  if nin ne nin_desired then begin
    use = 0
  endif else if n_zdata eq 0 then begin
    use = 1
  endif else if n_zdata lt nrec then begin
    use = 0
  endif else begin
    use = (total(zdata.pwr[nrec-nin:nrec-1] lt 1.0) eq 0)
  endelse
  if n_elements(hoptouse) gt 0 then begin
    use = use and (hoptouse eq hop or not hopped)
  endif

  ; If the hop cycle is OK then include it in the block in progress

  if use then begin
 
    ; Add the spectra to the sum for this block

    r = r + y

    ; Work on the tags for this block

    if n_zdata gt 0 then begin
      rcsta_min = rcsta_min < ( 1.0D*min(zdata.sec[nrec-nin:nrec-1]) + tzcorr*3600.0D )
      rcend_max = rcend_max > ( 1.0D*max(zdata.sec[nrec-nin:nrec-1]) + dwell $
                                                                     + tzcorr*3600.0D )
      az1_min = az1_min < min(zdata.az1[nrec-nin:nrec-1])
      az1_max = az1_max > max(zdata.az1[nrec-nin:nrec-1])
      az2_min = az2_min < min(zdata.az2[nrec-nin:nrec-1])
      az2_max = az2_max > max(zdata.az2[nrec-nin:nrec-1])
      rtt_min = rtt_min < min(zdata.rtt[nrec-nin:nrec-1])
      rtt_max = rtt_max > max(zdata.rtt[nrec-nin:nrec-1])
      rctime_block = rctime_block + total(zdata.sec[nrec-nin:nrec-1],/double) $
                     + nin*(dwell/2.0D) + nin*tzcorr*3600.0D
      za_block = za_block + total(zdata.za[nrec-nin:nrec-1])
      az1_block = az1_block + total(zdata.az1[nrec-nin:nrec-1])
      az2_block = az2_block + total(zdata.az2[nrec-nin:nrec-1])
      rtt_block = rtt_block + total(zdata.rtt[nrec-nin:nrec-1])
      pwr_block = pwr_block + total(zdata.pwr[nrec-nin:nrec-1])
      tsys_block = tsys_block + total(tsys[nrec-nin:nrec-1])
      gain_block = gain_block + total(zdata.grgt[nrec-nin:nrec-1])
      if radec then begin
        ra1_min = ra1_min < min(zdata.ra1[nrec-nin:nrec-1])
        ra1_max = ra1_max > max(zdata.ra1[nrec-nin:nrec-1])
        ra2_min = ra2_min < min(zdata.ra2[nrec-nin:nrec-1])
        ra2_max = ra2_max > max(zdata.ra2[nrec-nin:nrec-1])
        dec_min = dec_min < min(zdata.dec[nrec-nin:nrec-1])
        dec_max = dec_max > max(zdata.dec[nrec-nin:nrec-1])
        ra1_block = ra1_block + total(zdata.ra1[nrec-nin:nrec-1])
        ra2_block = ra2_block + total(zdata.ra2[nrec-nin:nrec-1])
        dec_block = dec_block + total(zdata.dec[nrec-nin:nrec-1])
      endif
    endif

    ; Work towards the mean value of sdev/rmsc for this block

    sdev_over_rmsc_block = sdev_over_rmsc_block + total(sdev_over_rmsc[nrec-nin:nrec-1])

    ; Increment the number of hop cycles included in this block so far

    iblock = iblock + 1L

  endif

  ; Time to process the block?
  ;
  ; Yes if a full block or else a partial block at the end of the file,
  ; or else (unless merge keyword is set) a partial block at the end of a scan

  fullBlock = iblock eq nblock
  if n_zdata gt 0 then begin
    endOfFile = eof(lun) or (nrec eq n_zdata)
    endOfScan = (endOfFile) ? 1 : (zdata.rec[nrec] eq 0)
  endif else begin
    endOfFile = eof(lun)
    endOfScan = endOfFile
  endelse

endrep until (fullBlock or endOfFile or (endOfScan and not keyword_set(merge)))

; Don't process an empty block at the end of the scan or at the end of the file

if iblock eq 0 then begin
  print,blocknumWithinScan,' -- no usable data in this block',format='(i8,a)'
  return
endif

; Define or extend four arrays: s (processed spectra for all blocks),
; sdev_block (sdev for all blocks), rmsc_block (calculated fractional rms scatter
; for all blocks), and bias_block (bias removed for all blocks).
; By not setting their dimensions in advance we don't have to preset the maximum
; number of blocks allowed.

if nout eq 0 then begin
  s = fltarr(1,nper)
  sdev_block = fltarr(1)
  rmsc_block = fltarr(1)
  bias_block = fltarr(1)
endif else begin
  s = [s, fltarr(1,nper)]
  sdev_block = [sdev_block, 0.0]
  rmsc_block = [rmsc_block, 0.0]
  bias_block = [bias_block, 0.0]
endelse

; Compute two parameters which we'll need to get the average sdev
; for this block

sdev_over_rmsc_block = sdev_over_rmsc_block/(nin_desired*iblock)
tau_block = float(dwell)*nin_desired*iblock

; If dehopping, do background subtraction for this block,
; then do normalization (division by fractional rms noise power)
;
; If the data were unhopped, just transfer the block spectrum to
; the s array for compatibility with dehopping routines

if n_elements(hoptouse) eq 0 then begin
  b = fltarr(nhops,npts)
  for n=0L,nhops-1 do $
    for k=0L,nhops-1 do $
      if k ne n then b[n,*] = b[n,*] + r[k,*]
  b = b/(nhops-1)

  rmsc_spec = sqrt( nhops / ( df*tau_block ) )
  rmsc_bg_raw = sqrt( nhops / ( (nhops-1)*df*tau_block ) )

  ; For each hop color, the ratio of spectrum to background is biased
  ; a bit above unity even in the signal-free region (see header comments
  ; to function ratioexpectation); compute the mean ratio so that the
  ; baselined, normalized block spectra have zero mean, and get the
  ; expected rms deviation about that mean ratio.
  ;
  ; ratioexpectation is somewhat slow, so if we're fitting a polynomial
  ; to the background spectra, we'll just assume a mean of 1 and will
  ; (later) compute the expected rms deviation for each hop color
  ; by adding in quadrature the signal and background rms values

  if nobgfit then begin
    rmsc_bg = rmsc_bg_raw
    if bias_allow then begin
      meanratio = 1
      rms = sqrt(rmsc_spec^2 + rmsc_bg^2)
    endif else begin
      meanratio = ratioexpectation(rmsc_spec, rmsc_bg, rms=rms)
    endelse
  endif else begin
    meanratio = 1   ;  Later we'll compute rms for each hop color
  endelse

  for k=0L,nhops-1 do begin
    spec = reform(r[k,lbin[k]:rbin[k]])
    bg_raw = reform(b[k,lbin[k]:rbin[k]])
    if nobgfit then begin
      bg = bg_raw
    endif else begin

      ; First do an unweighted low-order polynomial fit to get an idea
      ; of the proper baseline; then do a weighted fit assuming measurement
      ; errors which are proportional to the initial fit values.  Note that
      ; setting errors proportional to the *raw* values would put lower
      ; weight on positive noise fluctuations, resulting in a fit with
      ; negative bias.
      ;
      ; In order to use the lowest-degree fit possible, start with degree = 0,
      ; then increase the degree until t-tests show that the extra term isn't
      ; required by the data.  Actually, check for one beyond that, since it
      ; can happen (for example) that a quadratic fit is little improvement
      ; over a linear fit but that a cubic fit then provides significant
      ; improvement over the quadratic fit.
      ;
      ; Under no circumstances use a fit with degree > degreemax, as the risk
      ; of strange behavior in regions excluded from the fit (e.g., near DC)
      ; is too great.

      degree = 0
      repeat begin
        bg_coeffs = polymaskfit(xfit,bg_raw,degree,/spikeremove,yfit=bg_trial,/silent)
        bg_raw_err = bg_trial*rmsc_bg_raw
        bg_coeffs = polymaskfit(xfit,bg_raw,degree,measure_errors=bg_raw_err, $
                                /spikeremove,covar=covar,n_out=n_out,         $
                                yfit=bg_trial,yprederr=bg_sigma_trial,        $
                                probability=prob,/silent)
        if degree eq 0 then begin
          prev_degree_justified = 1
          degree_justified = 1
        endif else if isnull(bg_coeffs[0]) then begin

          ; Fit gave undefined results

          prev_degree_justified = degree_justified
          degree_justified = 0
        endif else begin
          prev_degree_justified = degree_justified
          t_stat = bg_coeffs[degree]/sqrt(covar[degree,degree])
          t_df = n_out - degree
          degree_justified = t_prob(t_df, t_stat) lt 0.05
        endelse
        if degree_justified then begin
          degree_use = degree
          prob_use = prob
          bg = bg_trial
          bg_sigma = bg_sigma_trial
        endif
        if degree_justified or prev_degree_justified then degree = degree + 1
      endrep until (not (degree_justified or prev_degree_justified) or degree gt degreemax)

      ; If desired, show the fit

      if plotfit then begin
        titlestring = chanstring + ' scan ' + string(nscan+1,format='(i0)')            $
                              + ', block ' + string(blocknumWithinScan,format='(i0)')  $
                              + ', hop color ' + string(k,format='(i0)')               $
                              + ':  fit order = ' + string(degree_use,format='(i0)')
        if nper le 10000 then begin
          plot,findgen(nper),bg_raw,title=titlestring,linestyle=2
        endif else begin
          plot,findgen(nper),bg_raw,title=titlestring,psym=3
        endelse
        oplot,findgen(nper),bg,thick=4
        pause,prompt='Press Enter to reduce the next hop: '
      endif

      ; Give a warning if chi-squared is too high for the fit to be credible

      if prob_use lt 0.001 then begin
        print,' '
        print,'WARNING in processRaw_doBlock1 for hop color ',k,' of block ', $
              blocknumWithinScan,':',format='(a,i0,a,i0,a)'
        print,'        The probability that baseline fit (degree = ',degree_use, $
              ') is valid is only ',prob_use,format='(a,i0,a,e9.2)'
        print,' '
      endif

      ; The error on the baseline fit should be so small that it's not worth
      ; calling ratioexpectation to get the rms (or mean) for signal/baseline

      rmsc_bg = mean(bg_sigma/bg)
      rms = sqrt(rmsc_spec^2 + rmsc_bg^2)

    endelse

    ; Now divide signal by background and subtract the expected mean ratio
    ; -- we'll sum over all hop colors, then normalize to unit variance

    weight_color = 1.0/rms^2
    s[nout,*] = s[nout,*] + weight_color*(spec/bg - meanratio)
    bias_block[nout] = bias_block[nout] + weight_color*(meanratio - 1)
    rmsc_block[nout] = rmsc_block[nout] + weight_color

  endfor

  ; Now we can complete the sum over hop colors and normalize to unit variance,
  ; yielding a single dehopped, baselined, normalized spectrum for this block.

  rmsc_block[nout] = 1.0 / sqrt(rmsc_block[nout])    ;  1/sqrt(sum of weights)
  s[nout,*] = s[nout,*]*rmsc_block[nout]
  bias_block[nout] = bias_block[nout]*rmsc_block[nout]

endif else begin

  ; Life's easy if we're not dehopping and baselining....

  rmsc_block[nout] = 1.0/sqrt(df*tau_block)
  s[nout,*] = (hopped) ? r[hoptouse,*] : r[0,*]
  bias_block[nout] = 0

endelse

; Work on various tags, constructing a weighted mean over blocks

sdev_block[nout] = sdev_over_rmsc_block*rmsc_block[nout]
weight_block = 1/sdev_block[nout]^2
nffts = nffts + iblock*nin_desired*floor(srate*dwell/npts + 1d-10)
tau = tau + tau_block
if n_zdata gt 0 then begin
  rctime_mean = rctime_mean + weight_block*(rctime_block/(nin_desired*iblock))
  za_mean = za_mean + weight_block*(za_block/(nin_desired*iblock))
  az1_mean = az1_mean + weight_block*(az1_block/(nin_desired*iblock))
  az2_mean = az2_mean + weight_block*(az2_block/(nin_desired*iblock))
  rtt_mean = rtt_mean + weight_block*(rtt_block/(nin_desired*iblock))
  pwr_mean = pwr_mean + weight_block*(pwr_block/(nin_desired*iblock))
  tsys_mean = tsys_mean + weight_block*(tsys_block/(nin_desired*iblock))
  gain_mean = gain_mean + weight_block*(gain_block/(nin_desired*iblock))
  if radec then begin
    ra1_mean = ra1_mean + weight_block*(ra1_block/(nin_desired*iblock))
    ra2_mean = ra2_mean + weight_block*(ra2_block/(nin_desired*iblock))
    dec_mean = dec_mean + weight_block*(dec_block/(nin_desired*iblock))
  endif
endif

; Compute the average background level for this block
; (avoiding DC if we're not dehopping)

if n_elements(hoptouse) eq 0 then begin
  back = 0.0
  for k=0L,nhops-1 do $
    back = back + total(b[k,lbin[k]:rbin[k]])/(1L*nper*iblock*nhops)
endif else begin
  bin2 = npts/2L - (npts/100L > 1L)
  bin3 = npts/2L + (npts/100L > 1L)
  back = mean( [ [s[nout,0:bin2]], [s[nout,bin3:npts-1]] ] )
endelse

; Display information on this block,
; then increment the number of blocks processed

print,chanstring,blocknumWithinScan,iblock,back,sdev_block[nout],rmsc_block[nout],tau_block, $
      format='(4x,a2,2x,2i8,e13.4,f13.3,f13.4,f13.1)'
nout = nout + 1L
blocknumWithinScan = blocknumWithinScan + 1L

; Reset block parameters so that a new block can be started

iblock = 0L
sdev_over_rmsc_block = 0.0
r = r*0.0
if n_zdata gt 0 then begin
  rctime_block = 0.0D
  za_block = 0.0
  az1_block = 0.0
  az2_block = 0.0
  rtt_block = 0.0
  pwr_block = 0.0
  tsys_block = 0.0
  gain_block = 0.0
  if radec then begin
    ra1_block = 0.0
    ra2_block = 0.0
    dec_block = 0.0
  endif
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro processRaw_sumBlocks1,chan,merge=merge,hoptouse=hoptouse

; Sum up all blocks of raw spectra to create an output spectrum;
; set the tags; and load the spectrum

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common target,tname,diameter,period,jd0,phase0
common setup,srate,hopbw,hopf1,txoff
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common loadedBlock,loadedi,loaded1,loaded
common fileInit,doperr,nper,gconst,tsys,sdev_over_rmsc,yin,y,lun,nrec,endOfScan,endOfFile,nscan
common scanInit,r,nout,iblock,nrec_scan,rcsta_min,rcend_max,az1_min,az1_max,az2_min, $
                az2_max,rctime_block,za_block,az1_block,az2_block,rtt_block,pwr_block, $
                tsys_block,gain_block,rctime_mean,za_mean,az1_mean,az2_mean,rtt_mean, $
                pwr_mean,tsys_mean,gain_mean,nffts,tau,sdev_over_rmsc_block,hop, $
                ra1_min,ra1_max,ra2_min,ra2_max,ra1_block,ra2_block,dec_block, $
                ra1_mean,ra2_mean,dec_mean,dec_min,dec_max,rtt_min,rtt_max,blocknumWithinScan
common blockSpectra,s,sdev_block,rmsc_block,bias_block

monthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
chanstrings=['OC','SC','RE','IM']

; Quit if there are no blocks to sum

if nout eq 0 then begin
  print,'No blocks to sum for this spectrum: no output produced'
  return
endif

; Create the zeroed tag structure to be filled in later

tags1 = blanktags1()
ntags = n_tags(tags1)
tagnames = strlowcase(tag_names(tags1))

; Compute sdev, the radar cross-section equivalent (km^2)
; of the rms noise power of the grand weighted sum

sdev_block = sdev_block[0:nout-1] ; drop extra elements
sumwts = total(1/sdev_block^2)
sdev = sqrt(1.0/sumwts)
tags1.sdev = sdev
print,'Input',nrec_scan,' spectra from ',infile
print,'sdev = ',sdev

; Compute rmsc, the calculated rms fractional noise power

rmsc_block = rmsc_block[0:nout-1] ; drop extra elements
rmsc = sqrt(1.0/total(1/rmsc_block^2))

; Take the weighted sum of the unnormalized block spectra s[k,*]*sdev_block[k].
; Note that multiplying s[k,*]*sdev_block[k] by weight[k] = 1/sdev_block[k]^2
; amounts to having s[k,*]/sdev_block[k] in each term of the sum.  So we have
;  
;       unnormalized sum = ( sum over k of s[k,*]/sdev_block[k] ) / sumwts
;
; Also note that the rms noise in this unnormalized sum is sdev = sqrt(1/sumwts),
; so dividing by this value, and realizing that sdev*sumwts = 1/sdev, we obtain
;
;         normalized sum = ( sum over k of s[k,*]/sdev_block[k] ) * sdev

sum = total( s[0:nout-1,*] / (sdev_block # replicate(1.0,nper)) , 1 ) * sdev

print,'Peak = ',max(sum,m),'  at bin #',m,format='(a,f0,a,i0)'

; Do a similar weighted sum to get the bias subtracted from the spectrum
; (nonzero if /nobgfit was set and /bias_allow wasn't set)

bias_block = bias_block[0:nout-1]
biasremoved = total(bias_block/sdev_block)*sdev

; Define the frequency vector and set the zero-frequency tag
;
; Note that for unhopped data, the frequencies will be a bit off if the 
; tx offset isn't an integer multiple of the frequency resolution

xjcen = (n_elements(hoptouse) gt 0) ? nper/2 + round(txoff/df) : nper/2
f = (findgen(nper) - xjcen)*df
tags1.xjcen = xjcen

; Set some time tags
;
; Note that I'm deviating from the tkplay philosophy and treating the
; weighted mean receive time (rchour rcmin rcsec rcnsec) as the
; important time, one which will be recomputed when this sum is combined
; with others.  Hence I'm forcing the date (iyy, imm, idd) to correspond 
; to that time.  Extra tag jdmean also holds this information.
;
; Tags rcsta and rcend are less important and needn't agree with the date
; once you start summing dehopped spectra from different dates.  However,
; extra tags jdstart and jdend will still be meaningful.

if n_zdata gt 0 then begin

  jdstart = (date - 0.5) + rcsta_min/86400.0
  jdend = (date - 0.5) + rcend_max/86400.0

  ; Have to work at it to round rctime to the nearest nsec, since
  ; IDL 5.3 (Arecibo's version as of Feb 2002) doesn't allow the /L64 keyword
  ; for the round function.

  rctime_mean = rctime_mean/sumwts
  rctime_frac = rctime_mean - floor(rctime_mean)
  rctime_nsec = round(1.0d9*rctime_frac) / 1.0d9
  rctime_mean = floor(rctime_mean) + rctime_nsec

  jdmean = (date - 0.5) + rctime_mean/86400.0
  caldat,floor(jdmean + 0.5),mon,day,year
  old_rctime_mean = rctime_mean
  rctime_mean = rctime_mean - 86400L*floor(rctime_mean/86400)  ;  range [0,86400)
  rctime_shift = round(rctime_mean - old_rctime_mean)
  rcsta_min = rcsta_min + rctime_shift
  rcend_max = rcend_max + rctime_shift
  hr = floor(rctime_mean/3600.0)
  min = floor( (rctime_mean - 3600*hr)/60.0 )
  fsec = rctime_mean - 3600*hr - 60*min
  sec = floor(fsec)
  nsec = round( (fsec - sec)*1.0d9 )
  tags1.rchour = hr
  tags1.rcmin = min
  tags1.rcsec = sec
  tags1.rcnsec = nsec
  tags1.rcsta = rcsta_min
  tags1.rcend = rcend_max

  ; The 'calmean' tag is a string listing the calendar date/time
  ; equivalent of jdmean, rounded to the nearest sec

  caldat_roundsec,jdmean,mon2,day2,year2,hr2,min2,sec2
  calmean = string(year2,format='(i4)') + ' '     $
            + monthnames[mon2-1] + ' '            $
            + string(day2,format='(i2.2)') + ' '  $
            + string(hr2,format='(i2.2)') + ':'   $
            + string(min2,format='(i2.2)') + ':'  $
            + string(sec2,format='(i2.2)') + ' '  $
            + outputTimeZone
endif else begin
  caldat,date,mon,day,year
endelse
tags1.iyy = year
tags1.imm = mon
tags1.idd = day

; Some tags can be set only if there's azel info;
; otherwise leave them at the "reset" value
;
; Set rotation phase to the phase at the mean rx time for the run, modulo 360;
; do NOT assign a phase modulo 360 to each block and then take the mean over blocks.
; Thus a run which covers 20 deg of rotation from 350 deg to 10 deg has its phase
; tag set to about 0.0, not 180.0.
;
; Similarly, use the azimuth range -- [0,360) vs. [180,540) -- which yields the most
; "compact" set of azimuth values.  Then put the result into the range [0,360).

if n_zdata gt 0 then begin
  tags1.elev = 90.0 - za_mean/sumwts
  az1_range = az1_max - az1_min
  az2_range = az2_max - az2_min
  if az2_range lt az1_range then $
    tags1.azim = (az2_mean/sumwts) mod 360 $
  else $
    tags1.azim = az1_mean/sumwts
  tags1.trpwr = pwr_mean/sumwts
  tags1.rttim = rtt_mean/sumwts
  tags1.tsys = tsys_mean/sumwts
  tags1.gain = (gain_mean/sumwts) / 1e12
  phase_mean = (jdmean - jd0)*(24.0/period)*360 + phase0
  tags1.phase = phase_mean - 360*floor(phase_mean/360)
endif

; Set some other tags

tags1.nffts = nffts
tags1.doppl = txoff
tags1.lfft = npts
tags1.igw = 1.0e6/float(srate)
tags1.dfreq = df
tags1.tau = tau
tags1.rmsc = rmsc
tags1.kpts = nper
tags1.jcp = long(chan)
tags1.nfreq = nhops
tags1.frstep = hopbw
tags1.color = (n_elements(hoptouse) gt 0) ? hoptouse : -1L
tags1.freq1 = hopf1
if notnull(filesuffix) then begin
  tags1.irun = long(strmid(filesuffix,1))
endif else if keyword_set(merge) and nscan gt 1 then begin
  tags1.irun = 0L
endif else begin
  tags1.irun = nscan
endelse

; Create extra tags which go in an rdf footer.
; We'll start here with the names and numbers of the regular cw tags;
; after we've loaded the spectrum we'll use addextra1 to add some more.

nextra = ntags
extratag = {format:'t', name:'', value:'', comment:''}
extratags = replicate(extratag, nextra)
extratags.name = tagnames
extratags.value = string(indgen(ntags),format='(i0)')

; Define a string that identifies the polarization

chanstring = chanstrings[chan-1]

; Put everything into a structure and load it

chanStruc = {freq:f, spec:sum, tags:tags1, extratags:extratags, $
             ndata:nper, ntags:ntags, nextra:nextra, tname:tname, $
             pol:chanstring}

*loaded1 = chanStruc

; Add some other extra tags

addextra1,'f','diameter',strtrim(diameter, 2)
addextra1,'f','period',strtrim(period, 2)
addextra1,'f','lambda',string(lambda,format='(f9.7)')
addextra1,'s','xmit_sta','Arecibo'
addextra1,'i','tzcorr',tzcorr
if notnull(outputTimeZone) then addextra1,'s','timezone',outputTimeZone
if n_zdata gt 0 then begin
  addextra1,'f','phase0',strtrim(phase0, 2)
  addextra1,'d','jd0',string(jd0,format='(d13.5)')
  addextra1,'d','jdstart',string(jdstart,format='(d13.5)')  ;  1 sec precision
  addextra1,'d','jdmean',string(jdmean,format='(d13.5)')
  addextra1,'s','calmean',calmean
  addextra1,'d','jdend',string(jdend,format='(d13.5)')
  rtt_1AU = 998.0095670D0   ;   seconds
  addextra1,'f','distmin',string((rtt_min/rtt_1AU),format='(f7.5)')
  addextra1,'f','distmean',string(( (rtt_mean/sumwts) / rtt_1AU ),format='(f7.5)')
  addextra1,'f','distmax',string((rtt_max/rtt_1AU),format='(f7.5)')
  if radec then begin
    ra1_range = ra1_max - ra1_min
    ra2_range = ra2_max - ra2_min
    if ra2_range lt ra1_range then begin
      addextra1,'f','ramin',string((ra2_min mod 24),format='(f8.5)')
      addextra1,'f','ramean',string(((ra2_mean/sumwts) mod 24),format='(f8.5)')
      addextra1,'f','ramax',string((ra2_max mod 24),format='(f8.5)')
    endif else begin
      addextra1,'f','ramin',string(ra1_min,format='(f8.5)')
      addextra1,'f','ramean',string((ra1_mean/sumwts),format='(f8.5)')
      addextra1,'f','ramax',string(ra1_max,format='(f8.5)')
    endelse
    addextra1,'f','decmin',string(dec_min,format='(f8.4)')
    addextra1,'f','decmean',string((dec_mean/sumwts),format='(f8.4)')
    addextra1,'f','decmax',string(dec_max,format='(f8.4)')
  endif
endif
if n_elements(hoptouse) eq 0 then addextra1,'f','biasremoved',biasremoved

; Now that some target parameters have been placed in the extra tags,
; use the siglim1 procedure to set the signal limit (jsnr1, jsnr2) tags
; according to the target bandwidth.

siglim1,0,chanstring=chanstring
lim1 = (*loaded1).tags.jsnr1
lim2 = (*loaded1).tags.jsnr2

; Last step:
;
; For hopped data, we can now measure the bias and rms of the noise baseline.
; -- Print the bias (in case dehopping didn't work and the baseline
;    removal was inadequate)
; -- Print the measured rms noise and also set the rmsm tag
;    (rmsm = measured rms noise divided by calculated rms noise)

if n_elements(hoptouse) eq 0 then begin
  ndata = (*loaded1).ndata
  in_mask = intarr(ndata) + 1
  in_mask[lim1:lim2] = 0
  n_noise = ndata - (lim2 - lim1 + 1L)
  if n_noise ge 2 then begin
    meanNoise = polymaskfit(f,sum,0,in_mask=in_mask,yerror=rmsNoise,/silent)
    print,'Measured/calculated noise: bias = ',meanNoise,' and rms = ', $
          rmsNoise,'  (',n_noise,' points)',format='(a,f6.3,a,f5.3,a,i0,a)'
    print,' '
    settag1,'rmsm',rmsNoise,/silent
  endif else begin
    print,"Entire spectrum is signal: Can't measure noise or set rmsm tag"
    print,' '
  endelse
endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro processRaw_combineChans,n_OC,n_SC,nstack1_start,n_RE,n_IM

; After processing hopped or unhopped date, decide what to do with the
; OC and SC output spectra sitting in the single-channel stack, then do it

common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

npol = n_params()-1

nmin = n_OC < n_SC
nmax = n_OC > n_SC
if npol gt 2 then begin
  nmin = nmin < n_RE
  nmax = nmax > n_RE
endif
if npol gt 3 then begin
  nmin = nmin < n_IM
  nmax = nmax > n_IM
endif

if n_OC eq 0 and n_SC eq 0 then begin

  ; Problem: No spectra were output!

  print,'No output spectra were produced'
  print,' '

endif else if nmin ne nmax then begin

  ; Problem:  Somehow there are more OC than SC output spectra (or vice versa)
  ; Solution: Just leave them all in the single-channel stack

  print,'Added ',n_OC,' OC and ',n_SC,' SC spectra to the single-channel stack,', $
        ' which now contains ',nstack1,' spectra',format='(a,i0,a,i0,2a,i0,a)'
  load1,(nstack1_start+1)
  print,'Loaded ',(*loaded1).pol,' spectrum #1 (stack1 #',nstack1_start+1,')', $
        format='(3a,i0,a)'
  print,' '

endif else begin

  ; It worked: Equal numbers of OC and SC output spectra, so combine them

  for k=0L,n_OC-1 do begin

    ; Splice together an OC/SC pair

    if npol eq 2 then splice,(nstack1_start + k + 1),(nstack1_start + k + 1 + n_OC),/keep,/silent
    if npol eq 3 then splice,(nstack1_start + k + 1),(nstack1_start + k + 1 + n_OC),(nstack1_start + k + 1 + 2*n_OC),/keep,/silent
    if npol eq 4 then splice,(nstack1_start + k + 1),(nstack1_start + k + 1 + n_OC),(nstack1_start + k + 1 + 2*n_OC),(nstack1_start + k + 1 + 3*n_OC),/keep,/silent
    

    ; If the extra tag for mean Julian date is set, just stick with the OC value.
    ; The two dates could differ to the extent that the two channels have different
    ; fractional variations in system temperature as a function of time, but such
    ; differences are unlikely to be large enough to be worth recording.
    ;
    ; Ditto for mean distance, RA, dec, and bias removed

    deleteextra,'jdmean_SC',/silent
    deleteextra,'calmean_SC',/silent
    deleteextra,'distmean_SC',/silent
    deleteextra,'ramean_SC',/silent
    deleteextra,'decmean_SC',/silent
    deleteextra,'biasremoved_SC',/silent
    setextraname,'jdmean_OC','jdmean',/silent
    setextraname,'calmean_OC','calmean',/silent
    setextraname,'distmean_OC','distmean',/silent
    setextraname,'ramean_OC','ramean',/silent
    setextraname,'decmean_OC','decmean',/silent
    setextraname,'biasremoved_OC','biasremoved',/silent
    if npol gt 2 then begin
      deleteextra,'jdmean_RE',/silent
      deleteextra,'calmean_RE',/silent
      deleteextra,'distmean_RE',/silent
      deleteextra,'ramean_RE',/silent
      deleteextra,'decmean_RE',/silent
      deleteextra,'biasremoved_RE',/silent
    endif
    if npol gt 3 then begin
      deleteextra,'jdmean_IM',/silent
      deleteextra,'calmean_IM',/silent
      deleteextra,'distmean_IM',/silent
      deleteextra,'ramean_IM',/silent
      deleteextra,'decmean_IM',/silent
      deleteextra,'biasremoved_IM',/silent
    endif

    ; Push the new spectrum onto the pair stack
    ; (especially important if there are more spectra to come)

    push,/silent

  endfor

  ; Now delete the single-channel spectra in stack1
  ; (Can't do it earlier: Deleting one spectrum changes the numbers of the others)

  for k=nstack1,nstack1_start+1,-1 do deletestack1,n=k,/silent

  ; Load pair #1, and report results to user

  print,'Added ',n_OC,' OC/SC pairs to the stack, which now contains ', $
        nstack,' spectra',format='(a,i0,a,i0,2a,i0,a)'
  load,(nstack - n_OC + 1)
  print,'Loaded pair #1 (stack #',nstack-n_OC+1,')',format='(a,i0,a)'
  print,' '

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro nohop_start,hoptouse=hoptouse,flipchan=flipchan,nostokes=nostokes,help=help,_extra=_ext

; Start processing an OC/SC spectral pair which wasn't hopped,
; or else a hopped pair for which you want the raw weighted sum for just one hop color
; If there are 4 files, do Stokes unless told not to.
;
; Produce a weighted spectral sum, which afterwards can be vignetted using vignette,
; then baseline-subtracted and normalized using nohop_finish

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common loadedBlock,loadedi,loaded1,loaded
common stackBlock,stacki,stack1,stack,nstacki,nstack1,nstack

if n_params() ne 0 or keyword_set(help) then begin
  print,' '
  print,'nohop_start[,hoptouse=hoptouse][,/noweight][,/merge OR ,/blocks]  $'
  print,'           [,/roundephcorr][,flipchan=flipchan][,/help]'
  print,' '
  print,'      hoptouse lets you get the raw weighted sum for ', $
        'just one hop color (start @ 0) of a hopped pair',format='(2a)'
  print,'           (-1 or omitted --> use all data)'
  print,'      /merge produces one output spectrum for the entire file, ', $
        'rather than one per scan',format='(2a)'
  print,'      /blocks produces one output spectrum per block rather than one per scan'
  print,'      /roundephcorr rounds Doppler ephemeris corrections so that the spectra'
  print,'            are shifted by an integer number of bins; otherwise no rounding'
  print,'            is carried out and interpolation (cubic convolution) is performed'
  print,'            rather than simple shifting'
  print,'      flipchan is used to correct a cabling error (switching I and Q):'
  print,'            flipping is done about the middle of the spectrum, and neither'
  print,'            xjcen, posfr, nor the frequency vector are changed.  Flipping is'
  print,'            done before any shifting is carried out due to /roundephcorr.'
  print,'            Permitted values are flipchan = 1, 2, or 3 (flip both channels).'
  print,'            Note: for Stonkes flipchan must be 0 or 3 if npol=4: you can''t'
  print,'            flip the crosses'
  print, '     /nostokes to skip Stokes processing even if .p3 and .p4 files are there'
  print,' '
  return
endif

if (notnull(zstem)) and (isnull(filesuffix)) then begin
  print,"ERROR in nohop_start: You forgot to set fsuf, so zfile isn't specified"
  return
endif

; Validate the hop color to be used if specified, or else set it to -1 (all)

if nhops le 1 then begin
  if nhops lt 1 then print,'Setting nhops = 1, hop1 = 0, assuming unhopped data'
  nhops = 1
  hop1 = 0
  hoptouse = -1
endif else if n_elements(hoptouse) gt 0 then begin
  if hoptouse lt -1 or hoptouse ge nhops then begin
    print,'ERROR in nohop_start: hoptouse must be in the range [0, ',nhops-1,'],', $
          format='(a,i0,a)'
    print,'                      or else -1 (or omit) to use all data'
    return
  endif
endif else begin
  hoptouse = -1
endelse

; Deal with the flipchan flag

if not keyword_set(flipchan) then begin
  iqerror_OC = 0
  iqerror_SC = 0
  flipchan = 0
endif else begin
  iqerror_OC = (flipchan eq 1 or flipchan eq 3) ? 1 : 0
  iqerror_SC = (flipchan eq 2 or flipchan eq 3) ? 1 : 0
endelse

; Find out how many spectra are in the single-channel stack at the outset

nstack1_start = nstack1

; Store the loaded single-channel spectrum so we can use that "space" and
; reload the spectrum later on

storeloaded1 = ptr_new(*loaded1)

; Process each polarization in turn, each time setting the datafile name
; to the appropriate value ('.p1' or '.p2' suffix for OC, SC, respectively).
; If necessary set the azel filename as well.
;
; Start with the OC data; the output spectra will be added to the
; single-channel stack

infile = filestem + '.p1' + filesuffix
if notnull(filesuffix) then zfile = zstem + filesuffix
nohop_start1,1,hoptouse=hoptouse,iqerror=iqerror_OC,/silent,_extra=_ext
n_OC = nstack1 - nstack1_start

; Process the SC data

infile = filestem + '.p2' + filesuffix
nohop_start1,2,hoptouse=hoptouse,iqerror=iqerror_SC,/silent,_extra=_ext
n_SC = nstack1 - (nstack1_start + n_OC)
npol = 2

if ~keyword_set(nostokes) then begin
; Look for and then Dehop the RE data

; iqerror_OC and SC are the same, can send in either.
  infile = filestem + '.p3' + filesuffix
  if file_test(infile, /read) then begin
    if flipchan ne 3 && flipchan ne 0 then begin
      print, "flipchan different for pols, skipping Stokes processing"
    endif else begin
      nohop_start1,3,hoptouse=hoptouse,iqerror=iqerror_SC,/silent,_extra=_ext
      n_RE = nstack1 - (nstack1_start + n_SC + n_OC)
      npol = 3

; Look for and Dehop the IM data

      infile = filestem + '.p4' + filesuffix
      if file_test(infile, /read) then begin
        nohop_start1,4,hoptouse=hoptouse,iqerror=iqerror_SC,/silent,_extra=_ext
        n_IM = nstack1 - (nstack1_start + n_RE + n_SC + n_OC)
        npol = 4
      endif
    endelse
  endif
endif


; Decide what to do with the two channels, then do it

if npol eq 3 then begin 
  print, "ERROR in Dehop: Found 3 pols, will process as OC/SC only"
  npol = 2
endif
if npol eq 2 then processRaw_combineChans,n_OC,n_SC,nstack1_start
if npol eq 4 then processRaw_combineChans,n_OC,n_SC,nstack1_start,n_RE,n_IM

; If everything worked correctly and all output channels were spliced into
; OC/SC pairs, reload the single-channel spectrum that was there at the start

if npol eq 2 then if n_OC eq n_SC then *loaded1 = *storeloaded1 else $
if npol eq 3 then if n_OC eq n_SC and n_OC eq n_RE then *loaded1 = *storeloaded1 else $
if npol eq 4 then if n_OC eq n_SC and n_OC eq n_RE and n_OC eq n_IM then *loaded1 = *storeloaded1

infile = ''
if notnull(filesuffix) then zfile = ''
ptr_free,storeloaded1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro nohop_start1,chan,hoptouse=hoptouse,noweight=noweight,merge=merge,blocks=blocks, $
                 roundephcorr=roundephcorr,iqerror=iqerror,silent=silent,help=help

; Start processing a single-channel spectrum which wasn't hopped,
; or else one hop color of a hopped single-channel spectrum
;
; Produce a weighted sum, then load the result
;
; Vignetting, baseline subtraction, and normalization are left for later
; processing (via procedures vignette1 and nohop_finish1)
;
; Modified by CM 2003 Dec 13: Check whether nout > 0 before pushing the spectrum onto the
;          single-channel stack, so that we don't get the last complete block pushed twice
;          when /blocks is set and only unused data follow this block

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda
common target,tname,diameter,period,jd0,phase0
common setup,srate,hopbw,hopf1,txoff
common azel,zfile,zstem,zdata,n_zdata,tzcorr,outputTimeZone,radec,deldopcorr
common loadedBlock,loadedi,loaded1,loaded
common fileInit,doperr,nper,gconst,tsys,sdev_over_rmsc,yin,y,lun,nrec,endOfScan,endOfFile,nscan
common scanInit,r,nout,iblock,nrec_scan,rcsta_min,rcend_max,az1_min,az1_max,az2_min, $
                az2_max,rctime_block,za_block,az1_block,az2_block,rtt_block,pwr_block, $
                tsys_block,gain_block,rctime_mean,za_mean,az1_mean,az2_mean,rtt_mean, $
                pwr_mean,tsys_mean,gain_mean,nffts,tau,sdev_over_rmsc_block,hop, $
                ra1_min,ra1_max,ra2_min,ra2_max,ra1_block,ra2_block,dec_block, $
                ra1_mean,ra2_mean,dec_mean,dec_min,dec_max,rtt_min,rtt_max,blocknumWithinScan

if not keyword_set(merge) then merge = 0
if not keyword_set(blocks) then blocks = 0
if not keyword_set(roundephcorr) then roundephcorr = 0
if not keyword_set(iqerror) then iqerror = 0

if n_params() ne 1 or keyword_set(help) then begin
  print,' '
  print,'nohop_start1,chan[,hoptouse=hoptouse][,/noweight][,/merge OR ,/blocks]  $'
  print,'                 [,/roundephcorr][,/iqerror][,/silent][,/help]'
  print,' '
  print,'       hoptouse lets you get the raw weighted sum for ', $
        'just one hop color (start @ 0) of a hopped pair',format='(2a)'
  print,'           (-1 or omitted --> use all data)'
  print,'       /merge produces one output spectrum for the entire file, ', $
        'rather than one per scan',format='(2a)'
  print,'       /blocks produces one output spectrum per block rather than one per scan'
  print,'       /roundephcorr rounds Doppler ephemeris corrections so that the spectrum'
  print,'             is shifted by an integer number of bins; otherwise no rounding'
  print,'             is carried out and interpolation (cubic convolution) is performed'
  print,'             rather than simple shifting'
  print,'       /iqerror is used to correct a cabling error (switching I and Q):'
  print,'             the spectrum is flipped about its middle, and neither xjcen,'
  print,'             posfr, nor the frequency vector are changed.  Flipping is'
  print,'             done before any shifting is carried out due to /roundephcorr.'
  print,' '
  return
endif else if chan lt 1 or chan gt 4 then begin
  print,' '
  print,'nohop_start1,chan[,hoptouse=hoptouse][,/noweight][,/merge OR ,/blocks]  $'
  print,'                 [,/roundephcorr][,/iqerror][,/silent][,/help]'
  print,'ERROR in nohop_start1: Must have chan = 1 (OC) or 2 (SC)'
  print,' '
  return
endif else if merge and blocks then begin
  print,' '
  print,'nohop_start1,chan[,hoptouse=hoptouse][,/noweight][,/merge OR ,/blocks]  $'
  print,'                 [,/roundephcorr][,/iqerror][,/silent][,/help]'
  print,"ERROR in nohop_start1: Can't set both /merge and /blocks"
  print,' '
  return
endif

; Validate the hop color to be used if specified, or else set it to -1 (all)

if nhops le 1 then begin
  if nhops lt 1 then print,'Setting nhops = 1, hop1 = 0, assuming unhopped data'
  nhops = 1
  hop1 = 0
  hoptouse = -1
endif else if n_elements(hoptouse) gt 0 then begin
  if hoptouse lt -1 or hoptouse ge nhops then begin
    print,'ERROR in nohop_start1: hoptouse must be in the range [0, ',nhops-1,'],', $
          format='(a,i0,a)'
    print,'                       or else -1 (or omit) to use all data'
    return
  endif
endif else begin
  hoptouse = -1
endelse

; Initialize various quantities before starting work

processRaw_fileInit1,chan,hoptouse=hoptouse,noweight=noweight

; Start reading spectra

repeat begin

  ; Initialize the next output spectrum

  processRaw_scanInit1

  ; Process the blocks of raw spectra which go into that output spectrum

  repeat begin
    processRaw_doBlock1,chan,merge=merge,hoptouse=hoptouse, $
                        roundephcorr=roundephcorr,iqerror=iqerror
  endrep until (endOfFile or blocks or (endOfScan and not merge))

  ; If any "good" blocks were processed, create the finished spectrum

  if nout gt 0 then begin

    ; Sum the blocks, set the tags, and load the output spectrum

    processRaw_sumBlocks1,chan,merge=merge,hoptouse=hoptouse

    ; Push the new spectrum onto the single-channel stack
    ; (especially important if there are more spectra to come)

    push1,silent=silent

  endif

endrep until (endOfFile)

; All raw spectra have been read in

close,lun
free_lun,lun

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro nohop_finish,degree,omitf=farray,help=help

; Finish processing an OC/SC spectral pair which wasn't hopped:
;
; Remove a polynomial baseline from the (presumably vignetted)
; summed spectrum output by nohop_start, then normalize by dividing
; by the baseline and by the value of the rmsc tag.
;
; Modified 8/9/02 to add the omitf=farray keyword option for omitting
; frequency ranges from the baseline fit (in addition to the signal
; region); this was spurred by the need to omit the imperfectly
; interpolated DC region when baselining some unhopped spectra.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded

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

if n_params() eq 0 then degree = -1  ;  just so it's defined

if n_params() ne 1 or degree lt 0 or not omitOK or keyword_set(help) then begin
  print,' '
  print,'nohop_finish,degree[,omitf=farray][,/help]'
  print,' '
  print,'        degree specifies a polynomial baseline'
  print,' '
  print,'        farray is an array of frequency pairs defining regions'
  print,'             (IN ADDITION to the signal region) which should be omitted'
  print,'             from the baseline fit.  Values are in Hz and must be'
  print,'             in increasing order even if the spectrum has frequency'
  print,'             increasing leftward.  Examples of farray are [-150,-120]'
  print,'             or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  return
endif

ndata = (*loaded).ndata
if ndata le 2 then begin
  print,'ERROR in nohop_finish: No pair is loaded'
  return
endif

; Do background subtraction and normalization

if n_elements(farray) gt 0 then begin
  bline,degree,omitf=farray,/divide
endif else begin
  bline,degree,/divide
endelse
pair = (*loaded).spec / ( (*loaded).tags.rmsc # replicate(1.0, ndata) )
(*loaded).spec = pair

; Mask out the frequency bin ranges which were omitted from the fit
; (other than the signal range, which could vary by polarization)

freq = (*loaded).freq
df = (*loaded).tags[0].dfreq
posfr = (*loaded).tags[0].posfr
bin0 = (*loaded).tags[0].xjcen
mask_partial = intarr(ndata) + 1
for k=0L,n_omitregions-2 do begin
  if posfr eq 1 then begin
    omit_left = bin0 + round(farray[0,k]/df)
    omit_right = bin0 + round(farray[1,k]/df)
  endif else begin
    omit_left = bin0 - round(farray[1,k]/df)
    omit_right = bin0 - round(farray[0,k]/df)
  endelse
  if omit_left lt ndata and omit_right ge 0 then begin
    omit_left = omit_left > 0
    omit_right = omit_right < (ndata - 1)
    mask_partial[omit_left:omit_right] = 0
  endif
endfor

; Print the measured rms noise and also set the rmsm tag
; (rmsm = measured rms noise divided by calculated rms noise)

npol = n_elements((*loaded).spec[*,0])

for ch=1,npol do begin

  ; Mask out the signal range for this channnel

  in_mask = mask_partial
  omit_left = (*loaded).tags[ch-1].jsnr1
  omit_right = (*loaded).tags[ch-1].jsnr2
  in_mask[omit_left:omit_right] = 0
  n_noise = round(total(in_mask))

  if n_noise ge 2 then begin
    spec = reform((*loaded).spec[ch-1,*])
    specmean = polymaskfit(freq,spec,0,in_mask=in_mask,yerror=rmsNoise,/silent)
    print,'Channel ',ch,': measured/calculated noise = ',rmsNoise,format='(a,i0,a,f5.3)'
    (*loaded).tags[ch-1].rmsm = rmsNoise
  endif else begin
    print,"Channel ",ch,": no baseline points, so can't measure noise or set rmsm tag", $
          format='(a,i0,a)'
  endelse
endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro nohop_finish1,degree,omitf=farray,help=help

; Finish processing an single-channel spectrum which wasn't hopped:
;
; Remove a polynomial baseline from the (presumably vignetted)
; summed spectrum output by nohop_start1, then normalize by dividing
; by the baseline and by the value of the rmsc tag.
;
; Modified 8/9/02 to add the omitf=farray keyword option for omitting
; frequency ranges from the baseline fit (in addition to the signal
; region); this was spurred by the need to omit the imperfectly
; interpolated DC region when baselining some unhopped spectra.
;
; 2006 Jun 25: frequency refers to center of bin, not left edge

common loadedBlock,loadedi,loaded1,loaded

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

if n_params() eq 0 then degree = -1  ;  just so it's defined

if n_params() ne 1 or degree lt 0 or keyword_set(help) then begin
  print,' '
  print,'nohop_finish1,degree[,omit=farray][,/help]'
  print,' '
  print,'        degree specifies a polynomial baseline'
  print,' '
  print,'        farray is an array of frequency pairs defining regions'
  print,'             (IN ADDITION to the signal region) which should be omitted'
  print,'             from the baseline fit.  Values are in Hz and must be'
  print,'             in increasing order even if the spectrum has frequency'
  print,'             increasing leftward.  Examples of farray are [-150,-120]'
  print,'             or [[-150,-120], [50,90.8], [205,215]].'
  print,' '
  return
endif

ndata = (*loaded1).ndata
if ndata le 2 then begin
  print,'ERROR in nohop_finish1: No single-channel spectrum is loaded'
  return
endif

; Do background subtraction and normalization

if n_elements(farray) gt 0 then begin
  bline1,degree,omitf=farray,/divide
endif else begin
  bline1,degree,/divide
endelse
spec = (*loaded1).spec / ( (*loaded1).tags.rmsc # replicate(1.0, ndata) )
(*loaded1).spec = spec

; Mask out the frequency bin ranges which were omitted from the fit

df = (*loaded1).tags.dfreq
posfr = (*loaded1).tags.posfr
bin0 = (*loaded1).tags.xjcen
in_mask = intarr(ndata) + 1
for k=0L,n_omitregions-2 do begin
  if posfr eq 1 then begin
    omit_left = bin0 + round(farray[0,k]/df)
    omit_right = bin0 + round(farray[1,k]/df)
  endif else begin
    omit_left = bin0 - round(farray[1,k]/df)
    omit_right = bin0 - round(farray[0,k]/df)
  endelse
  if omit_left lt ndata and omit_right ge 0 then begin
    omit_left = omit_left > 0
    omit_right = omit_right < (ndata - 1)
    mask_partial[omit_left:omit_right] = 0
  endif
endfor
omit_left = (*loaded1).tags.jsnr1
omit_right = (*loaded1).tags.jsnr2
in_mask[omit_left:omit_right] = 0
n_noise = round(total(in_mask))

; Print the measured rms noise and also set the rmsm tag
; (rmsm = measured rms noise divided by calculated rms noise)

if n_noise ge 2 then begin
  freq = (*loaded1).freq
  spec = (*loaded1).spec
  specmean = polymaskfit(freq,spec,0,in_mask=in_mask,yerror=rmsNoise,/silent)
  print,'Measured/calculated noise = ',rmsNoise,format='(a,f5.3)'
  (*loaded1).tags.rmsm = rmsNoise
endif else begin
  print,"No baseline points, so can't measure noise or set rmsm tag"
endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rmdc

; Uses interpolation to remove the DC value in a raw pair

common loadedBlock,loadedi,loaded1,loaded

ndata = (*loaded).ndata
if ndata le 2 then begin
  printf,'ERROR in rmdc: No pair is loaded'
  return
endif

pair = (*loaded).spec
n0 = ndata/2L
(*loaded).spec[*,n0] = (pair[*,n0-1] + pair[*,n0+1])/2

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rmdc1

; Uses interpolation to remove the DC value in a raw single-channel spectrum

common loadedBlock,loadedi,loaded1,loaded

ndata = (*loaded1).ndata
if ndata le 2 then begin
  printf,'ERROR in rmdc1: No single-channel spectrum is loaded'
  return
endif

spec = (*loaded1).spec
n0 = ndata/2L
(*loaded1).spec[n0] = (spec[n0-1] + spec[n0+1])/2

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rawsum,s,mean,ave=ave,help=help

; Output the unweighted mean of all raw spectra in the input file
;
; Also average ave spectra at a a time (default = 1) and output a vector
; containing the mean spectral values for the central 50% of the frequency
; range (but with DC excluded)

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda

if n_params() ne 2 or keyword_set(help) then begin
  print,' '
  print,'rawsum,spec,mean[,ave=ave,/help]'
  print,'    spec is the unweighted mean of all raw spectra in infile'
  print,'    mean is a vector containing the mean spectral values for the central', $
                 ' 50% of the frequency range (excluding DC)',format='(2a)'
  print,'    -- mean is computed for the unweighed mean of ave spectra at a time', $
                 ' (default = 1)',format='(2a)'
  print,' '
  return
endif

if n_elements(ave) eq 0 then ave = 1L

; Open the input file

err = openInfile(lun,infile,/get_lun)
if err ne 0 then return

; Initialize parameters

i = 0L
j = 0L
nspec = 0L
y = fltarr(npts)
s = fltarr(npts)
mean = fltarr(1)
bin1 = npts/4L
bin2 = npts/2L - (npts/100L > 1L)
bin3 = npts/2L + (npts/100L > 1L)
bin4 = 3L*npts/4
nbins = (bin2 - bin1 + 1) + (bin4 - bin3 + 1)

; Read the spectra: i counts the number of means, j counts within each mean

while not eof(lun) do begin
  readu,lun,y
  y = swap_endian(temporary(y), /swap_if_little_endian)
  s = s + y
  y=shift(y,npts/2L)
  mean[i] = mean[i] + ( total(y[bin1:bin2]) + total(y[bin3:bin4]) ) / nbins
  nspec = nspec + 1L
  j = j + 1L
  if j eq ave or eof(lun) then begin 
    mean[i] = mean[i]/j
    mean = [mean, 0.0]
    i = i + 1L
    j = 0L
  endif
endwhile

close,lun
free_lun,lun

; Report the results

mean = transpose(mean[0:i-1])
print,'Input ',nspec,' spectra from ',infile,'; output 1 mean spectrum and ',i, $
      ' mean spectral values',format='(a,i0,3a,i0,a)'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro getspec,s,n,removedc=removedc,help=help

; Read a single raw spectrum from the input file

common param,infile,filestem,filesuffix,npts,nhops,hop1,lbin,rbin,df,dwell,date,nblock,lambda

if n_params() ne 2 or keyword_set(help) then begin
  print,'getspec,spec,n[,/removedc,/help]'
  print,'    n is the number (starting from 1) of the desired raw spectrum'
  return
endif

; Open the input file

err = openInfile(lun,infile,/get_lun)
if err ne 0 then return

; Read spectra until you get to the one requested (or else the last one in the file)

i = 0L
s = fltarr(npts)
while not eof(lun) and i ne n do begin
  readu,lun,s
  s = swap_endian(temporary(s), /swap_if_little_endian)
  i = i + 1L
endwhile

close,lun
free_lun,lun

; Report which spectrum was obtained; if requested, remove DC as well

printstring = 'Read spectrum #' + string(i,format='(i0)')
if i lt n then printstring = printstring + ' (last spectrum in file)'
if keyword_set(removedc) then begin
  rmdc,s
  printstring = printstring + ' and removed DC'
endif
print,printstring

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
