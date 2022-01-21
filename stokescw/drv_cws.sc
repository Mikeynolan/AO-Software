#! /bin/tcsh -f
# new usage:
# drv_cw.sc controlfile firstscan lastscan
#
# old usage:
# drv_cw.sc fileNumfirstfile numfiles
#	..   
#	process cw information.. input setup info from $DRVSBCW file
# history:
#
# 19oct21   MCN adding in simple stokes processing. May be a sneaky way, but I'm
#           just going to make temp files.
#
# 14aug19 - PJP copied version from /pkg/aosft/common/bin and then replaces
#           keyval.sc< call with newer keyval.sc ..
# 10dec12 - PT updated all scripts to Linux-friendly versions in /pkg/aosoft,
#           otherwise same as ~cmagri/bin/drv_cw.sc
#
# 17jul07 - CM added "zerofill" parameter that (if set to "yes" or "true" or "1")
#           tells the zerofill routine to zero-fill the final fraction of a
#           transform rather than throwing it out.  This can be used to obtain
#           sufficient frequency resolution when the number of data points is
#           less than the fft length.  The default value is "no".
#
# 17oct04 - CM allowed fsuf values which are hyphenated numerical ranges
#           (such as "02-13"); if the firstscan and lastscan arguments are
#           omitted, just that range of scans is processed.
#
#           CM also allowed fsuf values which are neither a single number nor a
#           hyphenated numerical range (such as "3a" or "test"); if the firstscan
#           and lastscan arguments are omitted, just scan 1 is processed.
#
# 15feb04 - CM allowed fsuf values with leading zeroes (such as "03") for situations
#           when there are ten or more output files and you want an ls listing
#           to show them in the proper order.
#
# 12feb02 - CM added -n flag to the zerofill call (already there in Phil's version),
#           to throw out the final fraction of a transform in a scan rather than
#           zero-filling it.  This stops spurious extra spectra from being output.
#
# 27sep01 - CM added "fsuf" parameter to give output files an extra suffix (when
#           different scans in the same datafile need to be processed differently).
#           If fsuf is specified and no start or end scan numbers are given on the
#           command line, just the scan specified by fsuf is processed (instead of
#           just scan 1 as before).  One less thing to forget and mess up.
#
# 12jan01 - CM added "samplestoadd" option from Mike Nolan's version: can add a few
#           zero records to the start of an incorrectly received dwell (i.e., data
#           before 2000 Sep 10 UTC) rather than skipping most of dwell
#
# jan01 - added back in: multiple scan processing.  all scans output in same file.
#         added back in: numskip (why did I remove it anyway?)
#         also: uses my version of fftfilter to incorporate numskip
# nov99 - modified several options for my own use. uses my version of avgdata.
#         removed command line numbers and added controlfile - G. Black
#
#set verbose
#
set stem = "/pkg/aosoft/fedora4/x86_64/bin"
set parms=($*)              
if ( $#parms == 0 || $#parms > 3 ) then
  echo "Usage: drv_cw.sc <control file> fstart fend"
  exit(0)
endif 
set DRVSBCW = $1
set fstart = 1
set fend = 1
if ( $#parms > 1) then
  set fstart = $2
  set fend = $2
endif
if ( $#parms > 2) then
  set fend = $3
endif
#
set basefile=`keyval.sc< $DRVSBCW fbase`
set filesuffix = `keyval.sc< $DRVSBCW fsuf`
set inpfile=`keyval.sc< $DRVSBCW inpfile`
set bits=`keyval.sc< $DRVSBCW bits`
set fftlen=`keyval.sc< $DRVSBCW fftlen`
set numpol=`keyval.sc< $DRVSBCW numpol`
set numchan = $numpol
if ($numchan == 2) set numchan = 4 # Stokes
set firstrec=`keyval.sc< $DRVSBCW firstrec`
set lastrec=`keyval.sc< $DRVSBCW lastrec`
set ignorefirst=`keyval.sc< $DRVSBCW ignorefirst`
set nfftave=`keyval.sc< $DRVSBCW nfftave`
set numskip=`keyval.sc< $DRVSBCW numskip`
set samplestoadd=`keyval.sc< $DRVSBCW samplestoadd`
set zerofill=`keyval.sc< $DRVSBCW zerofill`
#
if ("$samplestoadd" != "") then
  # make sure we're not adding and skipping simultaneously
  if ("$numskip" != "") then
    echo "ERROR: Don't use both samplestoadd and numskip"
    exit(0)
  endif 
  # make a complex pair for each sample added
  dd if=/dev/zero bs=8 count=$samplestoadd  > tmp$$
else
  set samplestoadd = 0
  /bin/rm -f tmp$$
  touch tmp$$
endif
#
if ("$numskip" == "") then
  set numskip = 0
endif
#
if ("$nfftave" == "") then
  set nfftave = -1
endif
#
set ignore = ""
if ("$ignorefirst" == "yes" || "$ignorefirst" == "true" || "$ignorefirst" == "1") then
  set ignore = "-i"
endif
#
if ( "$firstrec" == "" ) then
  set firstrec=1
endif
if ( "$lastrec" == "" ) then
  set lastrec=2000000000
endif
#
if ("$zerofill" == "yes" || "$zerofill" == "true" || "$zerofill" == "1") then
  set nofillflag = ""
  set zerofill = "yes"
else
  set nofillflag = "-n"
  set zerofill = "no"
endif
#
@ torotate= $fftlen / 2 
#
set hdr=${basefile}_${fftlen}.hdr
if ( "$filesuffix" != "" ) then
  set hdr = ${hdr}.${filesuffix}
  if ( $#parms == 1 ) then
    set fsuf_is_one_number = `expr "$filesuffix" : '^0*\([0-9]\{1,\}\)$'`
    if ( "$fsuf_is_one_number" != "" ) then
      set fstart = $fsuf_is_one_number
      set fend = $fsuf_is_one_number
      echo "fsuf = ${filesuffix}: Processing scan ${fstart} only"
    else
      set fsuf_is_numerical_range = `expr "$filesuffix" : '^0*\([0-9]\{1,\}\)-[0-9]\{1,\}$'`
      if ( "$fsuf_is_numerical_range" != "" ) then
        set fstart = $fsuf_is_numerical_range
        set fend = `expr "$filesuffix" : '^[0-9]\{1,\}-0*\([0-9]\{1,\}\)$'`
        echo "fsuf = ${filesuffix}: Processing scans ${fstart} through ${fend}"
      else
        echo "fsuf = ${filesuffix}: Processing scan ${fstart} only"
      endif
    endif
  endif
endif
touch $hdr
echo "drv_cw.sc START       : `date`" >> $hdr
#     
#  
#
set tdir = `mktemp -d`
set fifoToUse=1
while ( $fifoToUse <= $numchan ) 
  set outfile=${basefile}_${fftlen}.p${fifoToUse}
  if ( "$filesuffix" != "" ) then
    set outfile = ${outfile}.${filesuffix}
  endif
  rm -f $outfile
  set scannum = $fstart
  while ( $scannum <= $fend )
    set ff=$tdir/data.$scannum.p${fifoToUse}
    switch ($fifoToUse)
    case 1:
    case 2:
      $stem/stripVme -h -o $scannum -n 1 -g "$firstrec $lastrec" < $inpfile |\
      $stem/unpriV -b $bits -i $numpol -f $fifoToUse -d 3 |\
      cat -s tmp$$ - |\
      $stem/convdatatype -q -i i4 -o f4 |\
      $stem/zerofill $nofillflag -b 8 -i $fftlen -o $fftlen |\
      $stem/fftfilter -d f -n $fftlen -s $numskip > $ff
      $stem/power_ao < $ff |\
      $stem/avgdata -d r4 -g $fftlen -h $nfftave $ignore |\
      $stem/rotate -i $fftlen -r $torotate >> $outfile
    breaksw
    case 3:
    case 4:
      @ which = $fifoToUse - 2
      set f1=$tdir/data.$scannum.p1
      set f2=$tdir/data.$scannum.p2
      $stem/stokes $f1 $f2  |\
      $stem/selectpnts -f $which -s 2 |\
      $stem/avgdata -d r4 -g $fftlen -h $nfftave $ignore |\
      $stem/rotate -i $fftlen -r $torotate >> $outfile
    endsw
      echo "Done with scan ${scannum}"
      @ scannum++
    
  end
  @ fifoToUse= $fifoToUse + 1
end
/bin/rm -r $tdir
#
rm tmp$$
#
echo "    length fft             : ${fftlen} " >> $hdr
echo "    numPol                 : ${numpol} " >> $hdr
echo "    data  file             : ${inpfile} " >> $hdr
echo "    file                   : $1 " >> $hdr
echo "    scans                  : $fstart-$fend" >> $hdr
echo "    records                : $firstrec-$lastrec" >> $hdr
echo "    nfftave,numskip        : $nfftave $numskip" >> $hdr
echo "    samplestoadd           : ${samplestoadd} " >> $hdr
echo "    zero-filled            : ${zerofill} " >> $hdr
echo "drv_cw.sc END         : `date`" >> $hdr
