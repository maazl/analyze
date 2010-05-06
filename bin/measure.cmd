/**/

PARSE ARG cfgfile opt

IF cfgfile = '' THEN
   CALL Error 40, 'usage: measure <config> [options]'

/* setup path */
CALL VALUE 'PATH', 'bin;'VALUE('PATH',,'OS2ENVIRONMENT'), 'OS2ENVIRONMENT'

/* read configs */
CALL ReadConfig 'setup.cfg'
CALL ReadConfig cfgfile'.cfg'

/* check required parameters and apply defaults */
/* setup */
CALL CfgQ 'device', 0
CALL CfgQ 'idevice', cfg.device
CALL CfgQ 'odevice', cfg.device
CALL CfgQ 'preexec', ''
CALL CfgQ 'initexec', ''
CALL CfgQ 'postexec', ''
CALL CfgQ 'xopt', ''

IF cfg.preexec \= '' THEN
   cfg.preexec

/* measurement type */
CALL CfgQ 'mtype'

CALL CfgQ 'fsamp', 48000

SELECT
 WHEN cfg.mtype = 'fft' THEN DO
   CALL CfgQ 'fftlen'
   CALL CfgQ 'fmin'
   CALL CfgQ 'fmax'
   CALL CfgQ 'scm', 1
   CALL CfgQ 'scale', 0
   CALL CfgQ 'rref', 1
   CALL CfgQ 'fbin', .05
   CALL CfgQ 'plotcmd', "l 'viewfft'"
   CALL CfgQ 'famin', cfg.fmin
   CALL CfgQ 'famax', cfg.fmax
   CALL CfgQ 'zerofile', ''
   CALL CfgQ 'harm', 0
   CALL CfgQ 'mst', 0
   CALL CfgQ 'finc', 1
   CALL CfgQ 'flog', 0 
   
   opt = ''
   IF cfg.zerofile \= '' THEN
      opt = opt" zr ""zf"cfg.zerofile""""

   /* generate reference file
   'ref.exe 'cfg.fftlen cfg.fmin cfg.fmax cfg.scale' 10 ref.wav'
   IF RC \= 0 THEN
      CALL Error RC, 'ref.exe failed.' */

   /* gnuplot environment */
   CALL STREAM 'gpenv', 'c', 'open write replace'
   CALL LINEOUT 'gpenv', 'fmin='cfg.fmin
   CALL LINEOUT 'gpenv', 'fmax='cfg.fmax
   CALL LINEOUT 'gpenv', 'rref='cfg.rref
   CALL LINEOUT 'gpenv', 'fftlen='cfg.fftlen
   CALL LINEOUT 'gpenv', 'fsamp='cfg.fsamp
   CALL STREAM 'gpenv', 'c', 'close'

   /*CALL SysSetPriority 4, 1*/
   '@pstat|grep -i NOISE.EXE -q'
   IF RC \= 0 THEN DO
      'start /MIN /C sp f3 o noise.exe bn'cfg.fftlen' fmin'cfg.fmin' fmax'cfg.fmax' harm'cfg.harm' finc'cfg.finc' flog'cfg.flog' scale'cfg.scale' mst'cfg.mst' loop wdnoise.dat wrref.dat ww- ^| sp t1 o buffer2 -p=60k -b=8M -h=100%% - \pipe\refplay.wav'
      CALL SysSleep 1
      END
   'start /C sp t2 o playrec \pipe\refplay.wav /bufcnt:8 /i:'cfg.odevice
   /*'CALL play FILE="\pipe\refplay.wav"'*/
   /*CALL SysSetPriority 2, -1*/
   CALL SysSleep 2

   IF cfg.initexec \= '' THEN
      cfg.initexec

   "playrec /f:"cfg.fsamp" /i:"cfg.idevice" /v:100 /r con | sp t1 o buffer2 -b=32M - - | analyze psa32000 loop fq"cfg.fsamp" rref"cfg.rref" scm"cfg.scm" mfft he bn"cfg.fftlen" mst"cfg.mst" wd ""plot"cfg.plotcmd""" fmax"cfg.fmax" fbin"cfg.fbin" fmin"cfg.fmin" harm"cfg.harm" finc"cfg.finc" flog"cfg.flog" famin"cfg.famin" famax"cfg.famax" "opt" "cfg.xopt"|gnuplot gpenv -"
   END

 WHEN cfg.mtype = 'sweep' THEN DO
   CALL CfgQ 'fftlen'
   CALL CfgQ 'fmin'
   CALL CfgQ 'fmax'
   /*CALL CfgQ 'scm', 1*/
   CALL CfgQ 'flog', 1.05946309
   /*CALL CfgQ 'plotcmd', "l 'viewfft'"*/
   CALL CfgQ 'sync'
   CALL CfgQ 'synclevel'
   CALL CfgQ 'syncphase'
   CALL CfgQ 'overlap', 0
   CALL CfgQ 'zerofile', ''

   opt = ''
   IF cfg.zerofile \= '' THEN
      opt = opt" zr ""zf"cfg.zerofile""""

   /* gnuplot environment */
   /*CALL STREAM 'gpenv', 'c', 'open write replace'
   CALL LINEOUT 'gpenv', 'fmin='cfg.fmin
   CALL LINEOUT 'gpenv', 'fmax='cfg.fmax
   CALL LINEOUT 'gpenv', 'rref='cfg.rref
   CALL LINEOUT 'gpenv', 'fftlen='cfg.fftlen
   CALL LINEOUT 'gpenv', 'fsamp='cfg.fsamp
   CALL STREAM 'gpenv', 'c', 'close'*/
   opt = opt||'n'cfg.fftlen' fmin'cfg.fmin' fmax'cfg.fmax' sync'cfg.sync' slvl'cfg.synclevel' sph'cfg.syncphase' sov'cfg.overlap

   /*CALL SysSetPriority 4, 1*/
   '@pstat|grep -i sweep\.exe -q'
   IF RC \= 0 THEN DO
      'start /MIN /C sp f1 o sweep.exe mr psa60000 'opt' ^| sp f2 o buffer2 -p=60k -b=8M - \pipe\refplay.wav'
      CALL SysSleep 1
      END
   'start /C sp f3 o playrec \pipe\refplay.wav /bufcnt:8 /i:'cfg.odevice
   /*CALL SysSetPriority 2, -1*/
   CALL SysSleep 1

   IF cfg.initexec \= '' THEN
      cfg.initexec

   "playrec /f:"cfg.fsamp" /i:"cfg.idevice" /v:100 /r con | sp f3 o buffer2 -b=32M - - | sweep psa10000 ma "opt" >log"
   END

/* WHEN cfg.mtype = 'hyst' THEN DO

   END*/

 WHEN cfg.mtype = 'gain' THEN DO
   CALL CfgQ 'fftlen'
   CALL CfgQ 'fmin', 5
   CALL CfgQ 'fmax', cfg.fsamp/2
   CALL CfgQ 'famin', cfg.fmin
   CALL CfgQ 'famax', cfg.fmax

   CALL SysSetPriority 4, 1
   '@pstat|grep -i noise.exe -q'
   IF RC \= 0 THEN DO
      'start /MIN /C noise.exe bn'cfg.fftlen' fmin'cfg.fmin' fmax'cfg.fmax' ln30 ww- ^| sp f2 o buffer2 -p=40k -b=16M - \pipe\refplay.wav'
      CALL SysSleep 1
      END
   'start /C playrec \pipe\refplay.wav /bufcnt:8 /i:'cfg.odevice
   /*'CALL play FILE="\pipe\refplay.wav"'*/
   /*CALL SysSetPriority 2, -1*/
   CALL SysSleep 1

   IF cfg.initexec \= '' THEN
      'CALL 'cfg.initexec

   tmax = (65536 + 20*cfg.fftlen) / cfg.fsamp + 1.0;
   "playrec /f:"cfg.fsamp" /i:"cfg.idevice" /v:100 /r /ct:"tmax" con 2>log | sp f3 o buffer2 -b=16M - - | analyze pte psa65536 bn"cfg.fftlen" scm1 mfft fq"cfg.fsamp" gg famin"cfg.famin" famax"cfg.famax" ln20 wd"
   END

 WHEN cfg.mtype = 'calibrate' THEN DO
   CALL CfgQ 'fftlen'
   CALL CfgQ 'fmin', 5
   CALL CfgQ 'fmax', cfg.fsamp/2
   CALL CfgQ 'famin', cfg.fmin
   CALL CfgQ 'famax', cfg.fmax
   CALL CfgQ 'loops', 20
   CALL CfgQ 'zerofile', 'zero.dat'

   /* gnuplot environment */
   CALL STREAM 'gpenv', 'c', 'open write replace'
   CALL LINEOUT 'gpenv', 'fmin='cfg.fmin
   CALL LINEOUT 'gpenv', 'fmax='cfg.fmax
   CALL LINEOUT 'gpenv', 'fftlen='cfg.fftlen
   CALL LINEOUT 'gpenv', 'fsamp='cfg.fsamp
   CALL STREAM 'gpenv', 'c', 'close'

   cyc = cfg.loops*2+10
   /*CALL SysSetPriority 4, 1*/
   '@pstat|grep -i noise.exe -q'
   IF RC \= 0 THEN DO
      'start /MIN /C sp f1 o noise.exe bn'cfg.fftlen' fmin'cfg.fmin' fmax'cfg.fmax' ln'cyc+10' ww- ^| sp f2 o buffer2 -p=40k -b=8M - \pipe\refplay.wav'
      CALL SysSleep 1
      END
   'start /C sp f3 o playrec \pipe\refplay.wav /bufcnt:8 /i:'cfg.odevice
   /*'CALL play FILE="\pipe\refplay.wav"'*/
   /*CALL SysSetPriority 2, -1*/
   CALL SysSleep 1

   IF cfg.initexec \= '' THEN
      cfg.initexec

   tmax = (65536 + cyc*cfg.fftlen) / cfg.fsamp + 1.0;
   IF opt = 'verify' THEN
      zmode = 'zd'
    ELSE
      zmode = 'zg'
   SAY
   SAY "***** Z = +Inf first, Z = 0 on request."
   SAY
   "playrec /f:"cfg.fsamp" /i:"cfg.idevice" /v:100 /r con 2>log | sp f3 o buffer2 -b=16M - - | analyze pte psa65536 bn"cfg.fftlen" scm1 mfft fq"cfg.fsamp" "zmode" zn ""zf"cfg.zerofile""" famin"cfg.famin" famax"cfg.famax" ln"cfg.loops" wd"
   END

 WHEN cfg.mtype = 'hyst' THEN DO
   CALL CfgQ 'fftlen'
   CALL CfgQ 'rref'
   CALL CfgQ 'fbase'
   CALL CfgQ 'plotcmd', "l 'viewhyst'"
   CALL CfgQ 'zerofile', 'zero.dat'
   CALL CfgQ 'shape', 'triangle'

   /* gnuplot environment */
   CALL STREAM 'gpenv', 'c', 'open write replace'
   CALL LINEOUT 'gpenv', 'rref='cfg.rref
   CALL LINEOUT 'gpenv', 'fftlen='cfg.fftlen
   CALL LINEOUT 'gpenv', 'fsamp='cfg.fsamp
   CALL STREAM 'gpenv', 'c', 'close'

   /* calculate harmonic closest to base frequency */
   harmonic = TRUNC(cfg.fsamp / cfg.fbase +.5)

   /*CALL SysSetPriority 4, 1*/
   '@pstat|grep -i REF.EXE -q'
   IF RC \= 0 THEN DO
      'start /MIN /C sp f1 o ref.exe 'cfg.fftlen'/'harmonic cfg.fsamp cfg.shape' ^| sp f2 o buffer2 -p=40k -b=8M - \pipe\refplay.wav'
      CALL SysSleep 1
      END
   'start /C sp f3 o playrec \pipe\refplay.wav /bufcnt:8 /i:'cfg.odevice
   /*'CALL play FILE="\pipe\refplay.wav"'*/
   /*CALL SysSetPriority 2, -1*/
   CALL SysSleep 1

   IF cfg.initexec \= '' THEN
      cfg.initexec

   "playrec /f:"cfg.fsamp" /i:"cfg.idevice" /v:100 /r con | sp f3 o buffer2 -b=16M - - | analyze ""zf"cfg.zerofile""" zr pdc2 psa32768 loop fq"cfg.fsamp" rref"cfg.rref" scm1 mxy bn"cfg.fftlen" har"2" wd"
   END

 OTHERWISE
   CALL Error 40, 'Unknown measurement method 'cfg.mtype
   END

IF cfg.postexec \= '' THEN
   'CALL 'cfg.postexec

EXIT

/* Read configuration file
   ARG(1)  filename
*/
ReadConfig: PROCEDURE EXPOSE cfg.
   IF STREAM(ARG(1), 'c', 'open read') \= 'READY:' THEN DO
      SAY 'Failed to open configuration file 'ARG(1)
      EXIT 10
      END
   DO WHILE STREAM(ARG(1), 's') = 'READY'
      l = STRIP(LINEIN(ARG(1)))
      SELECT
       WHEN LENGTH(l) = 0 | LEFT(l,1) = '#' THEN /* comment */
         NOP
       WHEN LEFT(l,1) = '@' THEN /* include */
         CALL ReadConfig SUBSTR(l,2)
       OTHERWISE
         PARSE VAR l key'='val
         CALL VALUE 'cfg.'key, val
         END
      END
   CALL STREAM ARG(1), 'c', 'close'
   RETURN

/* fetch config param */
CfgQ:
   IF VAR('cfg.'ARG(1)) THEN
      RETURN VALUE('cfg.'ARG(1))
   IF ARG(2,'e') THEN DO
      CALL VALUE 'cfg.'ARG(1), ARG(2)
      RETURN ARG(2)
      END
   CALL Error 30, 'Configuration parameter 'ARG(1)' is missing.'

Esc: PROCEDURE
   RETURN CHANGESTR("'", ARG(1), "^'")

Error: PROCEDURE
   CALL LINEOUT STDERR, ARG(2)
   EXIT ARG(1)


