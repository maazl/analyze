/**/

PARSE ARG cfgfile opt

IF cfgfile = '' THEN
   CALL Error 40, 'usage: measure <config> [options]'

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
   CALL CfgQ 'scale', 0
   CALL CfgQ 'rref'
   CALL CfgQ 'fbin', .05
   CALL CfgQ 'plotcmd', "l 'viewfft'"
   CALL CfgQ 'famin', cfg.fmin
   CALL CfgQ 'famax', cfg.fmax
   CALL CfgQ 'zerofile', 'zero.dat'

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
   'start /MIN /C sp f1 o ref.exe 'cfg.fftlen cfg.fmin cfg.fmax cfg.scale' -1 ^| sp f2 o buffer2 /p=40k /b=16M - \pipe\refplay.wav'
   CALL SysSleep 1
   'start /C sp f3 o playrec \pipe\refplay.wav /bufcnt:8 /i:'cfg.odevice
   /*'CALL play FILE="\pipe\refplay.wav"'*/
   /*CALL SysSetPriority 2, -1*/
   CALL SysSleep 1

   IF cfg.initexec \= '' THEN
      cfg.initexec

   "playrec /f:"cfg.fsamp" /i:"cfg.idevice" /v:100 /r con | sp f3 o buffer2 /b=16M - - | analyze ""zf"cfg.zerofile""" psa32768 loop fq"cfg.fsamp" rref"cfg.rref" scm1 mfft n"cfg.fftlen" wd ""plot"cfg.plotcmd""" fmax"cfg.fmax" fbin"cfg.fbin" fmin"cfg.fmin" famin"cfg.famin" famax"cfg.famax" |gnuplot gpenv -"
   END

 WHEN cfg.mtype = 'sweep' THEN DO
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
   'start /MIN /C ref.exe 'cfg.fftlen cfg.fmin cfg.fmax' 0 30 ^| sp f2 o buffer2 /p=40k /b=16M - \pipe\refplay.wav'
   CALL SysSleep 1
   'start /C playrec \pipe\refplay.wav /bufcnt:8 /i:'cfg.odevice
   /*'CALL play FILE="\pipe\refplay.wav"'*/
   /*CALL SysSetPriority 2, -1*/
   CALL SysSleep 1

   IF cfg.initexec \= '' THEN
      'CALL 'cfg.initexec

   tmax = (65536 + 20*cfg.fftlen) / cfg.fsamp + 1.0;
   "playrec /f:"cfg.fsamp" /i:"cfg.idevice" /v:100 /r /ct:"tmax" con 2>log | sp f3 o buffer2 /b=16M - - | analyze pte psa65536 n"cfg.fftlen" scm1 mfft fq"cfg.fsamp" gg famin"cfg.famin" famax"cfg.famax" ln20 wd"
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
   'start /MIN /C sp f1 o ref.exe 'cfg.fftlen cfg.fmin cfg.fmax' 0 'cyc+10' ^| sp f2 o buffer2 /p=40k /b=16M - \pipe\refplay.wav'
   CALL SysSleep 1
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
   "playrec /f:"cfg.fsamp" /i:"cfg.idevice" /v:100 /r con 2>log | sp f3 o buffer2 /b=16M - - | analyze pte psa65536 n"cfg.fftlen" scm1 mfft fq"cfg.fsamp" "zmode" zn ""zf"cfg.zerofile""" famin"cfg.famin" famax"cfg.famax" ln"cfg.loops" wd"
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
         IF val = '' THEN DO
            SAY 'Invalid line in configuration file'
            EXIT 10
            END
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


