/**/
CALL SysSetPriority 1,1
fname = STRIP(ARG(1), 'B', '"')
CALL STREAM fname, 'c', 'open read'
type = CHARIN(fname, , 4)
CALL STREAM fname, 'c', 'close'
SELECT
 WHEN type = 'OggS' THEN
   'oggdec "'fname'" -o - | ..\bin\analyzeE psa22 fq44100 scm1 win5 n262144 loop mfft ud'
 WHEN type = 'RIFF' THEN
   '..\bin\analyzeE psa22 fq44100 scm1 win5 n262144 loop mfft ud "in'fname'"'
 OTHERWISE
   'lame --mp3input --decode -t "'fname'" - | ..\bin\analyzeE fq44100 scm1 win5 n262144 loop mfft ud'
   END

