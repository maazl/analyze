/**/
fname = STRIP(ARG(1), 'B', '"')
CALL STREAM fname, 'c', 'open read'
type = CHARIN(fname, , 4)
CALL STREAM fname, 'c', 'close'
SELECT
 WHEN type = 'OggS' THEN
   'oggdec "'fname'" -o - | analyzeE psa22 fq44100 scm1 win5 n262144 loop mfft ud'
 WHEN type = 'RIFF' THEN
   'analyzeE psa22 fq44100 scm1 win5 n262144 loop mfft ud "in'fname'"'
 OTHERWISE
   'lame --mp3input --decode -t "'fname'" - | analyzeE fq44100 scm1 win5 n262144 loop mfft ud'
   END

