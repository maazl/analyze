/**/
CALL STREAM ARG(1), 'c', 'open read'
type = CHARIN(ARG(1), , 4)
CALL STREAM ARG(1), 'c', 'close'
IF type = 'OggS' THEN
   'rem oggdec 'ARG(1)' -o - | d:analyzeE psa22 fq44100 scm1 win5 n65536 loop mfft ud'
 ELSE
   'lame --mp3input --decode -t 'ARG(1)' - | d:analyzeE fq44100 scm1 win5 n65536 loop mfft ud'

