# build test executables
ADD_EXECUTABLE(test/moment moment.cpp $<TARGET_OBJECTS:utils>)
ADD_TEST(build_tests "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target test/moment)

# execute test executables
ADD_TEST(moment ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test/moment)
SET_TESTS_PROPERTIES(moment PROPERTIES DEPENDS build_tests)

# build main executables
ADD_TEST(build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target all)

# noise generation
ADD_TEST(noise_ref ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/analyze
  scale=-.5 fmin=20 fmax=20000 ln=2 fftlen=65536 mst wspec=noise.spec wref=noise.ref out=noise.wav)
SET_TESTS_PROPERTIES(noise_ref PROPERTIES DEPENDS build)

ADD_TEST(noise_ref_spec perl ${CMAKE_CURRENT_SOURCE_DIR}/datacheck.pl
  noise.spec 2,6 "{ $_[1] < 20 || $_[1] > 20000 ? (0, undef) : (($_[1]/20.5078) ** -.5, $_[0]&1 ? -1 : 1) }" .0001 )
SET_TESTS_PROPERTIES(noise_ref_spec PROPERTIES DEPENDS noise_ref)

# analyze noise
ADD_TEST(analyze_noise ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/analyze
  in=noise.wav fftlen=65536 mfft psa=22 diff wd=analyze.dat)
SET_TESTS_PROPERTIES(analyze_noise PROPERTIES DEPENDS noise_ref)

ADD_TEST(analyze_noise_dat perl ${CMAKE_CURRENT_SOURCE_DIR}/datacheck.pl
  analyze.dat 2,4,6,7,8,9 "{ $_[1] < 20 || $_[1] > 20000 ? (0, 0) : (($_[0]&1) * 3.67722 * ($_[1]/20.5078) ** -.5, 3.67722 * ($_[1]/20.5078) ** -.5, $_[0]&1, ($_[0]&1 ? -180 : undef), -($_[0]&1), 0) }" .0005 .0001 )
SET_TESTS_PROPERTIES(analyze_noise_dat PROPERTIES DEPENDS analyze_noise)

# generate sweep
ADD_TEST(sweep_ref ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/analyze
  fmin=20 fmax=20000 ln=2 flog=^/12 fftlen=65536 msweep wspec=sweep.spec wref=sweep.ref out=sweep.wav)
SET_TESTS_PROPERTIES(sweep_ref PROPERTIES DEPENDS build)

ADD_TEST(sweep_ref_spec perl ${CMAKE_CURRENT_SOURCE_DIR}/datacheck.pl
  sweep.spec 2,6 "{
    state $nf;
    return (0,undef) if $_[1] < 20 || $_[1] > 20000 || $_[0] < $nf;
    $nf = $_[0] * 2**(1/12) + .5;
    return (1,1);
    }" .0001 .0001 )
#  sweep.spec 2,6 "{ $_[1] < 20 || $_[1] > 20000 || abs(2 ** int(log($_[1]/20.5078)/log(2)*12+.5)/12 * 20.5078 - $_[1]) > .01 ? (0, undef) : (1, 1) }" .0001 .0001 )
SET_TESTS_PROPERTIES(sweep_ref_spec PROPERTIES DEPENDS sweep_ref)

