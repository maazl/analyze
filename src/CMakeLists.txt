FIND_PACKAGE(Threads REQUIRED)

# Work around for the insufficiency of cmake to handle common dependencies of targets w/o building them repeatedly.
# cmake will build each source file as often as it is referenced.
ADD_LIBRARY(utils OBJECT
	utils.cpp
	parser.cpp
	nrc.cpp
	moment.cpp
	interpolation.cpp
	pcminput.cpp
	pcmoutput.cpp
)

ADD_EXECUTABLE(analyze analyze.cpp analyzein.cpp analyzeout.cpp $<TARGET_OBJECTS:utils>)
TARGET_LINK_LIBRARIES(analyze fftw3f ${CMAKE_THREAD_LIBS_INIT})
ADD_EXECUTABLE(noise noise.cpp $<TARGET_OBJECTS:utils>)
TARGET_LINK_LIBRARIES(noise fftw3f)
ADD_EXECUTABLE(sweep sweep.cpp $<TARGET_OBJECTS:utils>)
TARGET_LINK_LIBRARIES(sweep fftw3f ${CMAKE_THREAD_LIBS_INIT})
ADD_EXECUTABLE(ref ref.cpp $<TARGET_OBJECTS:utils>)
ADD_EXECUTABLE(sampledist sampledist.cpp $<TARGET_OBJECTS:utils>)