INC_DOSL = ./dosl-3.0
TEST_DIR = ./test_dosl-3.0

# ================================

# TESTS
# -----

CC = g++
CFLAGS = -std=gnu++11 -O3 -g 
WARNS = -w

LIBS = -lm
LIBS_OPENCV = -lopencv_core -lopencv_highgui -lopencv_imgproc
LIBS_OPENGL = -lglut -lGLU -lGL -lXi -lXmu -lX11 -lXext
LIB_ARMADILLO = -larmadillo

INTRO_STRING = "\nNOTE: Including DOSL from $(INC_DOSL). \n      If you have installed DOSL and want to use the installation, this should be empty (change 'INC_DOSL' in makefile).\n\n"
OUTTRO_STRING = "\nTo run an executable: $(TEST_DIR)/<executable_name>\n\n"

tests: minimal_benchmark

.PHONY: minimal_benchmark
minimal_benchmark:
	@printf $(INTRO_STRING)
	$(CC) $(CFLAGS) $(WARNS) -I. -I$(INC_DOSL) -o $(TEST_DIR)/minimal_benchmark.o -c $(TEST_DIR)/minimal_benchmark.cpp
	$(CC) $(CFLAGS) $(WARNS) -I. -I$(INC_DOSL) -o $(TEST_DIR)/minimal_benchmark $(TEST_DIR)/minimal_benchmark.o $(LIBS)
	@printf $(OUTTRO_STRING)

clean:
	rm $(prog)
	rm $(prog).o

