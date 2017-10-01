# Discrete Optimal search Library (DOSL)
# Makefile for compiling examples

SHELL := /bin/bash

# Change these if needed or pass as input parameter
# PREFIX = /usr/local

# ================================

# INCLUDE FOLDER
# --------------

ifneq ("$(wildcard /usr/local/include/dosl)","")
    DOSL_FOLDER = /usr/local/include
else ifneq ("$(wildcard ../dosl)","")
    DOSL_FOLDER = ..
else
	DOSL_FOLDER = $(PREFIX)/include
endif

INC_DOSL = -I../$(DOSL_FOLDER) # relative to the cpp sources
INC_LOCAL = -I. -I..

# --------------------------------------------
# DOSL-specific

DOSL_ALGORITHM_LIST = $(shell find $(DOSL_FOLDER)/dosl/planners -type f -name '*.tcc' -exec basename {} .tcc \;)

DOSL_AVAILABLE_ALGORITHM_STRING = "\nSome example programs require to choose an algorithm during compilation. Choose algorithm for those...\nNote: not all examples are compatible with all the algorithms.\n\033[33m\033[1mAvailable algorithms:\033[21m $(DOSL_ALGORITHM_LIST)\033[39m \n"

DOSL_DEFAULT_ALGORITHM = AStar

DOSL_ALGORITHM_X = $(shell \
                        if [ -z $(DOSL_ALGORITHM) ] ; then \
                            printf $(DOSL_AVAILABLE_ALGORITHM_STRING) >$$(tty); \
                            printf "\033[33m\033[1mEnter algorithm name\033[21m\033[39m [\033[33m$(DOSL_DEFAULT_ALGORITHM)\033[39m]" >$$(tty); \
                            read -p ": " REPLY ; \
                            if [ -z $$REPLY ] ; then \
                                echo $(DOSL_DEFAULT_ALGORITHM) ; \
                            else \
                                echo $$REPLY ; \
                            fi \
                        else \
                            echo $(DOSL_ALGORITHM) ; \
                        fi )

DOSL_INTRO_STRING_DIR = "Including DOSL from '\033[33m$(DOSL_FOLDER)/\033[39m'.\n"
DOSL_INTRO_STRING_ALG = "$(eval DOSL_ALGORITHM ?= $(DOSL_ALGORITHM_X))Using algorithm '\033[33m$(DOSL_ALGORITHM)\033[39m' (\033[37m\033[2mto change algorithm use 'make <target> DOSL_ALGORITHM=<algorithm>'\033[21m\033[39m).\n"

TEST_SUCCESS_TEXT = "\033[32mTESTNAME compilation test successful!\033[39m\n"
#TEST_FAILURE_TEXT = "\033[31m\033[1mERROR:\033[21m TESTNAME compilation test failed! Please make sure that LIBNAME is installed and is in the compiler's search path.\033[39m\n"
TEST_FAILURE_TEXT = "\033[31m\033[1mERROR:\033[21m %s compilation test failed! Please make sure that %s is installed and is in the compiler's search path.\033[39m\n"

# --------------------------------------------
# common flags

CC = g++
CFLAGS = -std=gnu++11 -g -O3
PROFILE = -pg
WARNS = -w

LIBS = -lm
LIBS_OPENCV = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs
LIBS_OPENGL = -lglut -lGLU -lGL -lXi -lXmu -lX11 -lXext
LIBS_ARMADILLO = -larmadillo
LIBS_MATHEVAL = -lmatheval

LIB_FOLDERS = -L/usr/local/lib

BIN_FOLDER = bin

# --------------------------------------------
TEST_FOLDER=src/test
test: cv_test json_test matheval_test


#.PHONY: cv_test
cv_test:
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"
	@if $(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(TEST_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) ; \
	then printf ${TEST_SUCCESS_TEXT:TESTNAME=OpenCV} ; else printf $(TEST_FAILURE_TEXT) "OpenCV" "libopencv-3.x" ; false ; fi

#.PHONY: json_test
json_test:
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"
	@if $(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(TEST_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) ; \
	then printf ${TEST_SUCCESS_TEXT:TESTNAME=RSJp} ; else printf $(TEST_FAILURE_TEXT) "RSJp" "RSJparser" ; false ; fi

#.PHONY: matheval_test
matheval_test:
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"
	@if $(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(TEST_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_MATHEVAL) ; \
	then printf ${TEST_SUCCESS_TEXT:TESTNAME=MathEval} ; else printf $(TEST_FAILURE_TEXT) "MathEval" "libmatheval-1.x" ; false ; fi


# --------------------------------------------

VPATH = bin


# --------------------------------------------
SIMPLE_FOLDER=src/simple
simple:  minimal2d_PathPlanning_AStar  minimal2d_PathPlanning_AnyAlg  map2d_encapsulated_PathPlanning  homotopy2d_encapsulated_PathPlanning


.PHONY: minimal2d_PathPlanning_AStar
minimal2d_PathPlanning_AStar:
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); 
	$(CC) $(CFLAGS) $(PROFILE) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(SIMPLE_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS)

.PHONY: minimal2d_PathPlanning_AnyAlg
minimal2d_PathPlanning_AnyAlg:
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); 
	$(CC) $(CFLAGS) $(PROFILE) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(SIMPLE_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS)

.PHONY: map2d_encapsulated_PathPlanning
map2d_encapsulated_PathPlanning: cv_test
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	$(CC) $(CFLAGS) $(PROFILE) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(SIMPLE_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)

.PHONY: homotopy2d_encapsulated_PathPlanning
homotopy2d_encapsulated_PathPlanning: cv_test
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	$(CC) $(CFLAGS) $(PROFILE) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(SIMPLE_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)


# --------------------------------------------
ADVANCED_FOLDER=src/advanced
advanced:  demo2d_PathPlanning  map2d_PathPlanning  homotopy2d_PathPlanning  nonuniform2d_PathPlanning  robotArm_PathPlanning


.PHONY: demo2d_PathPlanning
demo2d_PathPlanning: cv_test
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)

.PHONY: map2d_PathPlanning
map2d_PathPlanning: cv_test json_test matheval_test
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)

.PHONY: homotopy2d_PathPlanning
homotopy2d_PathPlanning: cv_test json_test
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)

.PHONY: nonuniform2d_PathPlanning
nonuniform2d_PathPlanning: cv_test json_test matheval_test
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) $(LIBS_MATHEVAL) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)

.PHONY: robotArm_PathPlanning
robotArm_PathPlanning: cv_test json_test matheval_test
	@printf "\nNow compiling '\033[34m$@\033[39m'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) $(LIBS_MATHEVAL) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)


# --------------------------------------------

cv_examples:  demo2d_PathPlanning  map2d_PathPlanning  homotopy2d_PathPlanning  nonuniform2d_PathPlanning  robotArm_PathPlanning  map2d_encapsulated_PathPlanning  homotopy2d_encapsulated_PathPlanning

encapsulated_examples: map2d_encapsulated_PathPlanning  homotopy2d_encapsulated_PathPlanning

all: test simple advanced

# --------------------------------------------

run:
	@for filename in bin/* ; do \
         printf "\n\nNow running \033[34m./$$filename\033[39m ...\n" ; \
         printf "\033[37m\033[2m-----------\033[21m\033[39m\n\n" ; \
         ./$$filename ; \
    done


clean:
	rm bin/*

