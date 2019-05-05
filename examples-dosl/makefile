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
    INC_DOSL = -I$(DOSL_FOLDER)
else ifneq ("$(wildcard ../dosl)","")
    DOSL_FOLDER = ..
    INC_DOSL = -I../$(DOSL_FOLDER) # relative to the cpp sources
else
	DOSL_FOLDER = $(PREFIX)/include
	INC_DOSL = -I$(DOSL_FOLDER)
endif

INC_LOCAL = -I. -I..

# --------------------------------------------
# Macros

_BOLD = "\\033[1m"
_DIM = "\\033[2m"
_RED = "\\033[31m"
_GREEN = "\\033[32m"
_YELLOW = "\\033[33m"
_BLUE ="\\033[34m"
_GRAY = "\\033[37m"
_CLEAR = "\\033[0m"

# --------------------------------------------
# DOSL-specific

DOSL_ALGORITHM_LIST = $(shell find $(DOSL_FOLDER)/dosl/planners -type f -name '*.tcc' -exec basename {} .tcc \; | sort)

DOSL_AVAILABLE_ALGORITHM_STRING = "\nSome example programs require to choose an algorithm during compilation. Choose algorithm for those...\nNote: not all examples are compatible with all the algorithms.\n$(_YELLOW)$(_BOLD)Available algorithms:$(_CLEAR) $(_YELLOW)$(DOSL_ALGORITHM_LIST)$(_CLEAR) \n"

DOSL_DEFAULT_ALGORITHM = AStar

DOSL_ALGORITHM_X = $(shell \
                        if [ -z $(DOSL_ALGORITHM) ] ; then \
                            printf $(DOSL_AVAILABLE_ALGORITHM_STRING) >$$(tty); \
                            printf "$(_YELLOW)$(_BOLD)Enter algorithm name$(_CLEAR)$(_CLEAR) [$(_YELLOW)$(DOSL_DEFAULT_ALGORITHM)$(_CLEAR)]" >$$(tty); \
                            read -p ": " REPLY ; \
                            if [ -z $$REPLY ] ; then \
                                echo $(DOSL_DEFAULT_ALGORITHM) ; \
                            else \
                                echo $$REPLY ; \
                            fi \
                        else \
                            echo $(DOSL_ALGORITHM) ; \
                        fi )

DOSL_INTRO_STRING_DIR = "Including DOSL from '$(_YELLOW)$(DOSL_FOLDER)/$(_CLEAR)'.\n"
DOSL_INTRO_STRING_ALG = "$(eval DOSL_ALGORITHM ?= $(DOSL_ALGORITHM_X))Using algorithm '$(_YELLOW)$(DOSL_ALGORITHM)$(_CLEAR)' ($(_GRAY)$(_DIM)to change algorithm use 'make <target> DOSL_ALGORITHM=<algorithm>'$(_CLEAR)$(_CLEAR)).\n"


COMPILATION_SUCCESS_TEXT = "$(_GREEN)%s compilation successful!$(_CLEAR)\n"
COMPILATION_FAILURE_TEXT = "$(_RED)$(_BOLD)ERROR:$(_CLEAR) $(_RED)%s compilation failed!$(_CLEAR)\n"

COMPILE_PROGRAM = printf "$(_DIM)$(subst `,\`,$(strip $(2)))$(_CLEAR)\n" ; if $(2) ; then printf $(COMPILATION_SUCCESS_TEXT) "$(strip $(1))" ; else printf $(COMPILATION_FAILURE_TEXT) "$(strip $(1))" ; fi

# --------------------------------------------
# Dependencies


# --------------------------------------------
# common flags

CC = g++
CFLAGS = -std=gnu++11 -g -O3
PROFILE = -pg
WARNS = -w

LIBS = -lm
LIBS_OPENCV = `pkg-config --cflags --libs opencv` -L/usr/local/share/OpenCV/3rdparty/lib/  # The -L option fixes bug with libippicv
LIBS_OPENGL = -lglut -lGLU -lGL -lXi -lXmu -lX11 -lXext
LIBS_MATHEVAL = -lmatheval
#LIBS_ARMADILLO = -larmadillo  #`pkg-config --cflags --libs armadillo`

LIB_FOLDERS = -L/usr/local/lib

BIN_FOLDER = bin

# ================================
# --------------------------------

TEST_FOLDER=src/test
test: cv_test json_test matheval_test


#.PHONY: cv_test
cv_test:
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"
	@$(call COMPILE_PROGRAM,\
		$@ (test for libopencv-3.x),\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(TEST_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV)\
	)

#.PHONY: json_test
json_test:
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"
	@$(call COMPILE_PROGRAM,\
		$@ (test for RSJparser),\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(TEST_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS)\
	)

#.PHONY: matheval_test
matheval_test:
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"
	@$(call COMPILE_PROGRAM,\
		$@ (test for libmatheval-1.x),\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(TEST_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_MATHEVAL)\
	)


# --------------------------------------------

VPATH = bin


# --------------------------------------------
SIMPLE_FOLDER=src/simple
simple:  minimal2d_PathPlanning_AStar  minimal2d_PathPlanning_AnyAlg  map2d_encapsulated_PathPlanning  homotopy2d_encapsulated_PathPlanning


.PHONY: minimal2d_PathPlanning_AStar
minimal2d_PathPlanning_AStar:
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); 
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(PROFILE) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(SIMPLE_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS)\
	)

.PHONY: minimal2d_PathPlanning_AnyAlg
minimal2d_PathPlanning_AnyAlg:
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); 
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(PROFILE) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(SIMPLE_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS)\
	)

.PHONY: map2d_encapsulated_PathPlanning
map2d_encapsulated_PathPlanning: cv_test
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(PROFILE) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(SIMPLE_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)\
	)

.PHONY: homotopy2d_encapsulated_PathPlanning
homotopy2d_encapsulated_PathPlanning: cv_test
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(PROFILE) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(SIMPLE_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)\
	)


# --------------------------------------------
ADVANCED_FOLDER=src/advanced
advanced:  demo2d_PathPlanning  map2d_PathPlanning  map2d_VoronoiPartitioning  homotopy2d_PathPlanning  nonuniform2d_PathPlanning  robotArm_PathPlanning


.PHONY: demo2d_PathPlanning
demo2d_PathPlanning: cv_test
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)\
	)

.PHONY: map2d_PathPlanning
map2d_PathPlanning: cv_test json_test
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)\
	)

.PHONY: map2d_VoronoiPartitioning
map2d_VoronoiPartitioning: cv_test json_test
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)\
	)

.PHONY: homotopy2d_PathPlanning
homotopy2d_PathPlanning: cv_test json_test
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)\
	)

.PHONY: nonuniform2d_PathPlanning
nonuniform2d_PathPlanning: cv_test json_test matheval_test
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) $(LIBS_MATHEVAL) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)\
	)

.PHONY: robotArm_PathPlanning
robotArm_PathPlanning: cv_test json_test matheval_test
	@printf "\nNow compiling '$(_BLUE)$@$(_CLEAR)'...\n"; printf $(DOSL_INTRO_STRING_DIR); printf $(DOSL_INTRO_STRING_ALG);
	@$(call COMPILE_PROGRAM,\
		$@,\
		$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $(BIN_FOLDER)/$@ $(ADVANCED_FOLDER)/$@.cpp $(LIB_FOLDERS) $(LIBS) $(LIBS_OPENCV) $(LIBS_MATHEVAL) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)\
	)


# --------------------------------------------

cv_examples:  demo2d_PathPlanning  map2d_PathPlanning  homotopy2d_PathPlanning  nonuniform2d_PathPlanning  robotArm_PathPlanning  map2d_encapsulated_PathPlanning  homotopy2d_encapsulated_PathPlanning

encapsulated_examples: map2d_encapsulated_PathPlanning  homotopy2d_encapsulated_PathPlanning

all: test simple advanced

# --------------------------------------------

ALL_SAME_ALG_INTRO = "Using algorithm '$(_YELLOW)$(DOSL_ALGORITHM)$(_CLEAR)' for all ($(_GRAY)$(_DIM)Available algorithms: $(DOSL_ALGORITHM_LIST)$(_CLEAR)$(_CLEAR)). \n"

set_AStar:
	@printf "$(eval DOSL_ALGORITHM = AStar)\n"; printf $(ALL_SAME_ALG_INTRO);
all_AStar: set_AStar all


set_ThetaStar:
	@printf "$(eval DOSL_ALGORITHM = ThetaStar)\n"; printf $(ALL_SAME_ALG_INTRO);
all_ThetaStar: set_ThetaStar all


set_SStar:
	@printf "$(eval DOSL_ALGORITHM = SStar)\n"; printf $(ALL_SAME_ALG_INTRO);
all_SStar: set_SStar all

# --------------------------------------------

run:
	@for filename in bin/* ; do \
         printf "\n\nNow running $(_BLUE)./$$filename$(_CLEAR) ...\n" ; \
         printf "$(_GRAY)$(_DIM)-----------$(_CLEAR)$(_CLEAR)\n\n" ; \
         ./$$filename ; \
    done


clean:
	rm bin/*

