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

INC_DOSL = -I$(DOSL_FOLDER)
INC_LOCAL = -I. -Iinclude

# --------------------------------------------
# DOSL-specific

DOSL_ALGORITHM_LIST = $(shell find $(DOSL_FOLDER)/dosl/planners -type f -name '*.tcc' -exec basename {} .tcc \;) # AStar SStar
DOSL_ALGORITHM = AStar

DOSL_INTRO_STRING = "Including DOSL from '\033[33m$(DOSL_FOLDER)/\033[39m'\nUsing algorithm '\033[33m$(DOSL_ALGORITHM)\033[39m' [\033[33mAvailable algorithms: $(DOSL_ALGORITHM_LIST); To change algorithm use 'make <target> DOSL_ALGORITHM=<algorithm>'\033[39m].\n"

# --------------------------------------------
# common flags

CC = -g++
CFLAGS = -std=gnu++11 -g -Og 
WARNS = -w

LIBS = -lm
LIBS_OPENCV = -lopencv_core -lopencv_highgui -lopencv_imgproc
LIBS_OPENGL = -lglut -lGLU -lGL -lXi -lXmu -lX11 -lXext
LIB_ARMADILLO = -larmadillo

# --------------------------------------------

basic_examples: minimal2d_PathPlanning_AStar simple2d_PathPlanning
cv_examples:    demo2d_PathPlanning  map2d_PathPlanning  homotopy2d_PathPlanning

all_examples: basic_examples cv_examples

.PHONY: minimal2d_PathPlanning_AStar
minimal2d_PathPlanning_AStar:
	@printf "\nNow compiling '$@'...\n"
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $@ $@.cpp $(LIBS)

.PHONY: simple2d_PathPlanning
simple2d_PathPlanning:
	@printf "\nNow compiling '$@'...\n"
	@printf $(DOSL_INTRO_STRING)
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $@ $@.cpp $(LIBS) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)

.PHONY: demo2d_PathPlanning
demo2d_PathPlanning:
	@printf "\nNow compiling '$@'...\n"
	@printf $(DOSL_INTRO_STRING)
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $@ $@.cpp $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)

.PHONY: map2d_PathPlanning
map2d_PathPlanning:
	@printf "\nNow compiling '$@'...\n"
	@printf $(DOSL_INTRO_STRING)
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $@ $@.cpp $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)

.PHONY: homotopy2d_PathPlanning
homotopy2d_PathPlanning:
	@printf "\nNow compiling '$@'...\n"
	@printf $(DOSL_INTRO_STRING)
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) $(INC_DOSL) -o $@ $@.cpp $(LIBS) $(LIBS_OPENCV) -D_DOSL_ALGORITHM=$(DOSL_ALGORITHM)


# --------------------------------------------

all_tests: cv_test json_test

.PHONY: cv_test
cv_test:
	@printf "\nNow compiling '$@'...\n"
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) -o $@ $@.cpp $(LIBS) $(LIBS_OPENCV)

.PHONY: json_test
json_test:
	@printf "\nNow compiling '$@'...\n"
	$(CC) $(CFLAGS) $(WARNS) $(INC_LOCAL) -o $@ $@.cpp $(LIBS)

# --------------------------------------------

all: all_tests all_examples

clean:
	rm $(prog)
	rm $(prog).o

