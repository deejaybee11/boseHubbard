EXEC = exec
CC = icpc
SRCDIR = src/
MKLROOT = /opt/intel/mkl
MKLINCLUDE = $(MKLROOT)/include
LOCALINCLUDE = /usr/local/include
LOCALLIB = /use/local/lib

SRCEXT := cpp
SRC_FILES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))

CFLAGS = -c -g -fvar-tracking  -traceback -Wall -DMKL_ILP64 -openmp -fast -O3 -xhost -ip -qopt-report=5 -fbuiltin -ipo -no-ftz -static-intel -qopt-report-phase=par,vec,openmp -std=c++11 -mkl=parallel -I$(MKLINCLUDE) -I$(LOCALINCLUDE)
LFLAGS = -L$(MKLROOT)/lib/intel64 -L$(LOCALLIB) -parallel -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lcfitsio -lm

O_FILES = $(SRC_FILES:.cpp=.o)

print-%  : ; @echo $* = $($*)

$(EXEC): $(O_FILES)
	$(CC) -o $@ $(O_FILES) $(LFLAGS)
	
%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
