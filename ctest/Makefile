#config vars:
SHELL = /bin/sh
PYTHON = python3
PYTHONINCLUDE = -I/usr/include/python3.6m -I/usr/include/python3.6m
DIDEROTFEMPATH = /home/teocollin/gitcode/diderot
TEEMPATH = /home/teocollin/installs/teem/teem-build
TARGET = justTypes
NAMESPACE = $(TARGET)
CC = clang
CXX = clang++ -std=gnu++11
PYTHONFILE = run.py

DIDEROTFEM = $(DIDEROTFEMPATH)/bin/diderotc
DIDEROTINCLUDE = $(DIDEROTFEMPATH)/include
DIDEROTLIB = $(DIDEROTFEMPATH)/lib
TEEMLIB = $(TEEMPATH)/include

DIDEROTFLAGS = --debug --log --json --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree --namespace=$(NAMESPACE) 
INCLUDES = -I$(DIDEROTINCLUDE) -I/usr/local/include -I$(TEEMLIB) -I$(PYTHONINCLUDE)
LDFLAGS = -L$(DIDEROTLIB) -L/usr/local/lib -L/opt/local/lib -L/usr/lib
LIBS = -lteem -lpng -lz -lbz2
CPPFLAGS = $(INCLUDES) -g O2 -fPIC
CXXFLAGS = $(INCLUDES) -g -fPIC

ifeq ($(PARALLEL),yes)
PTHREAD_FLAGS =	-D_THREAD_SAFE 
PTHREAD_LIBS =	
else
PTHREAD_FLAGS =
PTHREAD_LIBS =
endif


$(TARGET).cxx $(TARGET).h $(TARGET).o $(TARGET).so $(TARGET).json: $(TARGET).diderot $(DIDEROTFEM)
	$(DIDEROTFEM) $(DIDEROTFLAGS) $(TARGET).diderot

build: $(TARGET).so

run: $(TARGET).so $(PYTHONFILE)
	$(PYTHON) $(PYTHONFILE)

clean:
	rm -f $(TARGET).cxx $(TARGET).h $(TARGET).o $(TARGET).so $(TARGET).json $(TARGET).h.cffi $(TARGET)_wrap.o $(TARGET)_wrap.c _$(TARGET).so $(TARGET).py $(TARGET).log
