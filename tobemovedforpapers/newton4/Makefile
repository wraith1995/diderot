s#config vars:
SHELL = /bin/sh
PYTHON = python3
PYTHONINCLUDE = -I/usr/include/python3.6m -I/usr/include/python3.6m
DIDEROTFEMPATH = /home/teocollin/gitcode/diderot
TEEMPATH = /home/teocollin/installs/teem/teem-build
TARGET = evalProg
NAMESPACE = $(TARGET)
CC = clang
CXX = clang++ -std=gnu++11
PYTHONDATA = runData.py
PYTHONSTAT = runStatic.py

DIDEROTFEM = $(DIDEROTFEMPATH)/bin/diderotc
DIDEROTINCLUDE = $(DIDEROTFEMPATH)/include
DIDEROTLIB = $(DIDEROTFEMPATH)/lib
TEEMLIB = $(TEEMPATH)/include

DIDEROTFLAGS = --debug --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree --double --namespace=$(NAMESPACE) 
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

run: $(TARGET).so $(PYTHONDATA)
	$(PYTHON) $(PYTHONDATA)
saveRun: $(TARGET).so $(PYTHONDATA) $(PYTHONSTAT)
	$(PYTHON) $(PYTHONDATA) -s
	$(PYTHON) $(PYTHONSTAT)
clean:
	rm -f $(TARGET).cxx $(TARGET).h $(TARGET).o $(TARGET).so

runJson: $(PYTHONDATA)
	$(PYTHON) $(PYTHONDATA)
