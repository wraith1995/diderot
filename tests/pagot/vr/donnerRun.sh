LD_LIBRARY_PATH=/home/teocollin/teem/build/bin
export LD_LIBRARY_PATH
make
/usr/bin/g++ -std=gnu++11 -shared -o evalProg.so evalProg.o /home/teocollin/diderot/lib/diderot-rt-seq-debug.o  -L/home/teocollin/teem/build/bin -lteem
