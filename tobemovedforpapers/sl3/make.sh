#!/bin/bash
cp evalProg.json s1/evalProg.json
cd s1
make clean
make
cd ..
cp evalProg.json s2/evalProg.json
cd s2
make clean
make
cd ..
cp evalProg.json s3/evalProg.json
cd s3
make clean
make
cd ..
cp evalProg.json s4/evalProg.json
cd s4
make clean
make
cd ..
cp evalProg.json s5/evalProg.json
cd s5
make clean
make
cd ..
cp evalProg.json s6/evalProg.json
cd s6
make clean
make
cd ..
cp evalProg.json s7/evalProg.json
cd s7
make clean
make
cd ..
cp evalProg.json s8/evalProg.json
cd s8
make clean
make
cd ..
