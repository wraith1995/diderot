#!/bin/bash
for i in {0..118}
do
    python3 runData.py $i
    ps -ef | grep py | awk '{print $2}' | xargs -I{} kill -9 {}
done
