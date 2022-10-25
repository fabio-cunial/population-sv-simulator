#!/bin/bash
#
cp ../src/*.java .
docker build -t fcunial/simulation .
docker push fcunial/simulation
rm -f *.java
