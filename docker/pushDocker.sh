#!/bin/bash
#
cp ../src/*.java .
docker build --progress=plain -t fcunial/simulation .
docker push fcunial/simulation
rm -f *.java
