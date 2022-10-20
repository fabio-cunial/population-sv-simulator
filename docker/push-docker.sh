#!/bin/bash
#
cp ../simulateHaplotypes/SimulateHaplotypes.java .
cp ../buildModel/BuildModel.java .
docker build -t fcunial/simulation .
docker push fcunial/simulation
rm -f SimulateHaplotypes.java BuildModel.java
