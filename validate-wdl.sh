#!/bin/bash
#
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l ./wdl/SimulatePopulation.wdl -i ./wdl/toy-input.json 
java -jar ${WOMTOOL_PATH} validate -l ./wdl/JointCalling.wdl
java -jar ${WOMTOOL_PATH} validate -l ./wdl/PavWrapper.wdl
