#!/bin/bash
#
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l SimulateHaplotypes.wdl -i toy-input.json 
java -jar ${WOMTOOL_PATH} validate -l JointCalling.wdl
