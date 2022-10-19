#!/bin/bash
#
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l SimulatePopulation.wdl -i toy-input.json 
java -jar ${WOMTOOL_PATH} validate -l JointCalling.wdl
java -jar ${WOMTOOL_PATH} validate -l PavWrapper.wdl
