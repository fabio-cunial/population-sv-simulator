#!/bin/bash
#
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l SimulateHaplotypes.wdl -i inputs/SimulateHaplotypes-i300-noPav.json
java -jar ${WOMTOOL_PATH} validate -l JointCalling.wdl
java -jar ${WOMTOOL_PATH} validate -l AnnotateVCFs.wdl -i inputs/AnnotateVCFs-i300-noPav.json
java -jar ${WOMTOOL_PATH} validate -l PerformanceMatrices.wdl -i inputs/PerformanceMatrices_svTypes.json
