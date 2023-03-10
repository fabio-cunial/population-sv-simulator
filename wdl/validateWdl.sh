#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l SimulateHaplotypes.wdl -i inputs/SimulateHaplotypes-i300-noPav.json
java -jar ${WOMTOOL_PATH} validate -l JointCalling.wdl
java -jar ${WOMTOOL_PATH} validate -l AnnotateVCFs.wdl -i inputs/AnnotateVCFs-i300-noPav.json
java -jar ${WOMTOOL_PATH} validate -l PerformanceMatrices.wdl -i inputs/PerformanceMatrices_svTypes.json
java -jar ${WOMTOOL_PATH} validate -l CallRealTrios.wdl -i inputs/CallRealTrios_noPav.json
java -jar ${WOMTOOL_PATH} validate -l ReadLengthDistribution.wdl -i inputs/private_ReadLengthDistribution.json
java -jar ${WOMTOOL_PATH} validate -l TriosCreateTruthVCFs.wdl -i inputs/private_TriosCreateTruthVCFs.json
java -jar ${WOMTOOL_PATH} validate -l TriosCreateTruthVCFs2.wdl -i inputs/private_TriosCreateTruthVCFs2.json
java -jar ${WOMTOOL_PATH} validate -l TriosPerformanceMatrices.wdl -i inputs/private_TriosPerformanceMatrices.json
java -jar ${WOMTOOL_PATH} validate -l TriosCreateSVLengthHistograms.wdl -i inputs/private_TriosCreateSVLengthHistograms.json
java -jar ${WOMTOOL_PATH} validate -l TriosDeleteLargeFiles.wdl -i inputs/private_TriosDeleteLargeFiles.json
java -jar ${WOMTOOL_PATH} validate -l TriosCoverageEffect.wdl -i inputs/private_TriosCoverageEffect.json
java -jar ${WOMTOOL_PATH} validate -l TriosCoverageEffectPerformanceMatrices.wdl -i inputs/private_TriosCoverageEffectPerformanceMatrices.json
