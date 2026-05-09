#!/bin/bash

scenarii=("Reference" "pDivided5" "pMultiplied5" "nDivided5" "nMultiplied5" "CaseControl1:7" "Equalp60" "Equalp240" "MainMO2SmallestOmics" "MainMO2LargestOmics" "HighMainMO" "HighConfMO" "HighConfMOIntersect" "HighConfSO" "HighConfSOIntersect" "OverlapAcrossEffects" "HighFractSignalFeat" "Noise" )

for sim in ${scenarii[@]}
do
	for rep in {1..40}
	do
		pimkl -fd ${sim}"view1Rep"${rep}.csv -nd ${sim}"view1Rep"${rep} -fd ${sim}"view2Rep"${rep}.csv -nd ${sim}"view2Rep"${rep} -fd ${sim}"view3Rep"${rep}.csv -nd ${sim}"view3Rep"${rep} --modelname EasyMKL "interaction"${sim}"Rep"${rep}.csv network "featureset"${sim}"Rep"${rep}.gmt genes preprocess output "labels"${sim}"Rep"${rep}.csv 0.2 5 5 
	done
done
