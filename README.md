# ES_Ensembles
Belonging to Willcock et al. Ecosystem service model ensembles

This repository contains codes and example codes for the Manuscript : Ensembles of ecosystem service models can improve accuracy and indicate uncertainty. By Willcock, Hooftman et al. 
All codes are written in Matlab (v7.14.0.739 ) 
Codes are for information purposes only and in this combination copyrighted to:
 D.A.P. Hooftman, Lactuca: environmental Data Analysis and Modelling.
This repository contains:
A)	For creating Ensemble models and testing and storing accuracy
-	An example of a main steering function (Ensemble_main_example); running ones with model and validation sets names included can be acquired upon request
-	The actual statistics function (Accuracy_statistics.m)
-	The Ensemble function itself, including collection of individual model results to create the ensemble
B)	For comparing among ensembles (Figs 3 MS and Table SI1-3 & Table SI-1-4) and comarping among per validation point accuracy vs uncertainty (Fig 4 of MS). Using input from A above, relabeled outside an coding environment.
-	Ensembles_comparison.m
C)	For creating maps of mean, median and sem per gridcell among models and hotspot and cold area calculations
-	An example function (G_to_meanG.m ) ; running ones with model and validation sets names included can be acquired upon request
Danny Hooftman, Februari-12 2020
