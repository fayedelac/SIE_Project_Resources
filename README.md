# SIE_Project_Resources

SIE Project
Faye Delacrétaz
12.01.2024

This repository contains all files relevant to the research SIE project. This project was conducted at the RIVER lab, EPFL, during the autumn semester of 2023-2024.

-------

## MATLAB code
This folder contains all relevant files to the derivation of the network's hierarchy and the nodes' estimated states.
The code is an adapated version of the algorithm developped by Durighetto et al. (2023)

### data_real.csv
CSV file containing the real observation states of the nodes recorded during the three conducted surveys.
1 denotes the wet state, 0 not observed and -1 dry state

### data_complete.csv
CSV file containing the artifical data created in the frame of the research project.
This dataset is completelly filled meaning there are no 0 values.

### hierarchyReconstructionAnalysis.mlx
Main code calling all necessary function to compute and display the network's hierarchy and estimated states.

### getHierarchy_p.m
Functions necessary to compute the networks hierarchy as well as the nodes' estimated states for different percentage of observed data.

### plotHierarchyGraph.m
Function displaying the hierarchy graph.

### plotStatus
Function displaying the confusion matrix which presents the observed and estimated states of the nodes.

-------

## QGIS files
Folder containing all necessary files for the creation of the stream network on QGIS.