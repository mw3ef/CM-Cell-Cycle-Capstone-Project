# CM-Cell-Cycle-Capstone-Project
This is all the code and documentation for the CM cell cycle capstone project. The goal for the CM cell cycle model is to predict reaction rate constants in order to determine the concentration of cells moving from one cell cycle phase to another. The model predicts the reaction rates through parameter estimation by using experimental data that has been stained with DAPi and ki67.

There are two python files in this documentation: thresholding.py and paramEst.py

## Thresholding.py
Thresholding.py runs through an inputted csv file of the raw DAPI and ki67 intensity values for the cells. It runs these values against the thresholded values in order to tabulate the cell count percentage for each cell cycle phases of G0, G1, S, G2/M, Aneuploidy, and G0(4c). Outputs scatter plots of the experimental data and bar graphs after tabulating the cell counts.

## paramEst.py
Refers to Thresholding.py to get the cell count percentage values. Uses least squares optimization to predict the reaction rate constants. Outputs the parameter estimation graphs for each of the reaction rate constants and the concentration versus time curves.
