# uLocalize
Matlab package for single molecule analysis
Matlab and Python package for simulating transaltion of single mRNAs with bursting dynamics
System Requirements: Tested on MATLAB R2019B, requires MATLAB's Image Processing Toolbox. Parallel Computing Toolbox recommended.
Installation Guide: Download files and add to path on MATLAB, with all subfolders included. Installation should take less than 5 minutes. 

Demo: Examples can be found in example folder. Driver scripts included to run uLocalize in 2D and 3D versons. Input appropriate folder paths and filenames. Sample output also included in Results folder for 2D and 3D uLocalize. Run time can take 1-15 minutes per image on regular desktop computer. Parallel computing recommended. 

Instructions: Sigma values for single RNA can be calibrated by first using GaussianFit3D algotithm to determine sigmaxy and z values. Then, these values can be plugged in to use the GaussianMask algorithm. uLocalize3D was used for calculation of single mRNA. uLocalize 2D was used for calculation of integrated intensity of RNA granules for images taken on our wide-field microscope. 
