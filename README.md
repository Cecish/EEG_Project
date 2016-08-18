# Application of data mining techniques on data obtained with an EEG-based Brain Computer Interface, for Hand Motor Rehabilitation
MSc Project

## Synopsis
Refined and tailored off-line EEG analysis, able to exploit characteristics inherent to each user.

## Dependencies
This EEG analysis uses the following toolboxes:
- EEGLAB: an open source MATLAB toolbox, designed for processing and analysing electro-physiological data. The latter needs to be downloaded from the Swartz Center for Computational Neuroscience's website (https://sccn.ucsd.edu/eeglab/), and then added the eeglab folder to the MATLAB path.
- Automatic Artifact Removal (AAR) toolbox v1.3, used in this project for removing EOG artefacts. The latter needs to be downloaded from http://germangh.github.io/aar/aardoc/aar.html. Please follow this link to install this toolbox (http://germangh.github.io/aar/aardoc/node4.html).
- mRMR MATLAB toolbox (accessible here: http://penglab.janelia.org/proj/mRMR/#publication): used in this project to select optimal EEG channels. There is no need to download it again, since it the latter available in the FeatureSelectionToolbox folder. You have have to ensure this folder is added to the MATLAB path.

## How to run
The EEG analysis proposed in this project is run with the command 'main'.
WARNING: You have approriately fill the 'input_file.txt' file for a successful run, as it for instance gathers information about the location of the acquired raw EEG signals, the EEG device used (EPOC or actiCAP)

## How to fill the input_file.txt file
- First line: Which hand should is (assumed to be) paralysed? (1: Right hand, 0: left hand)
- Second line: Data file path name
- Third line: Data file name
- Fourth line: Classifier used for the analysis (1: kNN and 2:SVM)
- Fifth line: Number of nearest neighbours (in case the kNN was chosen)
- Sixth line: Number of level for the discrete wavelet transform
- Seventh: Mother wavelet used for the discrete wavelet transform
- Eighth line: EEG system used (1: EPOC, 2: actiCAP)
