INTRODUCTION TO CODE

The MATLAB code is organized to have relatively few calls to auxiliary functions. The intention is to allow the reader/user to easily follow the processing and transformations the data undergoes from its raw form to the figures shown in the text, and relate specific lines of code to the equations in the text. We recommend using MATLAB's "cell mode" or "debug mode" to follow this step by step.

Each data set corresponds to a series of scripts: 
100's --- Retina
200's --- Thalamus (whisker)
300's --- Cortex (motor)

Each series is divided into two: "script" and "makefig". "script" files do the processing and analysis of the raw data and save the results. "makefig" files load these results and generate plots. Therefore it is necessary to first run the "script" for each analysis method and then run the "makefig". "makefig" files that produce figures which include results from multiple analysis methods will require all the corresponding analysis "script" files to be run beforehand. 

Except for the code that fits the Generalized Linear Model (GLM) and the L1 minimization mode (for fitting network GLM's), none of the code requires installation. 
To install the GLM we refer the user to the documentation of Jonathan Pillow's GLM package. 
To install the L1 minimization code we refer the user to Mark Schmidt's software website. (links to these appear in the Implementation section of the paper)
The rest of our code can be run after putting the data-files we supply in the appropriate directories. 

The raw data for each set is saved in slightly different formats. These are explained in the scripts that load the raw data files and do the necessary pre-processing (for example, script_101_RetinaData_BuildStimulus.m) and generate *.mat files that contain the stimulus (an N_X by N_M matrix) and the response (an N_M dimensional vector) for the training set (the portion of the data that is used to fit the model) and the test set (the portion of the model that is used to evaluate the predictions).

INSTALLATION

The code, data and saved results of the analysis are organized in seven directories

RetinaScripts
RetinaData
WhiskerScripts
WhiskerData
MotorScripts
MotorData
OtherScripts

1. To create these, simply unpack the file AnalysisOfNeuronalSpikeTrainsPrimer.zip under the desired directory.

2. Then copy the data files to the appropriate folder:

  2.1 Retina: 
  	  whitenoise.raw (checker board stimulus file)
  	  
  	  whitenoisec1.isk
  	  ...
  	  whitenoisec53.isk (spike trains of 53 cells)
  	  
  	  RetinaCellParameters_long.mat
  	  RetinaCellParameters_short2.mat
  	  RetinaCellParameters_short3.mat (coordinates of relevant window within the entire stimulus for each cell and corresponding time-lags for three choices of stimulus dimensionality)
  	  
  2.2 Whisker:
  	  VPM_cell_37_A2.mat
  	  VPM_cell_46_BC.mat
  	  VPM_cell_67_C4.mat
  	  VPM_cell_83_E2.mat
  	  VPM_cell_88_E1.mat
  	  VPM_cell_92_D2.mat
  	  VPM_cell_93_C4.mat (whisker position and corresponding spike trains for 7 rat thalamic cells. Description of data structure is in a separate file: Notes.txt)

  2.3 Motor:
	  mabel.mat
  	  
3. Add the folders where the auxiliary functions are saved to MATLAB's search path. These include the Chronux package, histcn.m, and the directories MotorGLM and L1Group in the MotorScripts directory.

4. Ensure the L1Group mex files have been compiled by running MotorScripts/L1Group/mexAll.m

5. Finally, install Jonathan Pillow's GLM package. 

LIST OF MAIN SCRIPTS

script_101_RetinaData_BuildStimulus.m
script_102_RetinaData_STA.m
script_103_RetinaData_STC_significance.m
script_104_RetinaData_STC_model.m
script_105_RetinaData_MNE_fitting.m
script_106_RetinaData_MNE_model.m
script_107_RetinaData_GLM.m
script_108_RetinaData_Prediction.m
script_109_RetinaData_Validation.m

makefig_101_RetinaData_StimulusPCA.m
makefig_102_RetinaData_STA.m
makefig_104_RetinaData_STC.m
makefig_105_RetinaData_MNE.m
makefig_106_RetinaData_GLM.m
makefig_107_RetinaData_PredictionValidation.m

script_201_WhiskingData_ExtractWhisking.m
script_202_WhiskingData_BuildStimulus.m
script_203_WhiskingData_STA.m
script_204_WhiskingData_WSTA.m
script_205_WhiskingData_STC.m
script_206_WhiskingData_WSTC.m
script_207_WhiskingData_MNE.m
script_208_WhiskingData_GLM.m
script_209_WhiskingData_TuningCurve.m
script_210_WhiskingData_Prediction.m
script_211_WhiskingData_Validation.m

makefig_201_WhiskingData_StimulusPCA.m
makefig_202_WhiskingData_IntervalDistributions.m
makefig_203_WhiskingData_WhiskingAutoCorrelation.m
makefig_204_WhiskingData_STA.m
makefig_205_WhiskingData_STC.m
makefig_206_WhiskingData_WSTA_WSTC.m
makefig_207_WhiskingData_MNE.m
makefig_208_WhiskingData_GLM.m
makefig_209_WhiskingData_TuningCurve.m
makefig_210_WhiskingData_PredictionValidation.m

script_301_MotorData_UncoupledGLM.m
script_302_MotorData_CoupledGLM.m
script_303_MotorData_SimCoupledGLM.m
script_304_MotorData_PredictionValidation.m

makefig_301_MotorData_Trials.m
makefig_302_MotorData_UncoupledGLM.m
makefig_303_MotorData_CoupledGLM.m
makefig_304_MotorData_PredictionValidation.m

script_301_MotorData_UncoupledGLM_crossval.m**
script_302_MotorData_CoupledGLM.m_crossval.m**
script_303_MotorData_SimCoupledGLM_crossval.m**
script_304_MotorData_PredictionValidation_crossval.m**
makefig_304_MotorData_PredictionValidation_crossval.m**

LIST OF OTHER FUNCTIONS USED

ReadFramev2.m
get_slow_var.m
phase_from_hilbert.m
herrorbar.m
addpath_recurse.m
subplot_pos.m
histcn.m

**These scripts perform the five-fold cross-validation used to select the optimal
  regularization for the coupled model, and the standard error estimates for the likelihood
  and mean coherence plots of Figure 14. script_303_MotorData_SimCoupledGLM_crossval
  is very computationally expensive, and will need to be modified by the user to 
  run on multiple CPUs/computers as availablilty dictates. These scripts are 
  not intended to run out of the box, but are instead provided as a template. 
