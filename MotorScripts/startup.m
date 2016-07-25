%Note: If MATLAB is started in MotorScripts directory then startup.m is run 
%      automatically. Otherwise, this file will need to be run to add the following
%      to the search path.

%Get the path of this file
fullpath = fileparts(mfilename('fullpath'));
cd(fullpath)

addpath([pwd '/MotorGLM']);
addpath_recurse([pwd '/L1Group/']);

%Note: add path to Pillow's GLM code here
%      not necessary if already added to MATLAB's search path
%      addpath_recurse('/home/lansdell/projects/neuroinf/GLM_Algorithm_Functions');
