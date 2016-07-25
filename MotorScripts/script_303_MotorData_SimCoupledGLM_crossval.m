lambdas = [.1 .3 1 3 10 30 100 300];
wd = '../MotorData/';

%Simulate coupled GLM with range of lambda values

%Note: this takes a long time, so must be split into many runs of a smaller number of
%      repetitions to divide between CPUs/computers as necessary. 
%      These runs should be moved to the one directory and combined and analyzed
%      afterwards.
%      
%      This script demonstrates how such a run would look like if ran serially
%      The particulars of how it is distributed over computers/CPUs depends on 
%      the user's setup and resources. 

running = [1,2,3,4,5];
for fold = running
	for idx = 1:length(lambdas)
		sim_coupled_GLM_L1_crossval(wd, 1, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 2, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 3, 83, idx, fold, nfolds);
	end
end

for fold = running
	for idx = 1:length(lambdas)
		sim_coupled_GLM_L1_crossval(wd, 4, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 5, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 6, 83, idx, fold, nfolds);
	end
end

for fold = running
	for idx = 1:length(lambdas)
		sim_coupled_GLM_L1_crossval(wd, 7, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 8, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 9, 83, idx, fold, nfolds);
	end
end