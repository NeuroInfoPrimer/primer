lambdas = [.1 .3 1 3 10 30 100 300];
wd = '../MotorData/';

%Simulate coupled GLM with range of lambda values
%This takes a long time, so split into many runs of a smaller number of
%repetitions to divide between CPUs/computers if necessary. 

%These runs should be moved to the one directory, and combined and analyzed afterwards
%by collapse_sim_coupled_GLM_L1.m
for idx = 1:length(lambdas)
	sim_coupled_GLM_L1(wd, 1, 83, idx);
	sim_coupled_GLM_L1(wd, 2, 83, idx);
	sim_coupled_GLM_L1(wd, 3, 83, idx);
	sim_coupled_GLM_L1(wd, 4, 83, idx);
	sim_coupled_GLM_L1(wd, 5, 83, idx);
	sim_coupled_GLM_L1(wd, 6, 83, idx);
	sim_coupled_GLM_L1(wd, 7, 83, idx);
	sim_coupled_GLM_L1(wd, 8, 83, idx);
	sim_coupled_GLM_L1(wd, 9, 83, idx);
end