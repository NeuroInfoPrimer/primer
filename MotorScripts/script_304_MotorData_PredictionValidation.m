%Set to working directory
wd = '../MotorData/';

%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%
lambdas = [.1 .3 1 3 10 30 100 300];
global RefreshRate;
RefreshRate = 100;              %Stimulus refresh rate
ds = 0.001;                     %Spike time resolution
dt = ds*RefreshRate;            %Time scale spike times are resovled at,
                                % relative to stim timescale
datafile = 'mabel.mat';
frames = 80;                    %no. stim frames 
binsize = 1/RefreshRate;
nRep = 747;                     %no. sim repetitions
standardize = 0;
[proc, proc_withheld] = preprocess([wd datafile], binsize, dt, frames, standardize);    

nU = size(proc.spiketrain,2);
logl_glm = zeros(length(lambdas), nU);

%Load data from the other runs and add it to Rt_glm...
for l = 1:length(lambdas)
    l
    lambda = lambdas(l);
    r1 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_1.mat']);
    r2 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_2.mat']);
    r3 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_3.mat']);
    r4 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_4.mat']);
    r5 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_5.mat']);
    r6 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_6.mat']);
    r7 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_7.mat']);
    r8 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_8.mat']);
    r9 = load([wd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_9.mat']);

    Rt_glm = {};
    for idx = 1:nU
        Rt_glm{idx} = sum(r1.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r2.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r3.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r4.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r5.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r6.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r7.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r8.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r9.Rt_glm{l,idx},1);
    end
    nB = size(Rt_glm{1}, 2);

    for i = 1:nU
        Rt = proc_withheld.spiketrain(1:nB,i);
        Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
        if size(Rt_glm{i},1)==1
            Rt_glm{i} = Rt_glm{i}';
        end
        %Compute log-likelihood:
        logl_glm(l, i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(1/RefreshRate)) ;
    end   
    save([wd '/preprocessed_networkglm_sims_lambda_' num2str(lambda) '.mat'], 'proc_withheld', 'nU', 'Rt_glm', 'RefreshRate')

    %coh_out = ['coherence_lambda_' num2str(lambda)];
    %jackknifecoherence(wd, ['/preprocessed_networkglm_sims_lambda_' num2str(lambda) '.mat'], coh_out)
end

%Compare likelihood to uncoupled likelihood:
logl_glm_uncoupled = [];
for idx = 1:nU
    uncoupled = load([wd '/GLM_cell_simulation_' num2str(idx) '.mat']);
    logl_glm_uncoupled(idx) = uncoupled.logl_glm;
end

save([wd '/script_304_MotorData_PredictionValidation.mat'])
