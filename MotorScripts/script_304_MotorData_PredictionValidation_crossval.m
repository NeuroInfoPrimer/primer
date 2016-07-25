%Note: This script demonstrates how 
%      The particulars of how it is distributed over computers/CPUs depends on 
%      the user's setup and resources. 

%Note: Directory containing output of script_303_MotorData_SimCoupledGLM_crossval
%      This will depends on how the user implements 
%inputd = '../MotorData/';
inputd = '../MotorCoupledSims/';

%Set to working directory
wd = '../MotorData/';
fn_out = 'GLM_coherence_compare';

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
nfolds = 5;
nU = 9;
logl_glm = zeros(nfolds, length(lambdas), nU);
logl_glm_null = zeros(nfolds, nU);

for fold = 1:nfolds
    display(['Loading data from fold ' num2str(fold)])
    [proc, proc_withheld] = preprocess_crossval([wd datafile], binsize, dt, frames, fold, nfolds, standardize);
    %[proc, proc_withheld] = preprocess([wd datafile], binsize, dt, frames, standardize);    
    nB = size(proc_withheld.grip,1);
    %Load data from the other runs and add it to Rt_glm...
    for l = 1:length(lambdas)
        lambda = lambdas(l);
        display(['lambda: ' num2str(lambda)])
        %GLM_coupled_simulation_L1_method_spg_lambda_0.3_ID_611_fold_3.mat
        Rt_glm = {};
        for j = 1:nU
            Rt_glm{j} = zeros(1,nB);
        end
        for idx = 1:nRep
            if exist([inputd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_' num2str(idx) '_fold_' num2str(fold) '.mat'], 'file')
                rs = load([inputd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_' num2str(idx) '_fold_' num2str(fold) '.mat']);
                for j = 1:nU
                    Rt_glm{j} = Rt_glm{j}+rs.Rt_glm{l,j};
                end
            else
                display(['Not found: ' inputd 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_' num2str(idx) '_fold_' num2str(fold) '.mat'])
            end
        end
        nB = size(Rt_glm{1}, 2);    
        for i = 1:nU
            Rt = proc_withheld.spiketrain(1:nB,i);
            Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
            if size(Rt_glm{i},1)==1
                Rt_glm{i} = Rt_glm{i}';
            end
            %Compute log-likelihood:
            logl_glm(fold, l, i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(1/RefreshRate)) ;
        end   
        save([wd '/preprocessed_networkglm_sims_lambda_' num2str(lambda) '_fold_' num2str(fold) '.mat'], 'proc_withheld', 'nU', 'Rt_glm', 'RefreshRate')
        %coh_out = ['coherence_lambda_' num2str(lambda)];
        %jackknifecoherence_crossval(wd, ['/preprocessed_networkglm_sims_lambda_' num2str(lambda) '_fold_' num2str(fold) '.mat'], coh_out)
    end

    %%%%%%%%%%%%%%%%%%%%
    %Fit the null model%
    %%%%%%%%%%%%%%%%%%%%
    for i = 1:nU
        %Fit the firing rate to the training data
        meanrate = mean(proc.spiketrain(:,i))*ones(1,nB);
        %On the test data evaluate the likelihood
        logl_null = mean(Rt.*log(meanrate')-(meanrate')*(1/RefreshRate)) ;
        logl_glm_null(fold, i) = logl_null;
    end
end

%Compare likelihood to uncoupled likelihood:
logl_glm_uncoupled = [];
for fold = 1:nfolds
    for idx = 1:nU
        uncoupled = load([wd '/GLM_cell_simulation_' num2str(idx) '_fold_' num2str(fold) '.mat']);
        logl_glm_uncoupled(fold, idx) = uncoupled.logl_glm;
    end
end

mu_logl_glm = squeeze(mean(logl_glm,1));
mu_logl_glm_unc = squeeze(mean(logl_glm_uncoupled,1));
std_logl_glm = squeeze(std(logl_glm,1));
std_logl_glm_unc = squeeze(std(logl_glm_uncoupled,1));

save([wd '/script_304_MotorData_PredictionValidation_crossval.mat'])
