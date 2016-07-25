%Set to working directory
wd = '../MotorData/';

%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%

global RefreshRate;
RefreshRate = 100;              %Stimulus refresh rate
ds = 0.001;                     %Spike time resolution
dt = ds*RefreshRate;            %Time scale spike times are resovled at,
                                % relative to stim timescale
datafile = 'mabel.mat';
nS = 4;                         %no. stim components
frames = 80;                    %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
standardize = 0;
[proc, proc_withheld] = preprocess([wd datafile], binsize, dt, frames, standardize);    
nU = size(proc.spiketrain, 2);
nB = size(proc.stim, 1);
trim = 1;
Dt = 20;
maxit = 10;
dt_glm = 0.1;
offset = 1;
method = 'spg';

%%%%%%%%%%%%%%%%%%%%%%%%%
%2 Fitting uncoupled GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%

ggs_cpl = {};
maxiter = 20;
lambdas = [.1 .3 1 3 10 30 100 300];
for l = 1:length(lambdas)
    lambda = lambdas(l);
    for icell = 1:nU
        disp(num2str(icell));
        stim = proc.stim;
        stim = stim/p;
        resp = proc.spikes{icell};
        sptrain = proc.spiketrain(:,icell);
        nicell = [(1:icell-1), (icell+1:nU)];
        %Add coupling to the other spike trains
        coupled = proc.spikes(nicell);
        for idx = 1:length(coupled)
            coupled{idx} = coupled{idx}';
        end
        stacked = proc.stacked;
        stacked = stacked/p;
        sta = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
        sta = reshape(sta,nF,[]);
        nspk(icell) = sum(sptrain);
        gg0 = makeFittingStruct_GLM_monkey(sta,dt,Dt);
        gg0.ihbas2 = gg0.ihbas;
        gg0.tsp = resp';
        gg0.tspi = 1;
        %Other spike trains
        gg0.tsp2 = coupled;
        %Add terms for other spike filters
        gg0.ih = zeros(size(gg0.ih,1),nU);
        opts = {'display', 'iter', 'maxiter', maxiter};
        [gg, negloglival] = MLfit_GLM_trim_L1_MLE(gg0,stim,opts,proc,trim, offset, lambda, method);
        ggs_cpl{l,icell} = gg;
    end
end
%Save all
save([wd '/all_units_network_L1_method_' method '.mat'], 'ggs_cpl', 'lambdas');

