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
datafile = '/mabel.mat';
nU = 9;                        %no. units
nS = 4;                         %no. stim components
frames = 80;                    %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
nRep = 747;                      %no. sim repetitions
standardize = 0;
[proc, proc_withheld] = preprocess([wd datafile], binsize, dt, frames, standardize);    
nB = size(proc.stim, 1);
trim = 1;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
offset = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%
%2 Fitting uncoupled GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%

ggs = {};
for icell = 1:nU
    disp(['Fitting unit ' num2str(icell)]);
    stim = proc.stim;
    stim = stim/p;
    resp = proc.spikes{icell};
    sptrain = proc.spiketrain(:,icell);

    stacked = proc.stacked;
    stacked = stacked/p;
    sta = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
    sta = reshape(sta,nF,[]);

    nspk(icell) = sum(sptrain);
    gg0 = makeFittingStruct_GLM_monkey(sta,dt_glm,Dt);
    gg0.tsp = resp';
    gg0.tspi = 1;

    opts = {'display', 'iter', 'maxiter', maxit};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, offset);
    ggs{icell} = gg;
end

%Save all
save([wd '/all_units.mat'], 'ggs', 'proc', '-v7.3');

%%%%%%%%%%%%%%%%%%%
%3 Simulate models%
%%%%%%%%%%%%%%%%%%%

maxBins = 1e20;
%maxBins = 20000;
nB = min(maxBins, size(proc_withheld.stim,1));

time_limit = 80;
units_conv = zeros(nU,1);
logl_glm_uncoupled = [];
rng('shuffle')
for icell = 1:nU
    %Simulation with test stim
    disp(num2str(icell));
    %Only within trial times...
    Tt = size(proc_withheld.stim(1:nB,:),1);
    Rt = proc_withheld.spiketrain(1:nB,icell);
    Rt_glm = zeros(1,Tt);
    nconverged = 0;
    gg = ggs{icell};

    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk, converged] = simGLM_monkey(gg, proc_withheld.stim(1:nB,:)/p, time_limit, 1);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
        nconverged = nconverged + converged;
    end
    units_conv(icell) = nconverged;
    Rt_glm = Rt_glm'/nRep + 1e-8;
    logl_glm = mean(Rt.*log(Rt_glm)-(Rt_glm)/RefreshRate) ;
    logl_glm_uncoupled(icell) = logl_glm;
    save([wd '/GLM_cell_simulation_' num2str(icell) '.mat'], 'Rt_glm', 'nconverged', 'logl_glm');
end

