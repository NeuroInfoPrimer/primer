function sim_coupled_GLM_L1(wd, id, nRep, l, maxBins)
    if nargin < 5
        maxBins = 1e20;
    end

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
    nB = min(size(proc_withheld.stim, 1), maxBins);
    method = 'spg';
    
    %Load fits (ggs_cpl, lambdas)
    load([wd '/all_units_network_L1_method_' method '.mat']);
    rng('shuffle')
    time_limit = 2400;
    stim = proc_withheld.stim(1:nB,:);
    stim = stim/p;
    
    for i = 1:nU
        ggs_cpl{l,i}.ihbas2 = ggs_cpl{l,i}.ihbas;
    end
    simstruct = makeSimStruct_GLMcpl(ggs_cpl{l,:});
    %Simulation with test stim
    Tt = size(stim,1);
    for i = 1:nU
        Rt_glm{l,i} = zeros(nRep,Tt);
    end
    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk] = simGLM_monkey(simstruct, stim, time_limit);
        for i = 1:nU
            Rt_glm{l,i}(ir, ceil(iR_glm{i})) = Rt_glm{l,i}(ir, ceil(iR_glm{i}))+1;
        end
        save([wd '/GLM_coupled_simulation_L1_method_' method '_lambda_' num2str(lambdas(l)) '_ID_' num2str(id) '.mat'], 'Rt_glm', 'nRep', 'ir');
    end
end