% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% this script fits GLM using a gradient descent algorithm on the likelihood
% function of the data given the parameters. It is based on the example
% scripts included in Jonathan Pillow's GLM code and uses routines from
% there to construct the stimulus and spike history filter basis set (see
% below) and the fitting routine itself. 

% the stimulus history and spike history basis set is defined in the
% function makeFittingStruct_GLM_Retina.m
% The basis set that is used for the stimulus history is the standard one
% (i.e. the "pixel basis" where every stimulus coordinate is independent).
% The basis set that is used for the spike history was unchanged from the
% original application in Jonathan Pillow's code (Equations 46 and 47).

% a useful tutorial of how to use this code, including options that were
% not used here can be found in the examples included with the original
% implementation

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
addpath(genpath(workdir)) ;
datadir = 'RetinaData/' ;

stim_length = {'short2','short3','long'} ;

dt = 1/30 ; % delta t of stimulus and response
sr = 1/dt ; % sampling rate

% the GLM provides a predicted spike train given the stimulus, and not a
% predicted firing rate. Each spike train generated given a single stimulus
% is different because the response depends on its own stochastic history.
% Therefore to obtain a predicted spike train we generate a large number
% (rep) of predicted spike trains and average them.
rep = 1000 ;

global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = sr ;

nJK = 5 ;
for icell = 3:3
    for iL = 2:3
        for iJK = 1:nJK
            
            load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            % rescaling the stimulus. We found that if the absolute value
            % of the stimulus is too large, fitting suffers from numerical
            % precision issues due to the non-saturating nonlinearity of
            % the GLM we implemented (exponential)
            S = S/N ;
            St = St/N ;

            % stimulus filter is initialized at the STA. 
            sta = S'*R/sum(R) - mean(S,1)' ;
            sta = reshape(sta,NX,NT)' ;
            
            % stimulus is reshaped to fit structure expected by GLM fitting routine 
            S1 = S(:,NX*(NT-1)+1:NX*NT) ;
            S1t = St(:,NX*(NT-1)+1:NX*NT) ;
            
            % spike times
            iR = find(R>0)' ;
            
            %  Initialize parameters for fitting
            gg0 = makeFittingStruct_GLM_Retina(sta,min(NX,NT),dt);
            gg0.tsp = iR ;
            gg0.tspi = 1 ;
             % Compute logliklihood of initial params
            [logli0,rr0,tt] = neglogli_GLM(gg0,S1);
            
            % Do ML estimation of model params
            opts = {'display', 'iter', 'maxiter', 100};
            [gg, negloglival] = MLfit_GLMbi(gg0,S1,opts); 
            
            Rt_glm = zeros(1,Tt) ;
            
            % Simulate the GLM rep times on the test set and average to
            % obtain a predicted firing rate.
            for ir = 1:rep
                [iR_glm, vmem,Ispk] = simGLM(gg, S1t) ;
                Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1 ;
            end
            Rt_glm = Rt_glm'/rep + 1e-8 ;
            
            save([datadir 'Retina_cell_' num2str(icell) '_glm_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'gg','Rt_glm') ;
        end
    end
end