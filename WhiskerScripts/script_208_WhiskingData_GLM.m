% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% this script fits GLM using a gradient descent algorithm on the likelihood
% function of the data given the parameters. It is based on the example
% scripts included in Jonathan Pillow's GLM code and uses routines from
% there to construct the spike history filter basis set (see
% below) and the fitting routine itself. 

% the stimulus history and spike history basis set is defined in the
% function makeFittingStruct_GLM_Whisker.m
% The basis set that is used for the stimulus history is composed from the 
% 12 leading principal components of the stimulus. In the text we discuss
% the reason for not using the standard "pixel basis" where every stimulus 
% coordinate is independent). The basis set that is used for the spike 
% history was unchanged from the original application in Jonathan Pillow's 
% code (Equations 46 and 47).

% a useful tutorial of how to use this code, including options that were
% not used here can be found in the examples included with the original
% implementation

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
addpath(genpath(workdir)) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
N = 150 ;
sr = 1000 ;
ds = 2 ; 

k = 12 ; % Number of PCs to be used in basis for temporal kernel
logl_glm = zeros(1,7) ;
rep = 500 ;
    
global RefreshRate ;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = sr/ds ; 

dt = ds/sr ; 
tp = dt*(-N+1:0)';  % time relative to spike of stim filter taps

nJK = 5 ;
for i = 1:7
    for iJK = 1:nJK
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        
        
        % rescaling the stimulus. We found that if the absolute value
        % of the stimulus is too large, fitting suffers from numerical
        % precision issues due to the non-saturating nonlinearity of
        % the GLM we implemented (exponential)
        S = S/N ;
        St = St/N ;

        % the GLM fitting routine treats this stimulus (single "spatial"
        % dimension, 150 temporal points) as a one dimensional stimulus.
        % the expansion to take into account the history is done within the
        % code. thus the stimulus is taken at the current time point
        S1 = S(:,1) ;  
        S1t = St(:,1) ;
        
        % spike times
        iR = find(R==1)' ;
    
        % the stimulus basis set is chosen to be the k leading principal
        % components of the stimulus (independent of spiking). these are
        % the k eigenvectors of the stimulus covariance with largest
        % eigenvalues.
        C = cov(S-repmat(mean(S,1),[T,1])) ;
        [vC,eC] = eig(C) ;
        [eC,iC] = sort(diag(eC),'descend') ;
        vC = vC(:,iC) ;

        % stimulus filter is initialized at the STA. 
        sta = reshape(simpleSTC(S1,iR,N),N,[]) ;


        %  Initialize parameters for fitting
        gg0 = makeFittingStruct_GLM_whisker(sta,dt,vC(:,1:k)) ;  % projects sta into basis for fitting k
        gg0.tsp = iR ;
        gg0.tspi = 1 ;

        % Compute logliklihood of initial params
        [logli0,rr0,tt] = neglogli_GLM(gg0,S1); 
        

        % Do ML estimation of model params
        opts = {'display', 'iter', 'maxiter', 20} ;
        [gg, negloglival] = MLfit_GLM(gg0,S1,opts) ; 
        
        Rt_glm = zeros(1,Tt) ;
    
        % Simulate the GLM rep times on the test set and average to
        % obtain a predicted firing rate.
        for ir = 1:rep
            [iR_glm, vmem,Ispk] = simGLM(gg, S1t) ;
            Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1 ;
        end
        Rt_glm = Rt_glm'/rep + 1e-8 ;

        save([datadir 'VPM_cell_' Names{i} '_glm_JK_' num2str(iJK) '.mat'],'gg','tp','Rt_glm') ; 
    end
end