% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% The script uses routines written by Jeffery D. Fitzgerald in the 
% lab of Tatyana O. Sharpee. Further details can be found in
% Fitzgerald, J. D. Sincich, L. C. Sharpee, T. O. Minimal models of 
% multidimensional computations, PLoS Computational Biology, 7(3): 
% e1001111, 2011
% and 
% Fitzgerald, J. D Rowekamp, R. J. Sincich, L. C. Sharpee, T. O. Second 
% order dimensionality reduction using minimum and maximum mutual 
% information models, PLoS Computational Biology, 7(10): e1002249, 2011

% This script computes the Maximum Noise Entropy (MNE) LN model using a
% gradient ascent algorithm. In contrast to the STA and STC models, here
% the nonlinearity and the relevant subspace are computed together.
% Fitting is a convex problem with respect to the the training set, meaning
% that there are no local minima. However, convergence on the test set
% often leads to over-fitting. Therefore the trainig set is split into 4
% parts and fitting is repeated 4 times:
% in each repeated fitting ('jackknife') the gradient ascent runs and the
% likelihood function is computed on the training (3/4) and test (1/4)
% parts of the data. Because the problem is convex, the likelihood of the
% training portion always increases. The likelihood function computed on
% the test portion starts to decrease when the model is over-fit, so
% the algorithm stops at that point.

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
addpath(genpath(workdir)) ;

datadir = 'RetinaData/' ;

stim_length = {'short2','short3','long'} ;

order   = 2 ;   % order of MNE model to fit
njack   = 4 ;   % number of jackknives to run (also determines the size of each jackknives)
                % this is used to split the training part of the data
                
nJK = 5 ;       % number of jackknives used to test the prediction of the model

for icell = 3:3
    for iL = 1:3
        for iJK = 1:nJK
            load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;

            % the response is binarized because the MNE models the 
            % probability of a spike rather than a spike rate 
            R = sign(R) ; 
            
            % the MNE fitting routine produces a single parameter vector
            % for a second order model, there are N^2 + N + 1 parameters:
            % (see Equation 32 and Table 2)
            % a - scalar that allows model to satisfy zeroth order constraint
            % h - an N dimensional vector that allows model to satisfy first order constraint
            % J - an NxN matrix that allows model to satisfy second order constraint
            
            mne_model = zeros(njack,N*N+N+1) ;
            
            % number of frames in each fitting jackknife
            Tj = floor(T/njack);
            
            for j = 1:njack
                jtst = (1+(j-1)*Tj:j*Tj) ; % indices of part of data used for test
                jfit = setdiff(1:T,jtst) ; % indices of part of data used for training
                Sjtst = S(jtst,:) ;
                Rjtst = R(jtst) ;
                Sjfit  = S(jfit,:) ;
                Rjfit  = R(jfit) ;
                
                % Running the MNE fitting routine
                mne_model(j,:) = MNEfit_RetinaData(Sjfit, Rjfit, Sjtst, Rjtst, order) ;
                disp(['cell: ' num2str(icell) ', stimulus configuration: ' stim_length{iL} ', validation jackknife: ' num2str(iJK) ', fitting jackknife: ' num2str(j) ]) ; 
            end
            save([datadir 'Retina_cell_' num2str(icell) '_mne_model_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'mne_model','N') ;
        end
    end
end