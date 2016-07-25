% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

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
addpath(genpath(workdir))
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;

sr = 1000 ;        % (Hz) sampling rate
N = 150 ;          % stimulus dimensionality 
ds = 2 ;           % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;       % delta t of stimulus 

order   = 2 ;      % order of MNE model to fit
njack   = 4 ;      % number of jackknives to run (also determines the size of each jackknives)
                   % this is used to split the training part of the data

rep = 500 ;        % number of repeats used to compute the null distribution of covariance eigenvalues

a = 0 ;            % significance threshold for eigenvalues of J


nJK = 5 ;          % number of jackknives used to test the prediction of the model
for i = 1:7
    for iJK = 1:nJK
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ; 
    
        S = S/N ;
    
        S = S - repmat(mean(S),[T,1]);
        S = S./repmat(std(S),[T,1]);
    
        
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
            mne_model(j,:) = MNEfit_WhiskerData(Sjfit, Rjfit, Sjtst, Rjtst, order) ;

            disp(['cell: ' Names{i} ', validation jackknife: ' num2str(iJK) ', fitting jackknife: ' num2str(j) ]) ; 

        end
        save([datadir 'VPM_cell_' Names{i} '_mne_model_JK_' num2str(iJK) '.mat'],'mne_model','N') ;
    end
end
for i = 1:7
    for iJK = 1:nJK
        % loop over jackknives used to test predictions

        load([datadir 'VPM_cell_' Names{i} '_mne_model_JK_' num2str(iJK) '.mat'])

        % averaging over jackknives used for fitting
        m = mean(mne_model,1) ;

        % extracting the parameters a, h, and J returned by the MNE
        % fitting routine. See definitions in Equation 32 and Table 2 
        A = m(1) ;
        H = m(2:N+1) ;
        J = reshape(m(N+2:end),N,N) ;
    
        % eigenvalues and eigenvectors of the matrix J (Equations 34
        % and 35)
        [vJ,eJ] = eig(J) ; 
        [eJ,iJ] = sort(diag(eJ)) ;
        vJ = vJ(:,iJ) ;
    
        % variable that will hold null eigenvalue distribution 
        eJnull = zeros(1,N*rep) ;     
    
        % diagonal and off-diagonal elements of J are shuffled
        % separately
        dJ = diag(J) ;
        oJ = reshape(triu(J,1),1,N*N) ;
        [~,iJ] = sort(abs(oJ),'descend') ;
        iJ = iJ(1:N*(N-1)/2) ;
        oJ = oJ(iJ) ;
    
        % shuffling matrix elements of J and computing eigenvalues is
        % repeated rep times
        for ir = 1:rep
            dinull = randperm(N) ;
            oinull = randperm(N*(N-1)/2) ;
            Jnull  = zeros(N) ;
            Jnull(iJ(oinull)) = oJ ;
            Jnull = Jnull + Jnull' + diag(dJ(dinull)) ;
           eJnull((ir-1)*N+(1:N)) = eig(Jnull) ;
        end
        max_enull = prctile(eJnull(:),100-a) ;  % upper bound of null distribution  
        min_enull = prctile(eJnull(:),a)  ;     % lower bound of null distribution
    
        isig = [find(eJ>max_enull) ; find(eJ<min_enull)] ;
        nsig = length(isig) ;                   % number of significant quadratic features
    
        if nsig>0 
            [~,iis] = sort(abs(eJ(isig)),'descend') ;
            isig = isig(iis) ;                  % orders the quadratic features according to the absolute value of their corresponding eigenvalues   
        end

        save([datadir 'VPM_cell_' Names{i} '_mne_JK_' num2str(iJK) '.mat'],'mne_model','N','nsig','J','H','A','isig','eJ','vJ','eJnull','max_enull','min_enull','rep') ;
    end
end