% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script computes the Spike Triggered Average (STA) and the STA LN
% model (i.e. the 1 dimensional linear nonlinear model where the single
% dimension - or feature - is the STA). The nonlinearity is computed using
% Bayes' rule as explained in the manuscript (Equations 8 and 9). 

% The difference from script_203_WhiskingData_STA.m is that here the
% stimulus is whitened (decorrelated) using the pseudoinverse of the 
% stimulus covariance matrix. Doing so introduces an extra algorithm 
% parameter - L - the order of the pseudoinverse (i.e. how many
% eigenvalues of the stimulus covariance are not set to 0).

% The STA model is computed in three steps, in the exact same way as before, 
% but using the whitened stimulus rather than the original stimulus: 
% 1. computation of the spike triggered avereage (i.e. reducing the 
%    dimensionality of the neuron's stimulus dependence to 1). (Equation
%    14)
% 2. computation of the spiking nonlinearity (Equations 8 and 9)
% 3. Building an LN model for a general stimulus by interpolating the 
%    spiking nonlinearity.
    
% Given the 1d LN model and the test data set, the log likelihood of the
%  model is also computed. 

clear ; 
workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

if exist('histcn','file') == 0
    error('this script uses histcn.m, the multidimensional histogram function. please add it to the search path and rerun') ;
end

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ; % Cell indices (used in file names)

sr = 1000 ;        % (Hz) sampling rate
N = 150 ;          % stimulus dimensionality 
ds = 2 ;           % downsampling fakctor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;       % delta t of stimulus 
Linv = 2:1:(N/5) ; % range of pseudoinverse orders (Equations 40 and 41)


nbns = 10 ;          % number of bins used to construct the stimulus distributions and estimate the nonlinearity
tp = -dt*(N:-1:1) ;  % the time vector associated with the stimulus history

nJK = 5 ;            % number of jack-knives 

for i = 1:7 
    % looping over all cells
    
    for iJK = 1:nJK
    % loop over jack-knives. 
    
    load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ; 
    
    iR = R==1 ;        % binary response vector 
    
    S = S/N ;          % rescaled training stimulus. 
                       %  S is a (T x N) matrix meaning that 
                       %  the 1st dimension corresponds to time
                       %  the 2nd dimension corresponds to the stimulus 'pixels'
                       
    St = St/N ;        % rescaled test stimulus.
        
    Ps = mean(iR) ;                      % Probability of a spike independent of the stimulus
    
    C = cov(S) ;                         % the stimulus covariance matrix 
    [vC,eC] = eig(C) ;                   % eigenvalues and eigenvectors of C
    [eC,iC] = sort(diag(eC),'descend') ; % sorting eigenvalues in decreasing order
    vC = vC(:,iC) ;                      % sorting eigenvectors to match eigenvalues
    
    logl_staL  = zeros(length(Linv),1) ;   % loglikelihood of each pseudoinverse model 
    
    for L = Linv 
        % looping over the randge of pseudo inverse orders
        
        eCp = [eC(1:L).^(-1/2) ; zeros(N-L,1) ] ;  % implementation of Equation 30:
                                                   % constructing a vector of the eigenvalues of pseudoinverse covariance of order L. 
                                                   % the first L eigenvalues are the inverse of the first L eigenvalues of C
                                                   % the last N-L eigenvalues are 0 
        Cpinv = vC'*diag(eCp)*vC ;                 % the pseudoinverse is given by the product of eigenvector matrix and a
                                                   % diagonal matrix with the eigenvalues on the diagonal (see Equations 40, 41)
        
        Cpinv = diag(diag(Cpinv))+triu(Cpinv,1)+triu(Cpinv,1)' ; % the matrix is 'forced' to be symmetric to avoid issues caused by roundoff numerical errors
        

        SL = S*Cpinv ; % SL is the prewhitened stimulus matrix (i.e. the stimulus times the pseudoinverse covariance)
        
        % from here, the STA and the STA model will be computed as before, using SL instead of S
        staL = mean(SL(iR,:),1) - mean(SL,1) ; 
        staL = staL' / norm(staL) ;
        
        zL = SL*staL ;
        zxL = max(abs(zL(:)))+1e-3 ;
        bns = linspace(-zxL,zxL,nbns) ;
        dbns = bns(2)-bns(1) ;
        ctr = bns(1:nbns-1)+dbns/2 ;
    
        PfL  = histcn(zL,bns) ;
        PfL  = PfL/(T*dbns) ;          
        
        PfsL = histcn(zL(iR),bns) ;
        PfsL = PfsL/(T*Ps*dbns) ;
        
        PsfL = PfsL./PfL*Ps ;
        PsfL(isnan(PsfL)) = 0 ;
        
        PsfL = PsfL' ;
        
        sta_modelL = @(x)interp1(ctr,PsfL,x,'pchip') ;
        sta_model_rectL = @(x) heaviside(sta_modelL(x)).*sta_modelL(x)+1e-8 ;
        Ps_modelL = sum(sta_model_rectL(SL*staL)) ;
        sta_model_rect_normL = @(x) sta_model_rectL(x)*Ps*T/Ps_modelL ;
        save([datadir 'VPM_cell_' Names{i} '_wsta_L' num2str(L) '_JK_' num2str(iJK) '.mat'],'sta_model_rect_normL','staL','tp','N','Cpinv','Ps','PfL','PfsL','ctr') ;
    
        StL = St*Cpinv ; % the prewhitened test stimulus 
        logl_staL(L-1) = mean(Rt.*log(sta_model_rect_normL(StL*staL))-sta_model_rect_normL(StL*staL)) ;
    end
    [~,Lopt] = max(logl_staL) ; % the optimal pseudoinverse order is the one that gives the best prediction on the test set
    Lopt = Lopt + 1 ;
    save([datadir 'VPM_cell_' Names{i} '_wsta_logl_JK_' num2str(iJK) '.mat'],'logl_staL','Lopt') ;
    end
end