% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script computes the Spike Triggered Covariance (STC) and the STC LN
%  model (i.e. a >1 dimensional linear nonlinear model where the first 
%  dimension - or feature - is the STA and the the second feature is the STC).
% The nonlinearity is computed using
%  Bayes' rule as explained in the manuscript.

% The difference from script_205_WhiskingData_STC.m is that here the
%  stimulus is whitened (decorrelated) using the pseudoinverse of the 
%  stimulus covariance matrix. Doing so introduces an extra algorithm 
%  parameter - L - the order of the pseudoinverse (i.e. how many
%  eigenvalues of the stimulus covariance are not set to 0).

% The STC model is computed in five steps, in the exact same way as before, 
%  but using the whitened stimulus rather than the original stimulus: 
% 1. computation of the spike triggered average. (Equations 14 and 15)
% 2. computation of the spike triggered covariance (Cs) and the matrix of 
%    covariance differences (dC) (Equations 18-22) 
% 3. finding significant STC dimensions by comparing to a null eigenvalue
%    distribution (Equation 23)
% 4. computation of the spiking nonlinearity (Equations 8 and 9)
% 5. building an LN for general stimulus by interpolating the spiking 
%    nonlinearity.    

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
ds = 2 ;           % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;       % delta t of stimulus 
Linv = 2:1:(N/4) ; % range of pseudoinverse orders  


nbns = 10 ;          % number of bins used to construct the stimulus distributions and estimate the nonlinearity
nbns_plot = 25 ;

tp = -dt*(N:-1:1) ;  % the time vector associated with the stimulus history
rep = 500 ;          % number of repeats used to compute the null distribution of covariance eigenvalues

nJK = 5 ;            % number of jack-knives 

nspk = zeros(7,nJK) ;% a vector with the number of spikes for each cell 

a = 0.25 ;           % level of significance required to call an STC feature significant
                     % if an eigenvalue of dC is smaller than the a-th percentile of the null eigenvalue 
                     % distribution or larger than the 100-a percentile, it will be determined to be significant
for i = 1:7 
    % looping over all cells 
    
    for iJK = 1:nJK
        % loop over jack-knives.
    
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ; 
        nspk(i,iJK) = sum(R) ;
        iR = R==1 ;        % binary response vector 
    
    
    
        S = S/N ;          % rescaled training stimulus. 
                           %  S is a (T x p) matrix meaning that 
                           %  the 1st dimension corresponds to time
                           %  the 2nd dimension corresponds to the stimulus 'pixels'
                       
        St = St/N ;        % rescaled test stimulus.
        
        PsL = mean(iR) ;
    
        C = cov(S) ;                         % the stimulus covariance matrix 
        [vC,eC] = eig(C) ;                   % eigenvalues and eigenvectors of C
        [eC,iC] = sort(diag(eC),'descend') ; % sorting eigenvalues in decreasing order
        vC = vC(:,iC) ;                      % sorting eigenvectors to match eigenvalues
    
        logl_stcL  = -inf*ones(length(Linv),1) ;   % loglikelihood of each pseudoinverse model 
    
        for L = Linv 
            % looping over the randge of pseudo inverse orders
        
            eCp = [eC(1:L).^(-1/2) ; zeros(N-L,1) ] ; % constructing a vector of the eigenvalues of pseudoinverse covariance of order L. 
                                                      % the first L eigenvalues are the inverse of the first L eigenvalues of C
                                                      % the last N-L eigenvalues are 0 
            Cpinv = vC'*diag(eCp)*vC ;                % the pseudoinverse is given by the product of eigenvector matrix and a
                                                      % diagonal matrix with the eigenvalues on the diagonal (see Equations 40, 41)
                                                    
            Cpinv = diag(diag(Cpinv))+triu(Cpinv,1)+triu(Cpinv,1)' ; % the matrix is 'forced' to be symmetric to avoid issues caused by roundoff numerical errors

            SL = S*Cpinv ; % Sk is the prewhitened stimulus matrix (i.e. the stimulus times the pseudoinverse covariance)
            StL = St*Cpinv ; % Sk is the prewhitened stimulus matrix (i.e. the stimulus times the pseudoinverse covariance)
        
            % from here, the STA and the STC model will be computed as before, using SL instead of S  
            staL = mean(SL(iR,:),1) - mean(SL,1) ; 
            staL = staL' / norm(staL) ;
        
            CpL = cov(SL-repmat(mean(SL,1),[T,1])) ;  
        
            CsL = cov(SL(iR,:)-repmat(mean(SL(iR,:)),[nspk(i,iJK),1])) ; 
                                             
            dCL = CsL-CpL ;
            [vdCL,edCL] = eig(dCL) ;
            [edCL,idCL] = sort(diag(edCL)) ;
            vdCL = vdCL(:,idCL) ;
    
            edCnullL = zeros(1,N*rep) ;
        
            for ir = 1:rep                           
                % loop over number of repeats of shuffled spike trains 
            
                Rrnd = R(mod(randi(T)+(1:T),T)+1) ;
                iRrnd = Rrnd == 1 ;
            
                CrndL = cov(SL(iRrnd,:)-repmat(mean(SL(iRrnd,:)),[nspk(i,iJK),1])) ; 
        
                edCnullL((ir-1)*N+(1:N)) = eig(CrndL-CpL) ;
            end
    
            max_enullL = prctile(edCnullL(:),100-a) ;  % upper bound of null distribution  
            min_enullL = prctile(edCnullL(:),a)  ;     % lower bound of null distribution
            vprct = linspace(a,100-a,N) ;
            prctenull = prctile(edCnullL(:),vprct) ;   % percentiles of null distribution (all) 
    
            isigL = [find(edCL>max_enullL) ; find(edCL<min_enullL)] ;
                                    
            nsigL = length(isigL) ;                    % number of significant STC features
            if nsigL>0 
                [~,iis] = sort(abs(edCL(isigL)),'descend') ;
                isigL = isigL(iis) ;                   % orders the STC features according to the absolute value of their corresponding eigenvalues   
            end
        
            save([datadir 'VPM_cell_' Names{i} '_stcsig_L' num2str(L) '_JK_' num2str(iJK) '.mat'],'vdCL','edCL','edCnullL','max_enullL','min_enullL','isigL','nsigL','a','vprct','prctenull') ;
    
            if nsigL > 0
    
                stcL = vdCL(:,isigL) - staL*(vdCL(:,isigL)'*staL)' ;
                for iisig = 1:nsigL
                    stcL(:,iisig) = stcL(:,iisig)/norm(stcL(:,iisig)) ;
                end
                zL = zeros(T,2) ;                     % a [T x 2] array that holds the training stimulus projected on the STA (column 1) and on STC1 (column 2)
                zL(:,1) = SL*staL ;
                zL(:,2) = SL*stcL(:,1) ;
        
                ztL = zeros(Tt,2) ;                   % a [Tt x 2] array that holds the test stimulus projected on the STA (column 1) and on STC1 (column 2)
                ztL(:,1) = StL*staL ;
                ztL(:,2) = StL*stcL(:,1) ;
        
                % finding the range of z values to generate appropriate bins
                % (Equation 7)
               zx = max(abs(zL(:)))+1e-5 ;
            
                % cell array of 2 bin vectors (here they are identical, but they do not have to be)
                bns = cell(1,2) ;                   
                bns{1} = linspace(-zx,zx,nbns) ;
                bns{2} = linspace(-zx,zx,nbns) ; 
        
                % the size of each bin (for computing proper probability distributions)   
                dbns = bns{1}(2)-bns{1}(1) ;
            
                % same for high resolution bins used for plotting 
                bns_plot = cell(1,2) ;                 
                bns_plot{1} = linspace(-zx,zx,nbns_plot) ;
                bns_plot{2} = linspace(-zx,zx,nbns_plot) ; 
                dbns_plot = bns_plot{1}(2)-bns_plot{1}(1) ;
            
                PsL = nspk(i,iJK)/T ;            
                [PfL,ctr] = hist3(zL,bns) ;
                PfL = PfL/(T*dbns^2) ;     
        
                PfsL = hist3(zL(iR,:),bns) ;
                PfsL = PfsL/(nspk(i,iJK)*dbns^2) ;
        
                PsfL = PfsL./PfL*PsL ;
                PsfL(isnan(PsfL)) = 0 ;
            
        
                PfL_a  = histcn(zL(:,1),bns{1}) ;
                PfsL_a = histcn(zL(iR,1),bns{1}) ;
                PfL_c  = histcn(zL(:,2),bns{2}) ;
                PfsL_c = histcn(zL(iR,2),bns{2}) ;
                
                PfL_a_plot  = histcn(zL(:,1),bns_plot{1}) ;
                PfsL_a_plot = histcn(zL(iR,1),bns_plot{1}) ;
                PfL_c_plot  = histcn(zL(:,2),bns_plot{2}) ;
                PfsL_c_plot = histcn(zL(iR,2),bns_plot{2}) ;
                
                PfL_a  = PfL_a/(T*dbns) ;
                PfsL_a = PfsL_a/(T*PsL*dbns) ;
                PfL_a_plot  = PfL_a_plot/(T*dbns_plot) ;
                PfsL_a_plot = PfsL_a_plot/(T*PsL*dbns_plot) ;
                
                PfL_c  = PfL_c/(T*dbns) ;
                PfsL_c = PfsL_c/(T*PsL*dbns) ;
                PfL_c_plot  = PfL_c_plot/(T*dbns_plot) ;
                PfsL_c_plot = PfsL_c_plot/(T*PsL*dbns_plot) ;
                
                PsfL_a = PfsL_a./PfL_a*PsL ;
                PsfL_a(isnan(PsfL_a)) = 0 ;
                PsfL_c = PfsL_c./PfL_c*PsL ;
                PsfL_c(isnan(PsfL_c)) = 0 ;
                
                PsfL_a_plot = PfsL_a_plot./PfL_a_plot*PsL ;
                PsfL_a_plot(isnan(PsfL_a_plot)) = 0 ;
                PsfL_c_plot = PfsL_c_plot./PfL_c_plot*PsL ;
                PsfL_c_plot(isnan(PsfL_c_plot)) = 0 ;
                
                stc_modelL = @(x)interpn(ctr{1},ctr{2},PsfL,x(:,1),x(:,2),'linear') ;
                stc_model_rectL = @(x) heaviside(stc_modelL(x)).*stc_modelL(x)+1e-8 ;
                Ps_modelL = sum(stc_model_rectL(zL)) ;
                stc_model_rect_normL = @(x) stc_model_rectL(x)*nspk(i,iJK)/Ps_modelL ;
                
                save([datadir 'VPM_cell_' Names{i} '_wstc_L' num2str(L) '_JK_' num2str(iJK) '.mat'],'stc_model_rect_normL','staL','stcL','tp','N','PfsL','PsfL','PfL','nsigL','isigL','vdCL','PfsL_a','PfsL_a_plot','PsfL_a','PsfL_a_plot','PfL_a','PfL_a_plot','PfsL_c','PfsL_c_plot','PsfL_c','PsfL_c_plot','PfL_c','PfL_c_plot','bns','nbns','dbns','bns_plot','nbns_plot','dbns_plot','Cpinv') ; 
            
                Rt_stcL = stc_model_rect_normL(ztL) ;
                Rt_stcL(isnan(Rt_stcL)) = mean(Rt_stcL(~isnan(Rt_stcL))) ;
                logl_stcL(L-1) = mean(Rt.*log(Rt_stcL)-Rt_stcL) ;
            end
        end
        [~,Lopt] = max(logl_stcL) ; % the optimal pseudoinverse order is the one that gives the best prediction on the test set
        Lopt = Lopt + 1 ;
        save([datadir 'VPM_cell_' Names{i} '_wstc_logl_JK_' num2str(iJK) '.mat'],'logl_stcL','Lopt') ;
        
        disp(['cell: ' Names{i} ', validation jackknife: ' num2str(iJK) ]) ;
    end
end