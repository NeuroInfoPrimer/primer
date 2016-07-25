% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script and script_104_RetinaData_STC_model.m compute the 
% Spike Triggered Covariance (STC) and the STC LN model (i.e. the 
% multi-dimensional linear nonlinear model where the dimesnions 
% or features - are the STA and STC) (Equations 18-24). The nonlinearity is 
% computed by averaging the response conditioned on the stimulus being in a
% specific bin (Equations 6 and 7). 

% The STC model is computed in five steps:
% (steps 1,2,3 in script_104_RetinaData_STC_significance.m and 
%  steps 4,5   in script_104_RetinaData_STC_model.m) 

% 1. computation of the spike triggered average. (Equations 14 and 15)
% 2. computation of the spike triggered covariance (Cs) and the matrix of 
%    covariance differences (dC) (Equations 18-22) 
% 3. finding significant STC dimensions by comparing to a null eigenvalue
%    distribution (Equation 23)
% 4. computation of the spiking nonlinearity (Equations 6 and 7)
% 5. building an LN for general stimulus by interpolating the spiking 
%    nonlinearity.
    
clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

stim_length = {'short2','short3','long'} ;

rep = 1000 ;       % number of repeats used to compute the null distribution of covariance eigenvalues

a = 0 ;            % level of significance required to call an STC feature significant
                   % if an eigenvalue of dC is smaller than the a-th percentile of the null eigenvalue 
                   % distribution or larger than the 100-a percentile, it will be determined to be significant

nJK = 5 ;          % number of jackknives 

for icell = 3:3        
    % results for cell 3 are presented in paper. to compute model for 
    % other cells change to a vector including all the cell
    % numbers you are interested in (from 1 to 53).
    for iL = 2:3
        % for the STA model we consider two stimulus configurations: 'long'
        % and 'short3'.

        for iJK = 1:nJK
            % loop over jack-knives. the STA model is estimated nJK = 5 times 

            load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            
            % 1. computation of the spike triggered avereage:
            % -----------------------------------------------
    
            sta = S'*R/sum(R) - mean(S,1)' ;    % implementation of Equations 14 and 15:
                                                % the 1st term on the RHS is the stimulus average conditioned on a spike
                                                % the 2nd term on the RHS is the unconditioned stimulus average 
    
            sta = sta / norm(sta) ;             % the STA is normalized to 1 
            
            % 2. computation of the spike triggered covariance:
            % -------------------------------------------------
    
            Cp = cov(S-repmat(mean(S,1),[T,1])) ;  % implementation of Equations 18:
                                                   % Cp is the covariance matrix of the underlying stimulus
                                             
            Cs = 1/mean(R)*cov(S.*repmat(sqrt(R),[1 N])-repmat(mean(S.*repmat(sqrt(R),[1 N])),[T 1])) ; 
                                                   % implementation of Equations 21:
                                                   % Cs is the spike triggered covariance matrix
                                               
            % implementation of Equations 20, 22: 
            % the STC dimensions are found by diagonalizing dC - the matrix
            % of covariance differences
            dC = Cs-Cp ;                           
            [vdC,edC] = eig(dC) ;                       
            [edC,idC] = sort(diag(edC)) ;           
            vdC = vdC(:,idC) ;
        
            % 3. computation of the null eigenvalue distribution:
            % ---------------------------------------------------                                 
            
            % an array of covariance difference eigenvalues 
            % computed by assuming no association of the spike train
            % and the stimulus. (Equation 23)
            
            % the eigenvalues of dC (edC) will be compared to the 
            % distribution of null covariance eigenvalues to dermine 
            % which are significant 
            edCnull = zeros(1,N*rep) ;             
            
            for ir = 1:rep                         % loop over number of repeats of shuffled spike trains 
         
                % random spike train obtained by shifting the real spike train by a random amount
                % with periodic boundary conditions. This preserves the possible burst structure and 
                % statistics of spike counts
                Rrnd = R(mod(randi(T)+(1:T),T)+1) ;  
                
                % the null STC matrix is computed in the same way Cs is but conditioned on the 
                % random spike train rather than the real one (Equation 23)
                Crnd = 1/mean(R)*cov(S.*repmat(sqrt(Rrnd),[1 N])-repmat(mean(S.*repmat(sqrt(Rrnd),[1 N])),[T 1])) ; 
                                             
                % p null eigenvalues are computed in each iteration, so in sum
                % the distributino of null eigenvalues will be constructed from 
                % rep x N eigenvalues.
                edCnull((ir-1)*N+(1:N)) = eig(Crnd-Cp) ;
            end
    
            max_enull = prctile(edCnull(:),100-a) ;  % upper bound of null distribution  
            min_enull = prctile(edCnull(:),a)  ;     % lower bound of null distribution
            vprct = linspace(a,100-a,N) ;
            prctenull = prctile(edCnull(:),vprct) ;  % percentiles of null distribution (all) 
            
            % the significant STC filters corrspond to eigenvalues of dC outside of the null distribution
            isig = [find(edC>max_enull) ; find(edC<min_enull)] ;
                                                      
            % number of significant STC features:
            nsig = length(isig) ;                
            
            % if there is at least 1 significant STC feature, then features
            % are ordered features according to the absolute value of their 
            % corresponding eigenvalues
            if nsig>0 
                [~,iis] = sort(abs(edC(isig)),'descend') ;
                isig = isig(iis) ;                      
            end
            disp(['cell: ' num2str(icell) ', stimulus configuration: ' stim_length{iL} ', validation jackknife: ' num2str(iJK) ]) ; 

            save([datadir 'Retina_cell_' num2str(icell) '_stcsig_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'vdC','edC','edCnull','max_enull','min_enull','isig','nsig','a','vprct','prctenull') ;
        end
    end
end