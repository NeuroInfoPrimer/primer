% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the fitting results of
% script_105_RetinaData_MNE_fitting.m and averages over the jackknives used
% for fitting to give the final MNE model. It also computes random matrices
% with the same statistics as the actual matrix J of second order features 
% to determine which are significant.

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

stim_length = {'short2','short3','long'} ;

% number of repeats used to compute the null distribution of eigenvalues of J
rep = 500 ;

% significance threshold for eigenvalues of J
a = 0 ;

for icell = 3:3
    for iL = [1 3]
        % loop over jackknives used to test predictions
        for iJK = 1:nJK
            
            load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            load([datadir 'Retina_cell_' num2str(icell) '_mne_model_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            
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
                isig = isig(iis) ;                   % orders the quadratic features according to the absolute value of their corresponding eigenvalues
            end
            save([datadir 'Retina_cell_' num2str(icell) '_mne_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'mne_model','N','NT','NX','nsig','J','H','A','isig','eJ','vJ','eJnull','max_enull','min_enull','rep') ;
        end
    end
end