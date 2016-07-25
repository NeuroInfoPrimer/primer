% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039


% This script loads 
%   -- raw checkerboard stimulus file 
%   -- spike train file
%   -- cell parameter file 
% and uses these to generate the stimulus and response matlab variables for
% the training and test sets. 80% of the data is used for training and 20%
% for test. The portion of the data that is used as a test set is shifted
% to obtain a jackknife estimates of the fitted models and their prediction
% errors.

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
addpath(genpath(workdir)) ;

datadir = 'RetinaData/' ;
        % the code assumes data-set files are saved in a subfolder of the 
        % parent directory called 'RetinaData'. These files include the
        % checkerboard stimulus, responses of 53 cells and a matlab file
        % which specifies the stimulus patch within the entire stimulus 
        % each cell was sensitive to and the lag: how many stimulus frames 
        % into the past are relevant for a cell's response. The cell
        % parameters were found by looking at the Spike-Triggered-Average
        % of the entire stimulus over a large number of frames.

stim_length = {'short2','short3','long'} ;

nJK = 5 ; % number of jack-knives 

for icell = 3:3
        % results for cell 3 are presented in paper. to generate stimulus and
        % response for other cells change to a vector including all the cell
        % numbers you are interested in (from 1 to 53).
    
    for iL = 1:3
        % there are three stimulus configurations we worked with:
        % long:   N_X = 10 x 10 = 100 pixels, N_T = 6 frames, N = 600
        % short3: N_X = 14 x 14 = 196 pixels, N_T = 3 frames, N = 588
        % short2: N_X = 10 x 10 = 100 pixels, N_T = 2 frames, N = 200
        %         short2 is not shown in the paper. it was used to illustrate
        %         over-fitting effects seen in the MNE models
        
        % loading the cell parameters:
        load([datadir 'RetinaCellParameters_' stim_length{iL} '.mat']) ;
        lag = lagshifts(icell) ;
        fsize = Nv(icell)^2 ;
        
        % loading the spike train:
        fid = fopen([datadir '/whitenoisec' int2str(icell) '.isk'], 'rt') ;
        R = textscan(fid, '%u\n') ;
        fclose(fid) ;
        R = double(R{1,1}) ;
        
        T = length(R) ; % total number of stimulus frames
        
        % loading the stimulus:
        fid = fopen([datadir 'whitenoise.raw'], 'rb');
        S = ReadFramev2(fid,T,Nx,Nv(icell),cx,x0(icell),y0(icell));
        % cutting out the patch the specific cell was sensitive to.
        % patch dimensions vary for different stimulus configurations
        fclose(fid);
        
        NX = length(S)/T ; % Spatial dimension of stimulus 
        N = NX*NT ;        % Stimulus dimension is equal to product of spatial dimension and temporal dimension
                           % temporal dimension (NT) is saved in and loaded with cell parameter file 
        S = reshape(S, [NX T])';
        S = S - repmat(mean(S), [T 1]);
        S = S./repmat(std(S), [T 1]);
        
        R = circshift(R,-lag);
        R = R(1:end-lag);
        T = length(R);
        S = S(1:end-lag,:);
        % shifting the response and stimulus to align them (and take response lag into account)   

        if NT>1
            % if the cell is thought to depend on the stimulus at multiple
            % times, this loop duplicates and expands the stimulus such 
            % that the redefined stimulus includes the spatial pattern presented 
            % at multiple lags (see Equations 12 and 13)
            T = T - (NT-1);
            S1 = zeros(T,N);
            for i=1:NT
                S1(:,NX*(i-1)+1:NX*i) = S(i:T+i-1,:) ;
            end
            S = S1;
            clear S1
            R = R(NT:length(R));
        end
        Sa = S ;
        Ra = R ;
        Ta = T ;
        
        % This loop splits the aligned stimulus and response variables into
        % training and test sets. Throughout the code:
        % S  - Stimulus used for training
        % St - Stimulus used for test
        % R  - Response used for training
        % Rt - Response used for test
        % T  - number of stimulus frames used for training
        % Tt - number of stimulus frames used for test
        
        for iJK = 1:nJK
            itest = ( round(Ta*(iJK-1)/nJK)+1 ) : 1 : round(Ta*iJK/nJK) ;
            ifit  = setdiff(1:1:Ta,itest) ;
            St = Sa(itest,:) ;
            Rt = Ra(itest) ;
            S  = Sa(ifit,:) ;
            R  = Ra(ifit) ;
            Tt = length(itest) ;
            T  = length(ifit) ;
            save([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'N','NX','NT','T','Tt','S','R','St','Rt') ;
        end
    end
end