% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the variables extracted from recordings limited to 
% times of whisking bouts and builds the stimulus and response matrices
% that will be used to fit and test LN models 

% The stimulus is the time history of the of the whisker position relative
% to the set-point (for definitions and code see
% script_201_WhiskingData_ExtractWhisking.m).

% The stimulus has two parameters: 
% 1. the dimensionality of the stimulus space: how many time points into the past the neuron is sensitive to. 
% 2. the time resolution of the stimulus space: how long are intervals between the stimulus history points. 
% Choosing values for these parameters can have a strong influence on the
% outcome of the analysis, and the considerations one should make when
% choosing them is explained in the manuscript. 

clear ;
workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ; % Cell indices (used in file names) 
sr = 1000 ;     % (Hz) sampling rate
N = 150 ;       % stimulus dimensionality 
ds = 2 ;        % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
                % delta T for down sampled stimulus is ds/sr (here, 2ms)
                % the range of the stimulus is p/(ds*sr) (here, 300ms)

nJK = 5 ;       % number of jackknives
for i = 1:7     % loop over cells
    load([datadir 'VPM_cell_' Names{i} '_iswhisking.mat']) ; % loading data (filtered whisker position and spikes extracted for all whisking bouts) 
    
    for iJK = 1:nJK
        T  = 1 ;  % number of time points in training set  
        Tt = 1 ; % number of time points in test set 
    
        S  = zeros(70000,N) ;   % array holding stimulus for training
        St = zeros(70000,N) ;   % array holding stimulus for test 
        Sph  = zeros(70000,1) ; % array holding phase stimulus for training
        Spht = zeros(70000,1) ; % array holding phase stimulus for test 
                                % (phase information was used only to fit a
                                % tuning curve model)
        Sam  = zeros(70000,1) ; % array holding amplitude stimulus for training
        Samt = zeros(70000,1) ; % array holding amplitude stimulus for test
                                % (amplitude stimulus was not used directly
                                % to fit any model, but it could in
                                % principle be used, for example by
                                % defining a 2 dimensional tuning curve:
                                % amplitude x phase
                                
        R  = zeros(70000,1) ; % array holding response for training
        Rt = zeros(70000,1) ; % array holding response for test
    
        ktst = (round(k*(iJK-1)/nJK)+1) : 1 : round(k*iJK/nJK) ; % whisking bouts are split into training and test sets. Since whisking bouts are relatively short, they are not split in the middle
                                                                 % (training and test sets contain only full whisking bouts)
        ktrn = setdiff(1:k,ktst) ;
        for ik = ktrn                                        % looping over training set whisking bouts
            R_k = round((isw_ispk{ik}-isw_isam{ik}(1))/ds) ; % downsampled spike time indices for current whisking bout 
            R_k(R_k<N) = [] ;                                % spikes with indices less than p (stimulus dimensionality) are thrown (they depend on whisker position at time where there was no whisking)  
            R(T+R_k) = 1 ;                                   % response array set to 1 at spike times 
    
            S_k = isw_pos{ik}-isw_spt{ik} ;                  % stimulus defined to be whisker position relative to set-point
            S_k = S_k(1:ds:end) ;                            % downsampled stimulus for current whisking bout 
        
            Sph_k = isw_phs{ik} ;                            % phase stimulus (extracted using hilbert transform)
            Sph_k = Sph_k(1:ds:end) ;                        % downsampled phase stimulus for current whisking bout 
        
            Sam_k = isw_amp{ik} ;                            % amplitude stimulus (extracted using hilbert transform)
            Sam_k = Sam_k(1:ds:end) ;                        % downsampled amplitude stimulus for current whisking bout 
        
            Sph(T+(N:length(Sph_k))) = Sph_k(N:length(Sph_k)) ; 
            Sam(T+(N:length(Sam_k))) = Sam_k(N:length(Sam_k)) ; 

            for it = N:length(S_k) ;                         % looping over whisking bout duration: 
                S(T,:) = S_k((it-N)+(1:N)) ;                 %    at each time point, spiking is assumed to depend
                                                             %    on the stimulus history N time points into the past
                                                             %    so starting from time N*dt into the bout the stimulus is taken
                                                             %    to be the last N stimulus points 
            
                T = T+1 ;
            end
        end
        % after building the stimulus and response using all the training set 
        % whisking bouts, the variables S, R at times greater than T are thrown away.
        S = S(1:T-1,:) ;    
        Sph = Sph(1:T-1) ;
        Sam = Sam(1:T-1) ;
        R = R(1:T-1) ;   
    
        % Exact same steps are repeated to build the stimulus and response for
        % the training set (St, Rt, Tt - respectively)
        for ik = ktst
            R_k = round((isw_ispk{ik}-isw_isam{ik}(1))/ds) ;
            R_k(R_k<N) = [] ;
            Rt(Tt+R_k) = 1 ;

            S_k = isw_pos{ik}-isw_spt{ik} ;
            S_k = S_k(1:ds:end) ;
        
            Sph_k = isw_phs{ik} ;
            Sph_k = Sph_k(1:ds:end) ; 
        
            Sam_k = isw_amp{ik} ;
            Sam_k = Sam_k(1:ds:end) ; 
        
            Spht(Tt+(N:length(Sph_k))) = Sph_k(N:length(Sph_k)) ; 
            Samt(Tt+(N:length(Sam_k))) = Sam_k(N:length(Sam_k)) ; 

            for it = N:length(S_k) ;
                St(Tt,:) = S_k((it-N)+(1:N)) ;
                Tt = Tt+1 ;
            end
        end
        St   = St(1:Tt-1,:) ;
        Spht = Spht(1:Tt-1) ;
        Samt = Samt(1:Tt-1) ;
        Rt   = Rt(1:Tt-1) ;
        
        T  = T  - 1 ;
        Tt = Tt - 1 ;
        save([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat'],'N','T','Tt','S','Sph','Sam','R','St','Spht','Samt','Rt','sr','ds') ;
    end
end