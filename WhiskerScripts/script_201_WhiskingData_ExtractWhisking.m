% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the whisking data and extracts the parts of the 
% recording where whisking was detected  (see detailed explanation on data 
% format in WhiskerData_Readme.txt). 

% It was previously shown that the whisker location holds predictive power 
% over spikes in VPM thalamus only during active whisking bouts, and not 
% during times when the rat is not whisking. Hence, the first step of 
% constructing a feed-forward LN model for the cell is to extract the 
% whisker position and spikes that occured during the active whisking bouts 

clear ;
workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
addpath(genpath(workdir)) ;
datadir = 'WhiskerData/' ;
Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ; % Cell indices (used in file names) 

sr = 1000 ; % (Hz) sampling rate

% Filtering - The whisker position is described as an angle in body
% centered coordinates. The mid point of the whisker oscillation (also 
% referred to here as set point) varies considerably from bout to bout but
% has relatively little predictive power for spiking of neurons in VPM
% thalamus. Therefore the whisk position is band-pass filtered to remove 
% the set point and high frequency noise, leaving only the relative whisk 
% position that can be decomposed into an amplitude and a phase. 


bp = [1 2] ; % lower and upper limits of band pass parameters of filter 
setpt_func = inline( '(max(x) + min (x)) / 2') ; % function describing the setpoint location 
amp_func  = @range ;                             % function describing the magnitude of amplitude 
thr = 10 ;  % (degrees) threshold on whisker postion amplitude above which that part of a recording is taken to be a whisking bout 
a = 0.005 ; % parameter for amplitude smoothing filter such that whisking bout 'cutouts' are not too short

lag = 750 ;                 % number of lags used to compute whisker autocorrelation
tau = (-lag:1:lag)/sr ;     % autocorrelation temporal lag vector

% for each whisking bout recorded we will have a cell array that contains
% the following variables
isw_isam = cell(1,100) ; % sample numbers                         
isw_ispk = cell(1,100) ; % sample numbers during which a spike was recorded
isw_itop = cell(1,100) ; % sample numbers during the whisker position has a peak 
isw_pos  = cell(1,100) ; % whisker position
isw_phs  = cell(1,100) ; % whisker phase 
isw_amp  = cell(1,100) ; % whisker amplitude  
isw_spt  = cell(1,100) ; % whisker set-point 


for i = 1:7 % loop over neurons 
    load([datadir 'VPM_cell_' Names{i} '.mat']) ; % load data for a particular neuron 
    nrec = length(whiskparams) ;                  % number of recording intervals for specific neuron  
    k = 0 ;                                       % whisking bout index 
    
    % cell arrays that will hold relevant information for each recording
    % interval: 
    time       = cell(1,nrec) ;                        
    phase      = cell(1,nrec) ;
    amplitude  = cell(1,nrec) ;
    spikes     = cell(1,nrec) ;
    setpoint   = cell(1,nrec) ;
    position   = cell(1,nrec) ;
    iswhisking = cell(1,nrec) ;   % a vector of indices pointing to times where the animal was whisking

    ACisw = zeros(2*lag+1,100) ;  % autocorrelation of whisker position during whisking 
    ACall = zeros(2*lag+1,nrec) ; % autocorrelation of whisker position during all times  

    for j = 1:nrec % loop over recording intervals
        wp = whiskparams(j) ;                             
        pos = wp.position' ;                            % angle during recording
        sam = wp.samples' ;                             % sample number
        spk = wp.spikes' ;                              % samples during which spike was recorded 
        
        phs = phase_from_hilbert(pos,250,bp) ;          % extracting phase from whisker position using hilbert transform
        
        [amp,itop] = get_slow_var(pos,phs,amp_func) ;   % extracting amplitude, indices of maxima and set point from whisker position using hilbert transform
        spt = get_slow_var(pos,phs,setpt_func) ;
        
        ACall(:,j) = xcorr(amp.*cos(phs)/2,lag,'unbiased')' ;

        ampfilt = filtfilt(a, [1 a-1],amp) ;            % filtered amplitude variable 
        iisw = find(heaviside(ampfilt-thr)) ;           % indices of times where the animal was whisking  
        isam = 1+(sam-sam(1))/40 ;                      % conversion of samples to vector indices 
        ispk = intersect(round((spk-sam(1))/40),isam) ; % indices during which spike was recorded
        
        iisw1 = zeros(size(isam)) ;
        iisw1(iisw) = 1 ;                               % converting 'is whisking' indices to a vector of 0,1's: 1 if during that time the rat is whisking and 0 otherwise  
        
        isw_con = bwconncomp(iisw1) ;                   % connected components of whisking give whisking bouts
        
        for ic = 1:isw_con.NumObjects 
            % looping over connected components of the vector iisw1 (each connected component is a whisking bout)
            % the variables in each whisking bout (indexed by k) are put into the appropriate cell array 
            k = k+1 ;
            iisw_temp = isw_con.PixelIdxList{ic} ;        % temporary vector of indices for the current bout
            isw_isam{k} = isam(iisw_temp) ;
            isw_ispk{k} = intersect(ispk,iisw_temp) ;     % spike times for current bout are intersection of all spike times in recording interval and times of current whisking bout
            isw_itop{k} = intersect(itop,iisw_temp) ;     % maxima times for current bout are intersection of all maxima in recording interval and times of current whisking bout 
            isw_pos{k} = pos(iisw_temp) ;
            isw_phs{k} = phs(iisw_temp) ;
            isw_amp{k} = amp(iisw_temp) ;
            isw_spt{k} = spt(iisw_temp) ;
            ACisw(:,k) = xcorr(isw_amp{k}.*cos(isw_phs{k})/2,lag,'unbiased')' ;
        end
        % variables saved for entire recording (not only during whisking)
        time{j}       = (1:length(sam))/sr ;                 
        phase{j}      = phs  ;
        amplitude{j}  = amp  ;
        spikes{j}     = ispk ;
        setpoint{j}   = spt  ;
        position{j}   = pos  ;
        iswhisking{j} = iisw ;
    end
    % each cell array was initialized to have 100 cells. here we are keeping only the ones containing the variables for the k whisking bouts 
    isw_isam = isw_isam(1:k) ;
    isw_ispk = isw_ispk(1:k) ; 
    isw_pos = isw_pos(1:k) ;         
    isw_phs = isw_phs(1:k) ;
    isw_amp = isw_amp(1:k) ;         
    isw_spt = isw_spt(1:k) ;         
    ACisw   = ACisw(:,1:k) ;    

    save([datadir 'VPM_cell_' Names{i} '_iswhisking.mat'],'k','isw_isam','isw_ispk','isw_itop','isw_pos','isw_amp','isw_phs','isw_spt') ;
    save([datadir 'VPM_cell_' Names{i} '_hilbert.mat'],'time','phase','amplitude','spikes','setpoint','position','nrec','iswhisking') ;
    save([datadir 'VPM_cell_' Names{i} '_autocorrelation.mat'],'tau','ACisw','ACall') ;
end