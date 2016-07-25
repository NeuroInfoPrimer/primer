function [processed, processed_withheld] = preprocess_crossval(datafile, binsize, dt, frames, fold, nfolds, stddev)
	%Preprocess both spike data and stim data
	%
	%Usage:
	%		[processed, processed_withheld] = preprocess_crossval(datafile, binsize, dt, frames, stddev)
	%
	%Input:
	%		datafile = .mat file with data
	%		binsize = size in seconds of each time bin for stim
	%		dt = relative size of bins spike times are provided at		
	%		frames = number of previous stim frames to include
	%		fold = fold number of cross-validation
	%		nfolds = total number of cross-val folds
	%		stddev = (optional, default = 0) whether to normalize stim
	%	
	%Output:
	%		processed is a structure containing the following fields:
	%			spikes = [nB x nU] array with spikes from all channels binned 
	%				according to binsize. 
	%				nB = no. bins, nU = no. units.
	%			cursor = [nB x 3] array of stimuli presented to 'unitidx'th cell
	%			grip = [nB x 1] array of grip force
	%			binsize = binsize used
	%			unitnames = names of units loaded into spikes
	%			unitidx = index of unit to fit GLM to
	%			stacked = stacked stim for computing STA
	%		processed_withheld is a structure containing 20% of data withheld for
	%			testing
	%
	%Test code:
	%	datafile = '../MotorData/mabel.mat';
	%	binsize = 1/100;
	%	dt = 0.1;
	%	frames = 4;
	%	fold = 1;
	%	nfolds = 5;
	%	processed = preprocess_crossval(datafile, binsize, dt, frames, fold, nfolds);

	%Event legend
	if (nargin < 7) stddev = 0; end
	TRIALEND = 15;
	TARGETAPPEARED = 110;
	TARGETREACHED = 111;
	ds = 0.01;

    load(datafile);
    processed.binsize = binsize;
    processed.cursor = [Cursor_X, Cursor_Y, Cursor_Z];
    processed.grip = Grip_force;
    processed.stim = [processed.cursor, processed.grip];
    %Normalize
    if stddev == 1
    	for idx = 1:size(processed.stim,2)
    		idx
	    	processed.stim(:,idx) = (processed.stim(:,idx)-mean(processed.stim(:,idx)))/std(processed.stim(:,idx));
    	end
	end
    processed.stim = resample(processed.stim, 1000*ds, 1000*binsize);
    processed.dt = dt;
    unitnames = who('CSPIK*');
    processed.unitnames = unitnames;
    processed.frames = frames;
    nU = length(unitnames);
    nB = size(processed.stim,1);
    nS = size(processed.stim,2);

    processed.spikes = {};
    processed.spiketrain = zeros(nB, nU);
    for idx = 1:nU
    	%Get data for idxth spike train
    	eval(['spikes = ' unitnames{idx} ';']);
    	processed.spikes{idx} = spikes*dt;
    	bins = ceil(spikes*dt);
    	for b = bins
	    	processed.spiketrain(b,idx) = processed.spiketrain(b,idx)+1;
		end
	end

    %Stack stim to include past and future frames relative to spike time
    processed.stacked = zeros(nB, nS*(2*frames+1));
    offsets = -frames:frames;
	nx = nS;
	nt = length(offsets);
	stim = vertcat(zeros(frames, nx), processed.stim, zeros(frames, nx));
	nB = size(stim,1);
	for i = 1:nx
	    for j = 1:nt
	        offset = offsets(j);
	        jj = (frames+offset+1):(nB-frames+offset);
	        processed.stacked(:,(i)*nt-j+1) = stim(jj, i);
	    end
	end

	%Here trial start is determined by Est_Movement_Initiation variable...
	%Note which bins are inside a trial
	trialendidx = find(Events_Data(2,:)==TRIALEND);
	trialstartbins = ceil(Est_Movement_Initiation*dt);
	trialendbins = ceil(Events_Data(1,trialendidx)*dt);
	%Only take the trial ends that are immediately after a movement start
	tendbins = [];
	for idx = 1:length(trialstartbins)
		startbin = trialstartbins(idx);
		%Find next trial end
		tend_idx = find(trialendbins > startbin, 1);
		tendbin = trialendbins(tend_idx);
		if ~isempty(tendbin)
			%If this is before the _next_ movement initiation, then take the pair
			nextstart_idx = find(trialstartbins > startbin,1);
			nextstart = trialstartbins(nextstart_idx);
			if isempty(nextstart)
				tendbins = [tendbins, tendbin];
			else 
				if nextstart > tendbin
					tendbins = [tendbins, tendbin];
				end
			end
		end
	end
	trialendbins = tendbins;

	trialchanges = zeros(size(processed.grip));
	trialchanges(trialstartbins) = 1;
	trialchanges(trialendbins) = -1;
	processed.intrial = cumsum(trialchanges);
	%Start outside of of a trial
	processed.intrial(1) = 0;
	processed.trialstartend = [trialstartbins', trialendbins'];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Split into train and test sets%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Split into 20% test and 80% training set
	nB = size(processed.spiketrain,1);
	blocksz = floor(nB/nfolds);
	preidx = (fold-1)*blocksz;
	postidx = fold*blocksz;
	trainingidx = [1:preidx, (postidx+1):(blocksz*nfolds)];
	testingidx = (preidx+1):postidx;

	processed_withheld = processed;
	processed.cursor = processed.cursor(trainingidx,:);
	processed.grip = processed.grip(trainingidx,:);
	processed.stim = processed.stim(trainingidx,:);
	processed.spiketrain = processed.spiketrain(trainingidx,:);
	processed.stacked = processed.stacked(trainingidx,:);
	processed.intrial = processed.intrial(trainingidx);

	processed_withheld.cursor = processed_withheld.cursor(testingidx,:);
	processed_withheld.grip = processed_withheld.grip(testingidx,:);
	processed_withheld.stim = processed_withheld.stim(testingidx,:);
	processed_withheld.spiketrain = processed_withheld.spiketrain(testingidx,:);
	processed_withheld.stacked = processed_withheld.stacked(testingidx,:);
	processed_withheld.intrial = processed_withheld.intrial(testingidx);

	%Update spike time matrices
	for idx = 1:nU
	    sp = processed.spikes{idx};
    	sptrain = sp(sp<preidx);
	    sptest = sp(sp>preidx & sp < postidx);
	    sptrainpost = sp(sp > postidx)-postidx+preidx;
    	processed.spikes{idx} = [sptrain, sptrainpost];
	    processed_withheld.spikes{idx} = sptest-preidx;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Update trial info matrices
	trialstartend = [];
	processed_withheld.trialstartend = [];
	%For all trials
	for idx = 1:size(processed.trialstartend,1)
		%For trials during the pretraining
		if processed.trialstartend(idx,2) < preidx 
			tstart = processed.trialstartend(idx,1);
			tend = processed.trialstartend(idx,2);
			trialstartend = [trialstartend; tstart, tend];
		end
		%For trials during the testing period
		if processed.trialstartend(idx,1) > preidx & processed.trialstartend(idx,2) < postidx
			tstart = processed.trialstartend(idx,1)-preidx;
			tend = processed.trialstartend(idx,2)-preidx;
			processed_withheld.trialstartend = [processed_withheld.trialstartend; tstart, tend];
		end
		%For trials during the posttraining period
		if processed.trialstartend(idx,1) > postidx
			tstart = processed.trialstartend(idx,1)-postidx+preidx;
			tend = processed.trialstartend(idx,2)-postidx+preidx;
			trialstartend = [trialstartend; tstart, tend];
		end
	end
	processed.trialstartend = trialstartend;
end


