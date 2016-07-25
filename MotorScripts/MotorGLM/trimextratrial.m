function [MSTM, SPNDS, SPNDS2] = trimextratrial(MSTM, SPNDS, processed, SPNDS2)
	%Chop out parts from MSTM, SPNDS, SPNDS2 that are outside of trial movement times
	%(for GLM training purposes)
	nB = size(MSTM,1);
	pt = [];
	dt = processed.dt;
	for idx = 1:size(processed.trialstartend,1)
		tstart = processed.trialstartend(idx,1);
		tend = processed.trialstartend(idx,2);
		if tend < nB
			pt(idx,1:2) = [tstart, tend];
		else
			break
		end
	end
	%Chop out spikes past end of stim
	%sp = SPNDS;
	SPNDS(SPNDS*dt>nB) = [];
	%Make list of intervals to be excluded from data
	trialendstart = [[1; pt(:,2)], [pt(:,1); nB]];
	%Chop out spikes that are outside of a trial
	for idx = 1:size(trialendstart,1)
		trimstart = trialendstart(idx,1);
		trimend = trialendstart(idx,2);
		%display([num2str(idx), ': ', num2str(trialendstart(end,2)), ' size MSTM: ', num2str(size(MSTM,1))]) 
		MSTM(trimstart:trimend,:) = [];
		Dt = (trimend-trimstart+1)/processed.dt;
		SPNDS((SPNDS*dt>trimstart) & (SPNDS*dt<trimend)) = [];
		SPNDS(SPNDS*dt>trimend) = SPNDS(SPNDS*dt>trimend)-Dt;
		trialendstart = trialendstart-(trimend-trimstart+1);
		%pause
	end

	if nargin == 4
		nC = length(SPNDS2);
		for ii = 1:nC
			%Chop out spikes past end of stim
			SPNDS2{ii}(SPNDS2{ii}*dt>nB) = [];
			%Make list of intervals to be excluded from data
			trialendstart = [[1; pt(:,2)], [pt(:,1); nB]];
			%Chop out spikes that are outside of a trial
			for idx = 1:size(trialendstart,1)
				trimstart = trialendstart(idx,1);
				trimend = trialendstart(idx,2);
				Dt = (trimend-trimstart+1)/processed.dt;
				SPNDS2{ii}((SPNDS2{ii}*dt>trimstart) & (SPNDS2{ii}*dt<trimend)) = [];
				SPNDS2{ii}(SPNDS2{ii}*dt>trimend) = SPNDS2{ii}(SPNDS2{ii}*dt>trimend)-Dt;
				trialendstart = trialendstart-(trimend-trimstart+1);
			end
		end
	else
		SPNDS2 = {};
	end
end

