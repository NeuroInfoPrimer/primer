wd = '../MotorData/';

icell = 6;
ttimes = [655.64, 675.64;
			635.64, 655.64;
			615.64, 635.64;
			595.64, 615.64;
			575.64, 595.64;
			555.64, 575.64;
			535.64, 555.64;
			515.64, 535.64;
			495.64, 515.64;
			475.64, 495.64;
			675.64, 695.64;
			695.64, 715.64;
			715.64, 735.64;
			735.64, 755.64;
			755.64, 775.64;
			775.64, 795.64];
clf
for ii = 1:size(ttimes,1)
	datafile = 'mabel.mat';
	binsize = 0.01;
	dt = .1;
	frames = 1;
	load([wd datafile]);
	processed = preprocess([wd datafile], binsize, dt, frames);
	spikes = processed.spiketrain(:,icell);
	%Legend:
	TRIALSTART = 10;
	TRIALEND = 15;
	ORIGINAPPEARS = 100;
	ORIGINREACHED = 101;
	TARGETAPPEARS = 110;
	TARGETREACHED = 111;
	GOSIGNAL = 103;
	RELEASEGRIP = 104;
	GRIPPRESSED = 368;
	GRIPRELEASED = 369;
	
	figure
	%Whole recording lasts 6000s...
	tstart = ttimes(ii,1); tend = ttimes(ii,2);
	tidx = (tstart*100):(tend*100);
	tt = tidx/100;
	
	%Plot grip force
	subplot(5,1,1)
	plot(tt, Cursor_X(tidx));
	xlim([tstart, tend])
	ylabel('cursor X (cm)')
	subplot(5,1,2)
	plot(tt, Cursor_Y(tidx));
	xlim([tstart, tend])
	ylabel('cursor Y (cm)')
	subplot(5,1,3)
	plot(tt, Cursor_Z(tidx));
	xlim([tstart, tend])
	ylabel('cursor Z (cm)')
	subplot(5,1,4)
	hold on
	h(1) = plot(tt, Grip_force(tidx), 'k');
	xlim([tstart, tend]);

	%h(1) = plot(tt, Cursor_X(tidx), '--', tt, Cursor_Y(tidx), '--', tt, Cursor_Z(tidx), '--');
	colors = colormap;
	
	%Plot Events
	%Only plots events within range
	idxstart = find(Events_Data(1,:)>=(tstart*1000), 1);
	idxend = find(Events_Data(1,:)>=(tend*1000), 1);
	evttimes = Events_Data(1,idxstart:idxend)/1000;
	evts = Events_Data(2,idxstart:idxend);
	h(2) = plot(evttimes, ones(size(evttimes))*.01, 'x');
	xlim([min(tt), max(tt)]);
	
	%Add trial start/end lines
	intrial = 0;
	inorigin = 0;
	intarget = 0;
	ingo = 0;
	ingrip = 0;
	trialstart = 0;
	originstart = 0;
	targetstart = 0;
	gostart = 0;
	gripstart = 0;
	
	for idx = 1:length(evts)
		evt = evts(idx);
		evttime = evttimes(idx);
		if evt == TRIALSTART;
			intrial = 1;
			trialstart = evttime;
		elseif evt == TRIALEND;
			if ~intrial
				trialstart = tstart;
			end
			%Draw line between here and trialstart
			h(3) = plot([trialstart, evttime], [-0.02, -0.02], 'r', 'LineWidth', 2);
			%line([trialstart, evttime], [-0.01, -0.01], 'LineWidth', 2)
			intrial = 0;
		elseif evt == ORIGINAPPEARS;
			inorigin = 1;
			originstart = evttime;
		elseif evt == ORIGINREACHED;
			if ~inorigin
				originstart = tstart;
			end
			%Draw line between here and trialstart
			h(4) = plot([originstart, evttime], [-0.06, -0.06], 'b', 'LineWidth', 2);
			%line([originstart, evttime], [-0.03, -0.03], 'LineWidth', 2)
			inorigin = 0;
		elseif evt == TARGETAPPEARS;
			intarget = 1;
			targetstart = evttime;
		elseif evt == TARGETREACHED;
			if ~intarget
				targetstart = tstart;
			end
			%Draw line between here and trialstart
			h(5) = plot([targetstart, evttime], [-0.1, -0.1], 'g', 'LineWidth', 2);
			%line([targetstart, evttime], [-0.05, -0.05], 'LineWidth', 2)
			intarget = 0;
		elseif evt == GOSIGNAL;
			ingo = 1;
			gostart = evttime;
		elseif evt == RELEASEGRIP;
			if ~ingo
				gostart = tstart;
			end
			%Draw line between here and trialstart
			h(6) = plot([gostart, evttime], [-0.14, -0.14], 'k', 'LineWidth', 2);
			%line([gostart, evttime], [-0.07, -0.07], 'LineWidth', 2)
			ingo = 0;
		elseif evt == GRIPPRESSED;
			ingrip = 1;
			gripstart = evttime;
		elseif evt == GRIPRELEASED;
			if ~ingrip
				gripstart = tstart;
			end
			%Draw line between here and trialstart
			h(7) = plot([gripstart, evttime], [-0.18, -0.18], 'y', 'LineWidth', 2);
			%line([gripstart, evttime], [-0.09, -0.09], 'LineWidth', 2)
			ingrip = 0;
		end
	end 
	
	%Close off any events unfinished:
	if intrial
		h(8) = plot([trialstart, tend], [-0.02, -0.02], 'r', 'LineWidth', 2);
		%line([trialstart, tend], [-0.01, -0.01], 'LineWidth', 2)
	end
	if inorigin
		h(9) = plot([originstart, tend], [-0.06, -0.06], 'b', 'LineWidth', 2);
		%line([originstart, tend], [-0.03, -0.03], 'LineWidth', 2)
	end
	if intarget
		h(10) = plot([targetstart, tend], [-0.1, -0.1], 'g', 'LineWidth', 2);
		%line([targetstart, tend], [-0.05, -0.05], 'LineWidth', 2)
	end
	if ingo
		h(11) = plot([gostart, tend], [-0.14, -0.14], 'k', 'LineWidth', 2);
		%line([gostart, tend], [-0.07, -0.07], 'LineWidth', 2)
	end
	if ingrip
		h(12) = plot([gripstart, tend], [-0.18, -0.18], 'y', 'LineWidth', 2);
		%line([gripstart, tend], [-0.09, -0.09], 'LineWidth', 2)
	end
	ylabel('grip force')
	hl = legend([h([1 3:7])], 'Grip force', 'Within trial', 'Origin seek', 'Target seek', 'Go signal', 'Grip pressed');
	
	subplot(5,1,5)
	%Smooth spikes
	RefreshRate = 100;
	sigma_fr = .1;
	sigma_fr = sigma_fr*RefreshRate;
	sz = sigma_fr*3*2;
	x = linspace(-sz/2, sz/2, sz);
	gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
	gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
	gftruesp = conv(spikes, gaussFilter_fr, 'same');
	
	%Plot spikes
	plot(tt, gftruesp(tidx)*RefreshRate);
	xlim([tstart, tend])
	xlabel('time (s)')
	ylabel('spikes/s')
	ylim([0 50])

	%Add ticks for spikes:
	hold on 
	truesp = processed.spikes{icell}*binsize;
	truesp = truesp(truesp < tend & truesp > tstart);
	for idx = 1:length(truesp)
		plot([truesp(idx) truesp(idx)], [0 5], 'k')
	end
	saveplot(gcf, [wd '/trialsample_cell_' num2str(icell) '_' num2str(ii) '.eps'], 'eps', [10 10])
end

%Plot autocorrelation too...
%Plot a bunch of preprocessing diagnostics
figure
%Compute auto- and cross-correlation in torque and example firing rate
samplerate = 100;
binsize = 1/samplerate;
maxlag = 50;

autotorqueX = xcov(Cursor_X(:)-mean(Cursor_X(:)),samplerate*maxlag);%, 'coeff');
autotorqueY = xcov(Cursor_Y(:)-mean(Cursor_Y(:)),samplerate*maxlag);%, 'coeff');
autotorqueZ = xcov(Cursor_Z(:)-mean(Cursor_Z(:)),samplerate*maxlag);%, 'coeff');
autotorqueGrip = xcov(Grip_force(:)-mean(Grip_force(:)),samplerate*maxlag);%, 'coeff');

tt = -maxlag:binsize:maxlag;
subplot(2,2,1)
plot(tt, autotorqueX/autotorqueX(5001));
title('Cursor X');		
xlabel('sec')
subplot(2,2,2)
plot(tt, autotorqueY/autotorqueY(5001));
title('Cursor Y');
xlabel('sec')
subplot(2,2,3)
plot(tt, autotorqueZ/autotorqueZ(5001));
title('Cursor Z')
xlabel('sec')
subplot(2,2,4)
plot(tt, autotorqueGrip/autotorqueGrip(5001));
title('Grip force');		
xlabel('sec')
saveplot(gcf, [wd '/stim_autocorrelations.eps'], 'eps', [6 10]);
