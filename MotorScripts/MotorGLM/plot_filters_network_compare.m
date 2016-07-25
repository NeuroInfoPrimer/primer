function plot_filters_network_compare(models, models_uncoupled, processed, fn_out)
	%Plot filters of a fitted GLM model, along with other fit statistics
	%     
	%Input:
	%	models = data structure output by function in ./fitting (containing fitted coefficients)
	%	fn_out = base filename to write plots to for each unit

	global RefreshRate;
	nM = length(models);
	nU = nM; %number of units
	names = {'spike history', 'curs x', 'curs y', 'curs z', 'grip'};
	nK = nU+5;

	%Prepare subplots
	plotheight=nU+5;
	plotwidth=nK+5;
	subplotsx=nK+1;
	subplotsy=nU+1;   
	leftedge=1.2;
	rightedge=0.4;   
	topedge=1;
	bottomedge=1.5;
	spacex=0.2;
	spacey=0.4;
	fontsize=5;    
	sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
	sub_pos = fliplr(sub_pos);

	%setting the Matlab figure
	f=figure;
	clf(f);
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperSize', [plotwidth plotheight]);
	set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

	ymincurs = 0;
	ymaxcurs = 0;
	ymingrip = 0;
	ymaxgrip = 0;
	yminsp = 0;
	ymaxsp = 0;
	for i = 1:nU
		model = models{i};
		stimfilt = model.k;
		const = model.dc;
		sphist = model.ihbas*model.ih;
		couples = model.ihbas*model.ih2;
		spfilters = couples;
		ymincurs = min(ymincurs, min(min(stimfilt(:,1:3))));
		ymaxcurs = max(ymaxcurs, max(max(stimfilt(:,1:3))));
		ymingrip = min(ymingrip, min(min(stimfilt(:,4))));
		ymaxgrip = max(ymaxgrip, max(max(stimfilt(:,4))));
		yminsp = min(yminsp, min(min(spfilters)));
		ymaxsp = max(ymaxsp, max(max(spfilters)));
	end

	for i = 1:nU
		model = models{i};
		uncoupled = models_uncoupled{i};
		idx = 1;
		stimfilt = model.k;
		const = model.dc;
		sphist = model.ihbas*model.ih;
		couples = model.ihbas*model.ih2;
		sphist_uncpl = uncoupled.ihbas*uncoupled.ih;
		stimfilt_uncpl = uncoupled.k;
		spfilters = couples;

		%Spike history filter
		ax=axes('position',sub_pos{1,i},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
		name = processed.unitnames{i};
		filt = sphist;
		filt = flipud(filt);
		filt_un = sphist_uncpl;
		filt_un = flipud(filt_un);
		dt_filt = model.dt*RefreshRate;
		tt = (0:length(filt)-1)*dt_filt;
		tt = (tt-max(tt))/RefreshRate/model.dt;
		ttu = (0:length(filt_un)-1)*dt_filt;
		ttu = (ttu-max(ttu))/RefreshRate/model.dt;
		plot(tt, filt, 'k', ttu, filt_un, 'k--');
		if length(tt) > 1
			xlim([min(tt) max(tt)]);
		end
		if i == 1
			title(name);
		end
		if i == nM
			xlabel('time (ms)');
		end

		%Coupling filters
		for j = 1:(nU)
			ax=axes('position',sub_pos{j+1,i},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			if j == i
				plot(0,0)
				axis off
				if j == 1
					name = processed.unitnames{j};
					title(['Unit: ' num2str(i) ' ' name]);
				end
			else
				if j > i
					jp = j-1;
				else
					jp = j;
				end
				name = processed.unitnames{j};
				filt = spfilters(:,jp);
				filt = flipud(filt);
				dt_filt = model.dt*RefreshRate;
				tt = (0:length(filt)-1)*dt_filt;
				tt = (tt-max(tt))/RefreshRate/model.dt;
				plot(tt, filt);
				ylim([yminsp, ymaxsp])
				if length(tt) > 1
					xlim([min(tt) max(tt)]);
				end
				if i == 1
					title(['Unit: ' num2str(j) ' ' name]);
				end
				if i == nM
					xlabel('time (ms)');
				end
			end
		end

		%Stim filters
		for j = nU+(1:4)
			ax=axes('position',sub_pos{j+1,i},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
			name = names{j-nU+1};
			filt = stimfilt(:,j-nU);
			filt_uncpl = stimfilt_uncpl(:,j-nU);
			dt_filt = model.dt*RefreshRate;
			tt = (0:length(filt)-1)*dt_filt;
			tt = (tt-max(tt)/2);
			ttu = (0:length(filt_uncpl)-1)*dt_filt;
			ttu = (ttu-max(ttu)/2);
			hold on
			if j == nU+4
				plot(tt, filt, 'r', ttu, filt_uncpl, 'r--');
				ylim([ymingrip, ymaxgrip])		
			else
				plot(tt, filt, 'Color', [0 .5 0]);
				hold on 
				plot(ttu, filt_uncpl, '--', 'Color', [0 .5 0])
				ylim([ymincurs, ymaxcurs])		
			end
			if length(tt) > 1
				xlim([min(tt) max(tt)]);
			end
			if i == 1
				title([name]);
			end
			if i == nM
				xlabel('time (ms)');
			end
		end

		if ischar(processed.unitnames)
			name = processed.unitnames;
		else
			name = processed.unitnames{i};
		end
		%Plot information about each subplot
		ax=axes('position',sub_pos{nK+1,i},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
		str1(1) = {['Unit no.: ' num2str(i) ' Unit: ' name]};
		%str1(2) = {['Deviance: ' num2str(dev)]};
		%str1(3) = {['Degrees of freedom: ' num2str(stats.dfe)]};
		str1(2) = {['Binsize: ' num2str(processed.binsize)]};
		str1(3) = {['Seconds of training: ' num2str(size(processed.cursor,1)*processed.binsize)]};
		str1(4) = {['Number of spikes: ' num2str(sum(processed.spiketrain(:,i)))]};
		text(0.1,0.8,str1, 'FontSize', 5, 'Interpreter', 'none')
		axis off
	end
	saveplot(gcf, fn_out, 'eps', [plotwidth, plotheight]);	
end