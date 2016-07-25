function jackknifecoherence(wd, fn_in, fn_out)

	load([wd fn_in])

	h = figure;
	g = figure;
	coh_05 = zeros(nU,1);
	coh_cpl_05 = zeros(nU,1);
	coh_SE_05 = zeros(nU,1);
	coh_SE_cpl_05 = zeros(nU,1);

	maxt = 20000;
	for icell = 1:nU
	    truesp = {};
	    simsp = {};
	    simsp_cpl = {};
	    runs = size(proc_withheld.trialstartend,1);
	    K = [];
	    N_sample_max = 0;
	    uncoupled = load([wd '/GLM_cell_simulation_' num2str(icell) '.mat']);
	    nspks = 0;
	    for i = 1:runs
	        ts = proc_withheld.trialstartend(i, 1);
	        te = proc_withheld.trialstartend(i, 2);
	        if ts > maxt | te > maxt
	        	break
	        end
	        truesp{i} = proc_withheld.spiketrain(1:maxt,icell);
	        truesp{i} = truesp{i}(ts:te);
	        simsp_cpl{i} = Rt_glm{icell}(:);
	        simsp_cpl{i} = simsp_cpl{i}(ts:te);
	        simsp{i} = uncoupled.Rt_glm(:);
	        simsp{i} = simsp{i}(ts:te);
	        N_sample_max = max(te-ts+1, N_sample_max);
	        nspks = nspks + sum(truesp{i});
	    end
	    runs = i-1;
	    Pad = 2^(1+nextpow2(N_sample_max));  
	    f=linspace(0, 1, Pad)*RefreshRate;

	    Coh_cpl = (zeros(runs, Pad));
	    Coh = (zeros(runs, Pad));
	    Power = (zeros(runs, Pad));
	    Power_cpl = (zeros(runs, Pad));
	    Power_true = (zeros(runs, Pad));
	
	    for i = 1:runs
	        %Compute coherence of all runs
	        [coh_cpl, pwr_true, pwr_cpl, k] = coherence(truesp{i}, simsp_cpl{i}, N_sample_max);
    	    [coh, pwr1, pwr]  = coherence(truesp{i}, simsp{i}, N_sample_max);
        	Coh_cpl(i,:) = coh_cpl;
	        Coh(i,:) = coh;
    	    Power(i,:) = pwr;
        	Power_cpl(i,:) = pwr_cpl;
        	Power_true(i,:) = pwr_true;
	        K(i) = k;
	    end
	
		Power = mean(Power);
		Power_cpl = mean(Power_cpl);
		Power_true = mean(Power_true);

	    %Only compute for not NaN results
	    Coh = Coh(~isnan(Coh(:,1)),:);
	    if size(Coh,1) > 1
	        [Coh, Coh_SE] = jackknife_coh(Coh);
	    end
	
	    Coh_cpl = Coh_cpl(~isnan(Coh_cpl(:,1)),:);
	    if size(Coh_cpl,1) > 1
	        [Coh_cpl, Coh_cpl_SE] = jackknife_coh(Coh_cpl);
	    end
	    K = mean(K(~isnan(Coh(:,1))));
	    DF = runs*K-1; Null = sqrt(1-0.05^(1/DF));
	    L = 5; 
	    ii = 1:fix(Pad/L);
	    ii_s = 6;

	    coh1hz_cpl = abs(Coh_cpl);
	    coh1hz = abs(Coh);
	    meancoh1hz_cpl(icell) = mean(coh1hz_cpl(2:11));
	    meancoh1hz(icell) = mean(coh1hz(2:11));

		coh_05(icell) = Coh(ii_s);
		coh_cpl_05(icell) = Coh_cpl(ii_s);
		coh_SE_05(icell) = Coh_SE(ii_s);
		coh_SE_cpl_05(icell) = Coh_cpl_SE(ii_s);

	    figure(h);
	    subplot(nU,1,icell)
		hold on
		area(f(ii), abs(Coh_cpl(ii))+Coh_cpl_SE(ii), 'FaceColor', [0.7 .7 .7])
		area(f(ii), abs(Coh_cpl(ii))-Coh_cpl_SE(ii), 'FaceColor', [1 1 1])
		area(f(ii), abs(Coh(ii))+Coh_SE(ii), 'FaceColor', [0.7 .7 .7])
		area(f(ii), abs(Coh(ii))-Coh_SE(ii), 'FaceColor', [1 1 1])
	    h1 = plot(f(ii), abs(Coh(ii)), 'Color', [0 0 0.6]);
	    h2 = plot(f(ii), abs(Coh_cpl(ii)), 'Color', [0 .6 0]);
	    plot(f(ii), Null*ones(size(ii)), 'r');
	    legend([h1; h2], {'Uncoupled', 'Coupled'})
	    title(['Unit: ' num2str(icell) ' no. spikes: ' num2str(nspks)])
	    ylabel('<|coh|>')
	    ylim([0 1])
	    figure(g)
	    subplot(nU,1,icell)
   	    plot(f(ii), [Power(ii); Power_cpl(ii); Power_true(ii)]);
   	    legend('Uncoupled', 'Coupled', 'True')
	    title(['Unit: ' num2str(icell) ' no. spikes: ' num2str(nspks)])
   	    ylabel('power')
	end
	saveplot(h, [wd fn_out '.eps'], 'eps', [6 12])
	saveplot(g, [wd fn_out '_power.eps'], 'eps', [6 12])

	%Plot scatter plot:
	clf
	xlim([0 1])
	ylim([0 1])
	plot([0 1], [0 1], 'k')
	hold on
	errorbar(coh_05, coh_cpl_05, coh_SE_cpl_05, '.')
	herrorbar(coh_05, coh_cpl_05, coh_SE_05, '.')
	xlabel('Coherence b/w true spikes and uncoupled GLM')
	ylabel('Coherence b/w true spikes and coupled GLM')
	saveplot(gcf, [wd fn_out '_scatter.eps'])

	%Plot scatter plot mean zero to one Hz coherence:
	clf
	xlim([0 1])
	ylim([0 1])
	plot([0 1], [0 1], 'k')
	hold on
	scatter(meancoh1hz, meancoh1hz_cpl)
	xlabel('Mean (0-1hz) Coherence b/w true spikes and uncoupled GLM')
	ylabel('Mean (0-1hz) Coherence b/w true spikes and coupled GLM')
	saveplot(gcf, [wd fn_out '_scatter_mean1hz.eps'])
end