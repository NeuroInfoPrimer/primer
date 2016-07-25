function jackknifecoherence_crossval(wd, fn_out, nF, lambda)
	%nF = n folds
	nU = 9;
    meancoh1hz_cpl_all = zeros(nU, nF);
    meancoh1hz_all = zeros(nU, nF);
	for f_idx = 1:nF
		display(['Loading fold: ' num2str(f_idx)])
		h = figure;
		g = figure;
		coh_05 = zeros(nU,1);
		coh_cpl_05 = zeros(nU,1);
		coh_SE_05 = zeros(nU,1);
		coh_SE_cpl_05 = zeros(nU,1);
		fn_in = ['/preprocessed_networkglm_sims_lambda_' num2str(lambda) '_fold_' num2str(f_idx) '.mat'];
		load([wd fn_in])
		maxt = 100000;
		for icell = 1:nU
		    truesp = {};
		    simsp = {};
		    simsp_cpl = {};
		    runs = size(proc_withheld.trialstartend,1);
		    K = [];
		    N_sample_max = 0;
		    %fn = [wd '/GLM_cell_simulation_' num2str(icell) '_fold_' num2str(f_idx) '.mat'];
			uncoupled = load([wd '/GLM_cell_simulation_' num2str(icell) '_fold_' num2str(f_idx) '.mat'])
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
	
		 	size(Coh_cpl)
		    coh1hz_cpl = abs(Coh_cpl);
		    coh1hz = abs(Coh);
		    meancoh1hz_cpl(icell) = mean(coh1hz_cpl(2:11));
		    meancoh1hz(icell) = mean(coh1hz(2:11));

		    meancoh1hz_cpl_all(icell,f_idx) = mean(coh1hz_cpl(2:11));
		    meancoh1hz_all(icell,f_idx) = mean(coh1hz(2:11));
	
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

	saveplot(h, [wd fn_out '_lambda_' num2str(lambda) '_fold_' num2str(f_idx) '.eps'], 'eps', [6 12])
	saveplot(g, [wd fn_out '_lambda_' num2str(lambda) '_fold_' num2str(f_idx) '_power.eps'], 'eps', [6 12])

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
	saveplot(gcf, [wd fn_out '_lambda_' num2str(lambda) '_scatter_fold_' num2str(f_idx) '.eps'])

	%Plot scatter plot mean zero to one Hz coherence:
	clf
	xlim([0 1])
	ylim([0 1])
	plot([0 1], [0 1], 'k')
	hold on
	scatter(meancoh1hz, meancoh1hz_cpl)
	xlabel('Mean (0-1hz) Coherence b/w true spikes and uncoupled GLM')
	ylabel('Mean (0-1hz) Coherence b/w true spikes and coupled GLM')
	saveplot(gcf, [wd fn_out '_lambda_' num2str(lambda) '_scatter_mean1hz_fold_' num2str(f_idx) '.eps'])
	end

	%Make the final plot...	
	clf
	hold on 
	plot([0 0], [1 1], 'r')
	xlim([0 1])
	ylim([0 1])
	%Compute circle of points
	nP = 200;
	pts = zeros(2, nP);
	pts(1,:) = cos((1:nP)/nP*2*pi);
	pts(2,:) = sin((1:nP)/nP*2*pi);
	for idx = 1:nU
		coh_cpl = meancoh1hz_cpl_all(idx,:);
		coh_uncpl = meancoh1hz_all(idx,:);
		data = [coh_uncpl', coh_cpl'];
		%Transform data: log(|c|^2)/(1-|c|^2|))
		d2 = data.*data;
		tdata = log(d2./(1-d2));
		%Compute statistics
		coh_mu_orig = mean(data);
		coh_cov_orig = cov(data);
		coh_mu = mean(tdata);
		coh_cov = cov(tdata);
		%Compute curves
		[uo, so, vo] = svd(coh_cov_orig);
		[u, s, v] = svd(coh_cov);
		rpts = repmat(coh_mu_orig', 1, nP)+1*uo*sqrt(so)*pts;
		tpts = repmat(coh_mu', 1, nP)+1*u*sqrt(s)*pts;
		%transform back to [0..1]
		rtpts = sqrt(exp(tpts)./(1+exp(tpts)));
		%Plot
		plot(rtpts(1,:), rtpts(2,:))
		%plot_gaussian_ellipsoid(coh_mu, coh_cov_orig, 1);
		%plot(rpts(1,:), rpts(2,:), '--r')
		%plot(coh_mu_orig(1), coh_mu_orig(2),'+');
	end
	plot([0 0], [1 1], 'r')
	xlabel('Mean (0-1hz) coherence b/w true spikes and uncoupled GLM')
	ylabel('Mean (0-1hz) coherence b/w true spikes and coupled GLM')
	saveplot(gcf, [wd fn_out '_lambda_' num2str(lambda) '_scatter_mean1hz_all_1sigma.eps'])
	%savefig([wd fn_out '_lambda_' num2str(lambda) '_scatter_mean1hz_all.fig'])
	saveas(gcf, [wd fn_out '_lambda_' num2str(lambda) '_scatter_mean1hz_all_1sigma.fig'], 'fig')
end

