function [Coh_all, Coh_SE] = jackknife_coh(Coh)
	nS = size(Coh,1);
	Coh_all = abs(mean(Coh));
	Coh_t = repmat(Coh_all, nS, 1);
	Coh_dropped = [];
	for i = 1:nS
		ii = [1:(i-1), (i+1):nS];
		Coh_dropped(i,:) = abs(mean(Coh(ii, :)));
	end
	Coh_SE = sqrt((nS-1)/nS*sum((Coh_dropped - Coh_t).^2,1));
end