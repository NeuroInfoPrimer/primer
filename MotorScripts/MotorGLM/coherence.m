function [Coh, Pwr1_tot, Pwr2_tot, K] = coherence(truesp, simsp, N_sample_max)
	N_sample=length(truesp) ; 
	F_sample=100; % Sampling rate
	Pad = 2^(1+nextpow2(N_sample_max));  
	% pad to > 2-times data length -
	f=linspace( 0, 1, Pad )*F_sample; %Frequency base
	NW= N_sample*0.5*(1/50); %Time-bandwidth product with bandwidth as a fraction (chose as 1/25) of f_Nyquist (0.5 in computer units) on a sample-by-sample basis
	K=1;fix(2*NW-1); % Degrees of freedom or 2*p-1
	[E,V] = dpss(N_sample,NW,K); % Family of multiple taper functions
	FT1=zeros(Pad,K);
	FT2=zeros(Pad,K);
	Pwr1=zeros(Pad,K);
	Pwr2=zeros(Pad,K);
	XPwr=zeros(Pad,K);
	Pwr1_tot=zeros(Pad,1);
	Pwr2_tot=zeros(Pad,1);
	XPwr_tot=zeros(Pad,1);
	Coh=zeros(Pad,1);
    z=1.0;
	for k=1:K; % sum over individual estimates of Power, with normalization corrections
	        %FT1(:,k)=fft(E(:,k).*(truesp-mean(truesp)),Pad);
	        %FT2(:,k)=fft(E(:,k).*(simsp-mean(simsp)),Pad);
	        FT1(:,k)=fft(E(:,k).*(truesp-z*mean(truesp)),Pad);
	        FT2(:,k)=fft(E(:,k).*(simsp-z*mean(simsp)),Pad);
	        Pwr1(:,k)=FT1(:,k).*conj(FT1(:,k));
	        Pwr2(:,k)=FT2(:,k).*conj(FT2(:,k));
	        XPwr(:,k)=FT1(:,k).*conj(FT2(:,k));
	        Pwr1_tot(:,1)=Pwr1_tot(:,1)+Pwr1(:,k);
	        Pwr2_tot(:,1)=Pwr2_tot(:,1)+Pwr2(:,k);
	        XPwr_tot(:,1)=XPwr_tot(:,1)+XPwr(:,k);
	end;
	Coh=XPwr_tot./sqrt(Pwr1_tot.*Pwr2_tot);
    Pwr1_tot=Pwr1_tot/K;
	Pwr2_tot=Pwr2_tot/K;
end