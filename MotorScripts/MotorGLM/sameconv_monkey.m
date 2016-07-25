function G = sameconv(A, B, offset);
	%  G = sameconv(A, B);
	%   
	%  Causally filters A with B, giving a column vector with same height as
	%  A.  (B not flipped as in standard convolution).
	%
	%  Convolution performed efficiently in (zero-padded) Fourier domain.
	
	[am, an] = size(A);
	[bm, bn] = size(B);
	nn = am+bm-1;
	
	G = ifft(sum(fft(A,nn).*fft(flipud(B),nn),2));
	%This will make stim before stim
	if offset == 0
		G = G(1:am,:);
	elseif offset == 1
		%This will make 'stim' half before, half after spikes
		G = G(ceil(bm/2):ceil(bm/2)+am-1,:);	
	elseif offset == 2
		%This will make 'stim' after spikes
		G = G(end-am+1:end,:);
	end
