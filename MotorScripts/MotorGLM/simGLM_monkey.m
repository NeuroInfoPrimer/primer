function [tsp,Vmem,Ispk,conv] = simGLM(glmprs,Stim, time_limit, offset);
% [tsp, Vmem,Ispk] = simGLM(glmprs,Stim);
% 
% Compute response of glm to stimulus Stim.
%
% Uses time rescaling instead of Bernouli approximation to conditionally
% Poisson process
%
% Dynamics:  Filters the Stimulus with glmprs.k, passes this through a
% nonlinearity to obtain the point-process conditional intensity.  Add a
% post-spike current to the linear input after every spike.
if (nargin < 3) time_limit = 20; end
if (nargin < 4) offset = 1; end

if size(glmprs.k,3) > 1  % Run "coupled" GLM model if multiple cells present
    if nargout <= 2
        [tsp,Vmem,a, conv] = simGLMcpl_monkey(glmprs,Stim, time_limit, offset);
    else
        [tsp,Vmem,Ispk, conv] = simGLMcpl_monkey(glmprs,Stim, time_limit, offset);
    end
else
    if nargout <= 2
        [tsp,Vmem, a, conv] = simGLMsingle(glmprs,Stim, time_limit, offset);
    else
        [tsp,Vmem,Ispk, conv] = simGLMsingle(glmprs,Stim, time_limit, offset);
    end
end
    
% ======================================================
function [tsp,Vmem,Ispk, conv] = simGLMsingle(glmprs,Stim, time_limit, offset);
% Sub-function for simGLM.m
% Simulates the GLM point process model for a single (uncoupled) neuron
% using time-rescaling
    
% --------------- Check Inputs ----------------------------------
global RefreshRate;
conv = 1;
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dt;
% Check that dt evenly divides 1
if (mod(1,dt) ~= 0), dt = 1/round(1/dt);
    fprintf(1, 'glmrunmod: reset dtsim = %.3f\n', dt);
end
if ~isfield(glmprs, 'nlfun'), nlfun = @exp;
else, nlfun = glmprs.nlfun;
end

[slen,swid] = size(Stim); % length of stimulus
rlen = round(slen/dt);  % length of binned response

% -------------  Compute filtered resp to signal ----------------
k = glmprs.k;  % stim filter
if swid == size(k,2) % convolve if k is the same width as stimulus
    %Convolve components individually and add them at the end...
    Vstm = zeros(size(Stim,1),1);
    for idx = 1:swid
        v = sameconv_monkey(Stim(:,idx),k(:,idx), offset);
        Vstm = Vstm + v; 
    end
else 
    error('Mismatch between stimulus and filter size -- glmrunmod');
end
Vstm = [0; Vstm; 0]+glmprs.dc;

% -------------- Compute interpolated h current ----------------------
if isfield(glmprs,'ihbas');  % Check if h-kernel parametrized by a basis
    if isempty(glmprs.ih), ih = [];
    else, ih = glmprs.ihbas*glmprs.ih;
    end
else, ih = glmprs.ih;
end
if ~isempty(ih) % Interpolate ih current
    ihthi = [dt:dt:max(glmprs.iht)]';  % time points for sampling
    ihhi = interp1(glmprs.iht, ih, ihthi, 'linear', 0);
    hlen = length(ihhi);
else % No ih current
    hlen = 1; ihhi = 0; ihthi = dt;
    if ~isempty(glmprs.iht)
        fprintf('Warning: no post-spike kernel supplied but iht nonempty!!\n');
    end
end
    
% -------------  Static nonlinearity & spiking -------------------
nsp = 0;
tsp = zeros(round(slen/5),1);  % allocate space for spike times
Vmem = interp1([0:slen+1]',Vstm,[.5+dt:dt:slen+.5]', 'linear');
if (nargout > 2), Ispk = Vmem*0;
end

jbin = 1; % current time bin
nsp = 0; % number of spikes
tspnext = exprnd(1);  % time of next spike (in rescaled time)
rprev = 0;  % Integrated rescaled time up to current point
t0=clock;
prevspbins = zeros(300,1);
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    rrnxt = nlfun(Vmem(iinxt))*dt/RefreshRate; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike!
        ispk = iinxt(min(find(rrcum>=tspnext))); % time bin where spike occurred
        prevspbins = [ispk; prevspbins(1:end-1)];
        nsp = nsp+1;
        tsp(nsp,1) = ispk*dt; % spike time
        mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
        if ~isempty(iiPostSpk)
            Vmem(iiPostSpk) = Vmem(iiPostSpk)+ihhi(1:mxi-ispk);
            if nargout == 3  % Record post-spike current
                Ispk(iiPostSpk) = Ispk(iiPostSpk)+ihhi(1:mxi-ispk);
            end
        end
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
    if etime(clock,t0) > time_limit
        disp('Ran out of time, bud.');
        conv = 0;
        break;
    end
    if all(diff(prevspbins)==-1)
        disp('300 spiking bins in a row, suspect unstable GLM parameters, skipping')
        conv = 0;
        break;
    end
end
tsp = tsp(1:nsp); % prune extra zeros
