function [gg, fval,H] = MLfit_GLM_trim(gg,Stim,optimArgs,processed, trim, offset, lambda, method);
%  [ggnew,fval,H] = MLfit_GLM(gg,Stim,optimArgs);
% 
%  Computes the ML estimate for GLM params, using grad and hessians.
%  Assumes basis for temporal dimensions of stim filter
%
%  Inputs: 
%     gg = param struct
%     Stim = stimulus
%     optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%     ggnew = new param struct (with ML params);
%     fval = negative log-likelihood at ML estimate

%Only include times within stimulus that are within trial

MAXSIZE  = 1e7;  % Maximum amount to be held in memory at once;
if (nargin < 6) offset = 1; end
if (nargin < 7) lambda = 1; end
if (nargin < 8) method = 'spg'; end

thresh = 1e-5;

% Set optimization parameters 
if nargin > 2
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on','display','iter');
end

% Set initial params
prs0 = extractFitPrs_GLM_trim(gg,Stim,MAXSIZE,processed, trim, offset);

% minimize negative log likelihood. The warm-up solution for L1 reg
[prs,fval] = fminunc(@(p) Loss_GLM_logli(p),prs0,opts);

groups = zeros(size(prs));
nG = size(gg.tsp2,2);
nP = size(gg.ihbas2,2);

%gOptions.maxIter = 1000;
%lambda = 30;
gOptions.maxIter = 500;
gOptions.verbose = 2; % Set to 0 to turn off output
gOptions.corrections = 10; % Number of corrections to store for L-BFGS methods
gOptions.norm = 2; % Set to inf to use infinity norm
options = gOptions;
istart = size(gg.kt,1)*size(gg.kt,2)+1+size(gg.ih,1);
for idx = 1:nG
	indices = istart + (((idx-1)*nP+1):(idx*nP));
	groups(indices) = idx; 
end
lambdaVect = lambda*ones(nG, 1);	

options.method = method;
if ismember(method, {'spg', 'opg', 'pqn'})
	w = L1GeneralGroup_Auxiliary(@Loss_GLM_logli,prs,lambdaVect,groups,options);
elseif ismember(method, {'bbst', 'qnst'})
	w = L1GeneralGroup_SoftThresh(@Loss_GLM_logli,prs,lambdaVect,groups,options);
else
	error('Invalid L1 optimizationn method. Please see help.')
end

% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLM(gg,w);

%Check which coupling terms are zero and mask these ones out, then redo the fitting with just MLE
%Anything less than 10^-5 counts as zero
gg.mask = any(abs(gg.ih2) > thresh,1);

% Set initial params
prs1 = extractFitPrs_GLM_trim(gg,Stim,MAXSIZE,processed, trim, offset);
[prs,fval] = fminunc(@(p) Loss_GLM_logli(p, gg.mask), prs1, opts);

if nargout > 2 % Compute Hessian if desired
    [fv,gradval,H] = Loss_GLM_logli(prs);
end

% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLM(gg,prs);
