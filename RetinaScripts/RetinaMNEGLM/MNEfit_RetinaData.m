function pfinal = MNEfit_RetinaData(stim, resp, teststim, testresp, order)


% If you are fitting a MNE model with constraints beyond first and second
% order correlations, you will have to change the following:

% Calculate constrained averages
[Nsamples,Ndim] = size(stim);
psp = mean(mean(resp));   % spike probability

avg = (stim'*resp)/Nsamples;  % Ndim x Nrep

avg = mean(avg,2);  % Ndim x 1
avgs = [psp;avg];
if order>1
    avgsqrd = stim'*(repmat(resp,[1,Ndim]).*stim)/Nsamples;  % Ndim x Ndim
    avgsqrd = reshape(avgsqrd,[Ndim^2,1]);
    avgs = [avgs;avgsqrd];
end
% Initialize parameters
pstart = log(1/avgs(1)-1); % bias constant
pstart(2:Ndim+1) = .001*(2*rand([1,Ndim])-1); % h field
if order>1
    temp = .0005*(2*rand([Ndim,Ndim])-1);
    pstart(Ndim+2:length(pstart)+Ndim^2) = reshape((temp+temp'),[1,Ndim^2]);
    clear temp;
end



% Run conjugate gradient algorithm
%pfinal = frprmn(pstart, @logloss, @dlogloss, stim, resp, teststim, testresp, order, avgs, fittype);
pfinal = frprmn_global_min(pstart, @logloss, @dlogloss, stim, resp, teststim, testresp, order, avgs, 0);
