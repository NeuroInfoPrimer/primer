wd = '../MotorData/';
datafile = 'mabel.mat';
method = 'spg';
binsize = 0.01;
dt = .1;
frames = 80;
processed = preprocess([wd datafile], binsize, dt, frames);
uncoupled = load([wd '/all_units.mat']);
load([wd '/all_units_network_L1_method_' method '.mat']);

for l = 1:length(lambdas)
    lambda = lambdas(l);
    plot_filters_network_all(ggs_cpl(l,:), processed,...
     [wd '/GLM_coupled_lambda_' num2str(lambda) '.eps']);
    plot_filters_network_compare(ggs_cpl(l,:), uncoupled.ggs, processed,...
     [wd '/GLM_coupled_compare_lambda_' num2str(lambda) '.eps']);
end