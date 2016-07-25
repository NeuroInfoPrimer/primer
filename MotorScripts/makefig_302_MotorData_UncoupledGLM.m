wd = '../MotorData/';
load([wd '/all_units.mat']);
datafile = 'mabel.mat';
method = 'spg';
binsize = 0.01;
dt = .1;
frames = 80;
proc = preprocess([wd datafile], binsize, dt, frames);
plot_filters(ggs, proc, [wd '/uncoupled_filters.eps']);