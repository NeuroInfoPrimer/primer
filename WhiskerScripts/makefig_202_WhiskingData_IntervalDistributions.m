% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_201 and generates plots that are
% shown in Figure 6F of the manuscript: histograms of inter-whisk-intervals
% and inter-spike-intervals (during whisking)

clear ; 

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
sr = 1000 ;
iwi = [] ;
isi = [] ;
for i = 1:7
    load([datadir 'VPM_cell_' Names{i} '_iswhisking.mat']) ;
    for ik = 1:k
        iwi = [iwi ; diff(isw_itop{ik})/sr] ;
        isi = [isi ; diff(isw_ispk{ik})/sr] ;
    end
end

isi(isi==0) = [] ;
bnw = linspace(0,0.4,50) ;
bns = linspace(0,0.4,50) ;
piwi = histc(iwi,bnw) ; 
piwi = piwi/sum(piwi)/(bnw(2)-bnw(1)) ;
pisi = histc(isi,bns) ; 
pisi = pisi/sum(pisi)/(bns(2)-bns(1)) ;

plot(bnw,piwi,'r','LineWidth',3) ; hold on ;
plot(bns,pisi,'b','LineWidth',3) ; hold on ;
legend('P(IWI) [Hz]','P(ISI) [Hz]') ;
xlabel('time (s)','FontSize',20) ;
set(gca,'FontSize',20) ;
