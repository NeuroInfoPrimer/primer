% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_107 and generates the plots shown
% in Figure 9 of the manuscript:
% panel A - stimulus filter (only GLM feature shown here, STA generated in
%           makefig_102
% panel B - spike history filter and exponentiated spike history filter

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

stim_length = {'short2','short3','long'} ;

sr = 30 ;          % (Hz) sampling rate
dt = 1/sr ;       % delta t of stimulus

fntsz = 15 ;

iJK = 5 ;

for icell = 3:3
    for iL = 3:3

        load([datadir 'Retina_cell_' num2str(icell) '_glm_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        
        figure ;

        subplot(3,1,1) ;
        cm = max(abs(gg.k(:))) ;
        Splot = zeros(sqrt(NX),sqrt(NX)*NT) ;
        for i = 1:NT
            Splot(1:sqrt(NX),(1:sqrt(NX))+(i-1)*sqrt(NX)) = reshape(gg.k(i,:),sqrt(NX),sqrt(NX)) ;
        end
        imagesc(1:sqrt(NX)*NT,1:sqrt(NX),Splot, [-cm cm]) ; hold on
        plot(0.5+[0 0 NT*sqrt(NX) NT*sqrt(NX) 0 ],0.5+[0 sqrt(NX) sqrt(NX) 0 0],'k') ; hold on ;
        for i = 1:NT-1
            plot(i*sqrt(NX)*[1 1]+0.5,0.5+[0 sqrt(NX)],'k') ; hold on ;
        end
        colorbar ; 
        colormap hot;
        axis image xy off ;
        set(gca,'FontSize',fntsz,'Yaxislocation','right') ;
        
        subplot(3,1,2) ;
        plot(dt*gg.iht,zeros(size(gg.iht)),'k') ; hold on ;
        plot(dt*gg.iht,gg.ihbas*gg.ih,'b','LineWidth',3) ; axis tight ; box off ;
        legend('GLM spike history filter') ;
        xlabel('time relative to spike (s)') ;
        set(gca,'FontSize',fntsz) ;
        
        subplot(3,1,3) ;
        plot(dt*gg.iht,ones(size(gg.iht)),'k') ; hold on ;
        plot(dt*gg.iht,exp(gg.ihbas*gg.ih),'b','LineWidth',3) ; axis tight ; box off ;
        legend('exp(GLM spike history filter)') ;
        xlabel('time relative to spike (s)') ;
        set(gca,'FontSize',fntsz) ;

    end
end