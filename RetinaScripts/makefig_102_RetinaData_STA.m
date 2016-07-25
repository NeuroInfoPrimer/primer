% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_102 and generates the plots shown
% in Figure 4 of the manuscript:
% panels A,B - STA feature using long and short stimulus configurations 
%              (NT = 3,6)
% panel  C   - distribution of stimulus projected on STA
% panel  D   - nonlinearity

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

stim_length = {'short2','short3','long'} ;
sr = 30 ;
dt = 1/sr ;
fntsz = 10 ;

iJK = 5 ; 
for icell = 3:3
    figure ; 
        
    for iL = 2:3
        load([datadir 'Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'RetinaCellParameters_' stim_length{iL} '.mat']) ;
        
        subplot(2,3,2+3*(iL-2)) ;
        plot(bns(1:end-1),g(1:end-1)/dt,'b','LineWidth',2) ; hold on ;
        axis tight square ; box off
        ylim([0 1.8/dt]) ;
        xlabel('S*sta','FontSize',fntsz) ;
        ylabel('firing rate [Hz]','FontSize',fntsz) ;
        set(gca,'FontSize',fntsz) ;
        if iL==2 
            title(['Cell #' num2str(icell)],'FontSize',fntsz) ;
        end
        
        subplot(2,3,3+3*(iL-2)) ;
        plot(ctr_plot,Pf_plot,'k','LineWidth',2) ; hold on ;
        axis tight square ; box off
        ylim([0 1]) ;
        xlabel('S*sta','FontSize',fntsz) ;
        ylabel('P(S*sta)','FontSize',fntsz) ;
        set(gca,'FontSize',fntsz) ;
        
        subplot(2,3,1+3*(iL-2)) ;
        cm = max(abs(sta)) ;
        imagesc(1:Nv(icell)*NT,1:Nv(icell),reshape(sta, [Nv(icell) Nv(icell)*NT]), [-cm cm]); hold on
        plot(0.5+[0 0 NT*Nv(icell) NT*Nv(icell) 0 ],0.5+[0 Nv(icell) Nv(icell) 0 0],'k') ; hold on ;
        for i = 1:NT-1
            plot(i*Nv(icell)*[1 1]+0.5,0.5+[0 Nv(icell)],'k') ; hold on ;
        end
        colorbar ;
        colormap hot;
        axis image xy off ;
     end
end    