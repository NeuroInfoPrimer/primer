% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_105, 106. MNE results for the
% retina dataset are not shown in manuscript. Generated here are plots
% showing:
% 1. eigenvalues of J and extent of null eigenvalue distribution
% 2. the vector H (first order feature with logistic nonlinearity)
% 3. all significant eigenvalues of J (features with quadratic-logistic 
%    nonlinearity 

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
    for iL = 2 
        load([datadir 'Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'Retina_cell_' num2str(icell) '_mne_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        
        ylm = 1.2*max(abs([eJ(:) ; eJnull(:)]))*[-1 1] ;
        
        figure ;
        
        subplot(nsig+1,2,[1 3]) ;
        fill([0 0 N N 0],[min_enull max_enull max_enull min_enull min_enull],0.8*ones(1,3),'EdgeColor','none') ; hold on ;
        plot(setdiff(1:N,isig),eJ(setdiff(1:N,isig)),'ok','MarkerFaceColor','k','MarkerSize',4) ; hold on ;
        plot(isig,eJ(isig),'or','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
        ylim(ylm) ;
        xlim([0 N]) ;
        set(gca, 'Layer', 'top','FontSize',fntsz) ;
        title([num2str(icell) ', ' stim_length{iL}],'FontSize',fntsz) ;
        ylabel('\lambda(\Delta C)','FontSize',fntsz) ;
        xlabel('eigenvalue #','FontSize',fntsz) ;
        axis square ;
        
        
        subplot(nsig+1,2,2) ;
        cm = max(max(abs(H)),max(abs(vJ(:))))/2 ;
        imagesc(1:sqrt(NX)*NT,1:sqrt(NX),reshape(H, [sqrt(NX) sqrt(NX)*NT]), [-cm cm]); hold on
        plot(0.5+[0 0 NT*sqrt(NX) NT*sqrt(NX) 0 ],0.5+[0 sqrt(NX) sqrt(NX) 0 0],'k') ; hold on ;
        for i = 1:NT-1
            plot(i*sqrt(NX)*[1 1]+0.5,0.5+[0 sqrt(NX)],'k') ; hold on ;
        end
        colorbar ;
        colormap hot;
        axis image xy off ;
    
   
        for j =  1:nsig
            subplot(nsig+1,2,2+2*j) ;
            
            imagesc(1:sqrt(NX)*NT,1:sqrt(NX),reshape(vJ(:,isig(j)), [sqrt(NX) sqrt(NX)*NT]), [-cm cm]) ; hold on
            plot(0.5+[0 0 NT*sqrt(NX) NT*sqrt(NX) 0 ],0.5+[0 sqrt(NX) sqrt(NX) 0 0],'k') ; hold on ;
            for i = 1:NT-1
                plot(i*sqrt(NX)*[1 1]+0.5,0.5+[0 sqrt(NX)],'k') ; hold on ;
            end
            title(num2str(eJ(isig(j)),2),'FontSize',fntsz) ;
            colormap hot;
            colorbar ;
            axis image xy off ;
        end
    end
end
