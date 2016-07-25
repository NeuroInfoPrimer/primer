% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_207 and generates
% plots that are shown in Figure 8 of the manuscript: 
% panel A - the MNE modes (first and second order).
%           h, the first order mode is compared to the STA
%           all modes are also shown filtered by the principal components
%           of the stimulus.
% panel B - the eigenvalue of J compared with the extent of the null
%           distribution

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
lgnd  = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;

sr = 1000 ;        % (Hz) sampling rate
N = 150 ;          % stimulus dimensionality 
ds = 2 ;           % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;       % delta t of stimulus 

a = 0 ;
k = 15 ;
tp = -dt*(N:-1:1) ;
fntsz = 10 ;

iJK = 5 ;
for i = 1:7
    load([datadir 'VPM_cell_' Names{i} '_mne_JK_' num2str(iJK) '.mat']) ; 
    load([datadir 'VPM_cell_' Names{i} '_sta_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
    
    S = S/N ;
    C = cov(S-repmat(mean(S,1),[T,1])) ;
    [vC,eC] = eig(C) ;
    [eC,iC] = sort(diag(eC),'descend') ;
    vC = vC(:,iC) ;
    
    Hpc = (H*vC(:,1:k))*vC(:,1:k)' ;
    vJpc = ((vJ'*vC(:,1:k))*vC(:,1:k)')' ;
    
    ylm = 1.2*max(abs([eJ(:) ; eJnull(:)]))*[-1 1] ;
    
    figure ; 
      
    subplot(max(nsig,2)+1,2,[1 3 5]) ;
    fill([0 0 N N 0],[min_enull max_enull max_enull min_enull min_enull],0.8*ones(1,3),'EdgeColor','none') ; hold on ;
    plot(setdiff(1:N,isig),eJ(setdiff(1:N,isig)),'ok','MarkerFaceColor','k','MarkerSize',4) ; hold on ;
    plot(isig,eJ(isig),'or','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
    ylim(ylm) ;
    xlim([0 N]) ;
    set(gca, 'Layer', 'top','FontSize',fntsz) ;
    title(lgnd{i},'FontSize',fntsz) ;
    ylabel('\lambda(\Delta C)','FontSize',fntsz) ;
    xlabel('eigenvalue #','FontSize',fntsz) ;
    axis square ;
    
    
    subplot(max(nsig,2)+1,2,2) ;
    h0 = plot(tp,zeros(size(tp)),'k') ; hold on ;
    h1 = plot(tp,-H/norm(H),'LineWidth',3) ; hold on ;
    h2 = plot(tp,-Hpc/norm(Hpc),'r','LineWidth',3) ; hold on ;
    h3 = plot(tp,sta/norm(sta),'g','LineWidth',3) ; axis tight ; box off
    set(gca,'FontSize',fntsz,'Yaxislocation','right','xtick',[]) ;
    legend([h1 h2 h3],'h','h filtered by PCA','STA')
    for j =  1:nsig
        subplot(max(nsig,2)+1,2,2+2*j) ;
        h0 = plot(tp,zeros(size(tp)),'k') ; hold on ;
        h1 = plot(tp,vJ(:,isig(j)),'LineWidth',3) ; hold on
        h2 = plot(tp,vJpc(:,isig(j)),'LineWidth',3) ; axis tight ; box off
        if j == 1
            legend([h1 h2],'eigenvector of J','filtered by PCA')
        end
        title(['\lambda(J)=' num2str(round(1e4*eJ(isig(j)))/1e4)],'FontSize',fntsz/1.5) ;
        if j<nsig
            set(gca,'FontSize',fntsz,'Yaxislocation','right','xtick',[]) ;    
        else
            set(gca,'FontSize',fntsz,'Yaxislocation','right') ;
            xlabel('time (s)','FontSize',fntsz) ;
        end
        
    end
end
