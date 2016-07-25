% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_101 and generates the plots shown
% in Figure 3 of the manuscript: 
% panel A - an illustration of the stimulus and response
% panel B - The spectrum of the stimulus covariance matrix and a comparison
%           to the Marcenko-Pastur distribution. The prediction of the MP
%           distribution is computed here.
% panel C - leading stimulus principal components (eigenvectors of stimulus
%           covariance

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

fntsz = 10 ;
stim_length = {'short2','short3','long'} ;
dt = 1/30 ;

iJK = 5 ;
for icell = 3:3
    for iL = 2:2
        
        figure ;
        
        load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        sig = 1 ;
        gam = T/N ;
        lCp = sig^2*(1+1/sqrt(gam))^2 ;
        lCm = sig^2*(1-1/sqrt(gam))^2 ;
        
        lv = linspace(lCm,lCp,N) ;
        dlv = lv(2)-lv(1) ;
        drho = 1/(2*pi*sig^2)*sqrt((lCp-lv).*(lv-lCm))./(lv./gam) ;
        crho = cumsum(drho*dlv) ;
        icrho = zeros(1,N) ;
        for i = 1:N
            [~,j] = min(abs(i/N - crho)) ;
            icrho(i) = lv(j) ;
        end
        icrho = sort(icrho,'descend') ;
        
        C = cov(S) ;
        [vC,eC] = eig(C) ;
        [eC,iC] = sort(diag(eC),'descend') ;
        vC = vC(:,iC) ;
        
        plot(1:N,eC,'ok','MarkerSize',8,'MarkerFaceColor','w') ; hold on 
        plot(1:N,icrho,'-r','LineWidth',2) ;
        xlabel('eigenvalue #','FontSize',fntsz) ;
        ylabel('\lambda','FontSize',fntsz) ;
        set(gca,'FontSize',fntsz) ;
        ylim([0 1.2]) ;
        
        figure ;
        for j = 1:12
            subplot(12,1,j) ;
            imagesc(reshape(vC(:,j), sqrt(NX), NT*sqrt(NX)), [-0.25 0.25]); hold on
            plot(0.5+[0 0 NT*sqrt(NX) NT*sqrt(NX) 0 ],0.5+[0 sqrt(NX) sqrt(NX) 0 0],'k') ; hold on ;
            for k = 1:NT-1
                plot(k*sqrt(NX)*[1 1]+0.5,0.5+[0 sqrt(NX)],'k') ; hold on ;
            end
            axis equal off ;
        end
        colormap('gray')
        
        figure ;
        j0 = 112 ;
        subplot(2,1,1) ;
        nplot = 20 ;
        Splot = zeros(sqrt(NX),sqrt(NX)*nplot) ;
        for j = 1:nplot
            Splot(1:sqrt(NX),(1:sqrt(NX))+(j-1)*sqrt(NX)) = reshape(S(j+j0,1:NX),sqrt(NX),sqrt(NX)) ;
        end
        imagesc(Splot, [-1 1]) ; hold on
        plot(0.5+[0 0 nplot*sqrt(NX) nplot*sqrt(NX) 0 ],0.5+[0 sqrt(NX) sqrt(NX) 0 0],'r','LineWidth',2) ; hold on ;
        for k = 1:nplot-1
            plot(k*sqrt(NX)*[1 1]+0.5,0.5+[0 sqrt(NX)],'r','LineWidth',2) ; hold on ;
        end
        axis equal off ;
        
        subplot(2,1,2) ;
        stem(dt*(1:T),R,'Color','k','Marker','none') ; hold on ;
        plot(dt*(j0+[1 nplot nplot 1 1]),[0 0 3 3 0],'r','LineWidth',2) ;
        xlabel('time [s]')
        ylabel('# of spikes')
        
        box off ;
        axis tight ;
        colormap('gray') ;
        xlim([2 12]) ;
        set(gca,'FontSize',20) ;
        
    end
end
