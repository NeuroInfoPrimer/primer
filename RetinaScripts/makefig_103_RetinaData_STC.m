% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_102, 103 and generates the plots 
% shown in Figure 5 of the manuscript: 
% panel A - STA and STC modes (after projecting out STA)
% panel B - eigenvalues of dC and extent of null eigenvalue distribution
% panel C - distribution of stimulus projected on STA and STC, full
%           nonlinearity and marginals

% Additionally this script generates plots that show: 
% 1. the STC features before and after projecting the STA out, and the
%    corresponding overlap
% 2. if kmodel (the maximal dimension of the model) is 3, a 3 dimensional
%    plot of the nonlinearity is shown, with surfaces for predicted firing
%    rates equal to 20% and 50% of the max
% 3. if kmodel is 2, the Signular Value Decomposition of the nonlinearity
%    is shown, together with the fraction of the full model that is
%    explained by it.

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

stim_length = {'short2','short3','long'} ;
sr = 30 ;
dt = 1/sr ;
fntsz = 10 ;

iJK = 5 ;

for icell =  3:3
    
    for iL = 3:3
        
        load([datadir 'Retina_cell_' num2str(icell) '_stc_'    stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'Retina_cell_' num2str(icell) '_stcsig_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'RetinaCellParameters_' stim_length{iL} '.mat']) ;
        
        N = NT*Nv(icell)^2 ;
        NX = N/NT ;
        
        switch kmodel
            case 2
                
                figure ;
                
                subplot(2,1,1) ;
                ylm = 1.2*max(abs([edC(:) ; edCnull(:)]))*[-1 1] ;
                fill([0 0 N N 0],[min_enull max_enull max_enull min_enull min_enull],0.8*ones(1,3),'EdgeColor','none') ; hold on ;
                plot(setdiff(1:N,isig),edC(setdiff(1:N,isig)),'ok','MarkerFaceColor','k','MarkerSize',4) ; hold on ;
                plot(isig,edC(isig),'or','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
                plot(1:N,prctenull,'Color',0.9*ones(1,3),'LineWidth',1) ;
                ylim(ylm) ;
                xlim([0 N]) ;
                set(gca, 'Layer', 'top','FontSize',fntsz) ;
                ylabel('\lambda(\Delta C)','FontSize',fntsz) ;
                xlabel('eigenvalue #','FontSize',fntsz) ;
                axis square ;
                title(num2str(icell),'FontSize',fntsz) ;
                
                subplot(2,1,2) ;
                imagesc(bns_plot,bns_plot,stc_model_plot/dt,[0 max(stc_model_plot(:))/dt]) ;
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('S*sta','FontSize',fntsz) ;
                colorbar ;
                xlim(ctr(end)*[-1 1]) ;
                ylim(ctr(end)*[-1 1]) ;
                axis square xy
                colormap('cool') ;
                set(gca,'FontSize',fntsz) ;
                
                figure ;
                
                subplot(2,2,1) ;
                plot(bns,g_a/dt,'b','LineWidth',2) ; hold on ;
                axis tight square ; box off
                ylim([0 1.5/dt]) ;
                xlabel('S*sta','FontSize',fntsz) ;
                ylabel('firing rate [Hz]','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                subplot(2,2,3) ;
                plot(bns,g_c/dt,'b','LineWidth',2) ; hold on
                axis tight square ; box off
                ylim([0 1.5/dt]) ;
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('firing rate [Hz]','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                subplot(2,2,2) ;
                plot(ctr_plot,Pf_a_plot,'k','LineWidth',2) ; hold on ;
                axis tight square ; box off
                ylim([0 1]) ;
                xlabel('S*sta','FontSize',fntsz) ;
                ylabel('P(S*sta)','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                subplot(2,2,4) ;
                plot(ctr_plot,Pf_c_plot,'k','LineWidth',2) ; hold on ;
                axis tight square ; box off
                ylim([0 1]) ;
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('P(S*stc_1)','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                figure ;
                
                imagesc(bns_plot,bns_plot,g_a*g_c'/dt,[0 1/dt]) ;
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('S*sta','FontSize',fntsz) ;
                colorbar ;
                xlim(ctr(end)*[-1 1]) ;
                ylim(ctr(end)*[-1 1]) ;
                axis square xy
                colormap('cool') ;
                
                f1 = reshape(f(:,1),Nv(icell)^2,NT) ;
                [f1_U,f1_S,f1_V] = svd(f1) ;
                f2 = reshape(f(:,2),Nv(icell)^2,NT) ;
                [f2_U,f2_S,f2_V] = svd(f2) ;
                
                cm = max(abs(f1_U(:))) ;
                
                subplot(2,3,1) ;
                plot(-dt*(0:1:NT-1),flipud(f1_V(:,1)),'k') ;
                xlim([-dt*(NT-1) 0]) ;
                title('temporal component') ;
                
                subplot(2,3,2) ;
                imagesc(1:sqrt(NX),1:sqrt(NX),reshape(f1_U(:,1),sqrt(NX),sqrt(NX)),[-cm cm]/2) ;
                axis xy square ;
                colormap('hot') ;
                title('spatial component') ;
                colormap hot;
                axis image xy ;
                
                subplot(2,3,3) ;
                imagesc(1:Nv(icell)*NT,1:Nv(icell),reshape(f1_U(:,1)*f1_V(:,1)', [Nv(icell) Nv(icell)*NT]), [-cm cm]); hold on
                plot(0.5+[0 0 NT*Nv(icell) NT*Nv(icell) 0 ],0.5+[0 Nv(icell) Nv(icell) 0 0],'k') ; hold on ;
                for i = 1:NT-1
                    plot(i*Nv(icell)*[1 1]+0.5,0.5+[0 Nv(icell)],'k') ; hold on ;
                end
                colormap hot;
                axis image xy off ;
                title(['outer product, ' num2str(f1_S(1,1)/sum(diag(f1_S)),3) ' of STA explained']) ;
                
                subplot(2,3,4) ;
                plot(-dt*(0:1:NT-1),flipud(f2_V(:,1)),'k') ;
                xlim([-dt*(NT-1) 0]) ;
                title('temporal component') ;
                
                subplot(2,3,5) ;
                imagesc(1:sqrt(NX),1:sqrt(NX),reshape(f2_U(:,1),sqrt(NX),sqrt(NX)),[-cm cm]/2) ;
                axis xy square ;
                colormap('hot') ;
                title('spatial component') ;
                colormap hot;
                axis image xy ;
                
                subplot(2,3,6) ;
                imagesc(1:Nv(icell)*NT,1:Nv(icell),reshape(f2_U(:,1)*f2_V(:,1)', [Nv(icell) Nv(icell)*NT]), [-cm cm]); hold on
                plot(0.5+[0 0 NT*Nv(icell) NT*Nv(icell) 0 ],0.5+[0 Nv(icell) Nv(icell) 0 0],'k') ; hold on ;
                for i = 1:NT-1
                    plot(i*Nv(icell)*[1 1]+0.5,0.5+[0 Nv(icell)],'k') ; hold on ;
                end
                colormap hot;
                axis image xy off ;
                title(['outer product, ' num2str(f2_S(1,1)/sum(diag(f2_S)),3) ' of STC_1 explained']) ;
                
            case 3
                
                subplot(2,1,1) ;
                ylm = 1.2*max(abs([edC(:) ; edCnull(:)]))*[-1 1] ;
                fill([0 0 N N 0],[min_enull max_enull max_enull min_enull min_enull],0.8*ones(1,3),'EdgeColor','none') ; hold on ;
                plot(setdiff(1:N,isig),edC(setdiff(1:N,isig)),'ok','MarkerFaceColor','k','MarkerSize',4) ; hold on ;
                plot(isig,edC(isig),'or','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
                plot(1:N,prctenull,'Color',0.9*ones(1,3),'LineWidth',1) ;
                ylim(ylm) ;
                xlim([0 N]) ;
                set(gca, 'Layer', 'top','FontSize',fntsz) ;
                ylabel('\lambda(\Delta C)','FontSize',fntsz) ;
                xlabel('eigenvalue #','FontSize',fntsz) ;
                axis square ;
                title(num2str(icell),'FontSize',fntsz) ;
                
                subplot(2,1,2) ;
                prct02 = max(stc_model_plot(:))*0.2 ;
                prct05 = max(stc_model_plot(:))*0.5 ;
                [x1,x2,x3] = ndgrid(bns_plot,bns_plot,bns_plot) ;
                iso02 = isosurface(x1,x2,x3,stc_model_plot/dt,prct02/dt) ;
                iso05 = isosurface(x1,x2,x3,stc_model_plot/dt,prct05/dt) ;
                p02 = patch(iso02); hold on ;
                p05 = patch(iso05); hold on ;
                set(p02,'FaceColor',[1 0 0],'EdgeColor','none','FaceAlpha',0.3) ;
                set(p05,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',0.3) ;
                legend([p02 p05],'20% of max FR','50% of max FR') ;
                box on ; axis square  ;
                xlabel('S*sta','FontSize',fntsz) ;
                ylabel('S*stc_1','FontSize',fntsz) ;
                zlabel('S*stc_2','FontSize',fntsz) ;
                view([20 -40]) ;
                set(gca,'FontSize',fntsz) ;
                
                figure ;
                
                subplot(3,2,1) ;
                plot(bns,g_a/dt,'b','LineWidth',2) ; hold on ;
                axis tight square ; box off
                ylim([0 1.5/dt]) ;
                xlabel('S*sta','FontSize',fntsz) ;
                ylabel('firing rate [Hz]','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                subplot(3,2,3) ;
                plot(bns,g_c1/dt,'b','LineWidth',2) ; hold on
                axis tight square ; box off
                ylim([0 1.5/dt]) ;
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('firing rate [Hz]','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                subplot(3,2,5) ;
                plot(bns,g_c2/dt,'b','LineWidth',2) ; hold on
                axis tight square ; box off
                ylim([0 1.5/dt]) ;
                xlabel('S*stc_2','FontSize',fntsz) ;
                ylabel('firing rate [Hz]','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                subplot(3,2,2) ;
                plot(ctr_plot,Pf_a_plot,'k','LineWidth',2) ; hold on ;
                axis tight square ; box off
                ylim([0 1]) ;
                xlabel('S*sta','FontSize',fntsz) ;
                ylabel('P(S*sta)','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                subplot(3,2,4) ;
                plot(ctr_plot,Pf_c1_plot,'k','LineWidth',2) ; hold on ;
                axis tight square ; box off
                ylim([0 1]) ;
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('P(S*stc_1)','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
                subplot(3,2,6) ;
                plot(ctr_plot,Pf_c2_plot,'k','LineWidth',2) ; hold on ;
                axis tight square ; box off
                ylim([0 1]) ;
                xlabel('S*stc_2','FontSize',fntsz) ;
                ylabel('P(S*stc_2)','FontSize',fntsz) ;
                set(gca,'FontSize',fntsz) ;
                
        end
        
        figure ;
        
        cm = max(abs(f(:))) ;
        for k = 1:kmodel
            subplot(kmodel,1,k) ;
            imagesc(1:Nv(icell)*NT,1:Nv(icell),reshape(f(:,k), [Nv(icell) Nv(icell)*NT]), [-cm cm]); hold on
            plot(0.5+[0 0 NT*Nv(icell) NT*Nv(icell) 0 ],0.5+[0 Nv(icell) Nv(icell) 0 0],'k') ; hold on ;
            for i = 1:NT-1
                plot(i*Nv(icell)*[1 1]+0.5,0.5+[0 Nv(icell)],'k') ; hold on ;
            end
            colormap hot;
            axis image xy off ;
            if k == 1
                title(num2str(icell),'FontSize',fntsz) ;
                ylabel('STA','FontSize',fntsz) ;
            else
                ylabel(['STC_' num2str(k-1)],'FontSize',fntsz) ;
            end
            colorbar ;
        end
        
        
        figure ;
        vdC_sta = zeros(N,nsig) ;
        for k = nsig:-1:1
            subplot(nsig,2,k*2-1) ;
            imagesc(1:Nv(icell)*NT,1:Nv(icell),reshape(vdC(:,isig(k))*sign(sum(f(:,1)))*sign(sum(vdC(:,isig(k)))), [Nv(icell) Nv(icell)*NT]), [-cm cm]); hold on
            plot(0.5+[0 0 NT*Nv(icell) NT*Nv(icell) 0 ],0.5+[0 Nv(icell) Nv(icell) 0 0],'k') ; hold on ;
            for i = 1:NT-1
                plot(i*Nv(icell)*[1 1]+0.5,0.5+[0 Nv(icell)],'k') ; hold on ;
            end
            colormap hot;
            axis image xy off ;
            title(['stc #' num2str(k) ' (\lambda(\Delta C) = ' num2str(edC(isig(k))) ', projection on sta = ' num2str(abs(vdC(:,isig(k))'*f(:,1)),4)],'FontSize',fntsz) ;
            
            vdC_sta(:,k) = vdC(:,isig(k)) - (vdC(:,isig(k))'*f(:,1))*f(:,1) ;
            vdC_sta(:,k) = vdC_sta(:,k)/norm(vdC_sta(:,k)) ;
            
            subplot(nsig,2,k*2) ;
            imagesc(1:Nv(icell)*NT,1:Nv(icell),reshape(vdC_sta(:,k), [Nv(icell) Nv(icell)*NT]), [-cm cm]); hold on
            plot(0.5+[0 0 NT*Nv(icell) NT*Nv(icell) 0 ],0.5+[0 Nv(icell) Nv(icell) 0 0],'k') ; hold on ;
            for i = 1:NT-1
                plot(i*Nv(icell)*[1 1]+0.5,0.5+[0 Nv(icell)],'k') ; hold on ;
            end
            colormap hot;
            axis image xy off ;
            title(['stc #' num2str(k) ', sta projected out'],'FontSize',fntsz) ;
            colorbar ;
        end
    end
end