% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_201 and generates plots that are
% shown in Figure 6E of the manuscript:
%   The eigenvalues of the stimulus covariance matrix
%   The stimulus principal components (eigenvectors of stimulus covariance
%   matrix

clear ;
workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
lgnd  = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;

iJK = 5 ;
for i = 1:7
    figure ; 
    
    load([datadir 'VPM_cell_' Names{i} '_iswhisking.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_hilbert.mat']) ;
    S1 = S/mean(sqrt(sum(S.^2,2))) ;
    C = cov(S1) ;
    [vC,eC] = eig(C) ;
    [eC,iC] = sort(diag(eC),'descend') ;
    vC = vC(:,iC) ;

    semilogy(1:N,eC,'ok','MarkerSize',8,'MarkerFaceColor','w') ;
    xlabel('eigenvalue #','FontSize',20) ;
    ylabel('\lambda','FontSize',20) ;
    set(gca,'FontSize',20) ;
    title([ 'Stimulus covariance spectrum, cell ' lgnd{i}]) ;

    figure ;
    for j = 1:12 
        subplot(12,1,j) ;
        plot(ds*(-N:1:-1)/sr,vC(:,j)) ; 
        ylim([-0.2 0.2]) ;
        if j == 1
            title([ 'Stimulus principal components, cell ' lgnd{i}]) ;
        end

        set(gca,'FontSize',10) ;
    end
    
    xlabel('time (s)','FontSize',10) ;
    set(gca,'FontSize',10) ;
end