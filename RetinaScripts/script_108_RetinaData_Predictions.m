% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads all the models that were fitted for each cell and each
% stimulus configuration and runs them on the test portion of the stimulus
% to give a predicted spike train. 

% the portion of the code that computes the model prediction for the MNE
% model is in the comments below. This model was not included in the paper,
% but the user can fit the model (scripts 105, 106) and generate
% predictions

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

stim_length = {'short2','short3','long'} ;

nJK = 5 ;            % number of jack-knives 

for icell = 3:3
    for iL = 3:3
        for iJK = 1:nJK
            clearvars -except datadir icell iL stim_length iJK nJK
        
%             % Predicted spike train from MNE model
%             load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
%             load([datadir 'Retina_cell_' num2str(icell) '_mne_'       stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
%             SHt  = St*H' ;
%             SJSt = sum(St.*(St*J),2) ;
%             Rt_mne = 1./(1+exp(SJSt+SHt+A)) ;

%             % Predicted spike train from MNE model
%             load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{1} '_JK_' num2str(iJK) '.mat']) ;
%             load([datadir 'Retina_cell_' num2str(icell) '_mne_'       stim_length{1} '_JK_' num2str(iJK) '.mat']) ;
%             SHt  = St*H' ;
%             SJSt = sum(St.*(St*J),2) ;
%             Rt_mne2 = 1./(1+exp(SJSt+SHt+A)) ;
%             Rt_mne2 = Rt_mne2(length(Rt_mne2)-length(Rt_mne)+(1:length(Rt_mne))) ;
        
            % Predicted spike train from STA model
            load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            load([datadir 'Retina_cell_' num2str(icell) '_sta_'       stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            Rt_sta = sta_model_rect_norm(St*sta) ;
        
            % Predicted spike train from STC model
            load([datadir 'Retina_cell_' num2str(icell) '_stc_'       stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            zt = St*f ;
            Rt_stc = stc_model_rect_norm(zt) ;
            Rt_stc(isnan(Rt_stc)) = min(~isnan(Rt_stc)) ;
            
            % Predicted spike train from glm model
            load([datadir 'Retina_cell_' num2str(icell) '_glm_'       stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;

            save([datadir 'Retina_cell_' num2str(icell) '_' stim_length{iL} '_Predictions_JK_' num2str(iJK) '.mat'],'Rt','Rt_sta','Rt_glm','Rt_stc') ;
            % if MNE model is also fit, add Rt_mne and Rt_mne2 to list of
            % saved variables
        end
    end
end