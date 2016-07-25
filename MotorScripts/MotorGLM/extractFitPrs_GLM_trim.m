function prs = extractFitPrs_GLM_trim(gg,Stim,MAXSIZE,processed, trim, offset);
% prs = extractFitPrs_GLM(gg,Stim,MAXSIZE);
%
% Set global variables for fitting and extract the 
% the parameters needed for fitting the (vanilla) GLM model

setupfitting_GLM_trim(gg,Stim,MAXSIZE,processed, trim, offset);  % Precompute quantities for optimization

prs = [gg.kt(:); gg.dc; gg.ih(:); gg.ih2(:)];

