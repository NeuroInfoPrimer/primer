load([wd '/script_304_MotorData_PredictionValidation.mat'])

for l = 1:length(lambdas)
    l
    lambda = lambdas(l);
    coh_out = ['coherence_lambda_' num2str(lambda)];
    jackknifecoherence(wd, ['/preprocessed_networkglm_sims_lambda_' num2str(lambda) '.mat'], coh_out)
end

clf
semilogx(lambdas(3:end), sum(logl_glm(3:end,:),2))
xlabel('\lambda');
ylabel('Total coupled log-likelihood')
saveplot(gcf, [wd '/GLM_loglikelihood_compare_semilog.eps'])

clf
maxlogl = max(logl_glm(3:end-1,:), [], 1);
uncoupledlogl = logl_glm_uncoupled;
plot(uncoupledlogl, maxlogl, 'o')
hold on
plot([-.4 0], [-.4 0], 'r')
xlabel('Uncoupled log-likelihood');
ylabel('Coupled log-likelihood')
saveplot(gcf, [wd '/GLM_loglikelihood_compare.eps'])