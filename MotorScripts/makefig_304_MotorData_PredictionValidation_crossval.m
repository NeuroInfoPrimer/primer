load([wd '/script_304_MotorData_PredictionValidation_crossval.mat'])

clf
semilogx(lambdas(3:end), sum(mu_logl_glm(3:end,:),2))
hold on 
n = size(mu_logl_glm, 2);
semilogx(lambdas(3:end), sum(mu_logl_glm(3:end,:),2) + sum(std_logl_glm(3:end,:),2)/sqrt(n), '-')
semilogx(lambdas(3:end), sum(mu_logl_glm(3:end,:),2) - sum(std_logl_glm(3:end,:),2)/sqrt(n), '-')
%semilogx(lambdas(3:end), sum(mu_logl_glm_unc(3:end),2)*ones(size(lambdas(3:end))), '-.')
xlabel('\lambda');
ylabel('Total coupled log-likelihood')
saveplot(gcf, [wd '/GLM_loglikelihood_compare_semilog_crossval.eps'])

clf
maxlogl = max(mu_logl_glm(3:end-1,:), [], 1);
maxstd_logl_glm = max(std_logl_glm(3:end-1,:), [], 1);
uncoupledlogl = mu_logl_glm_unc;
%plot(uncoupledlogl, maxlogl, 'o')
herrorbar(uncoupledlogl, maxlogl, std_logl_glm_unc/sqrt(n), '.');
hold on
errorbar(uncoupledlogl, maxlogl, maxstd_logl_glm/sqrt(n), '.');
plot([-.4 0], [-.4 0], 'r')
xlabel('Uncoupled log-likelihood');
ylabel('Coupled log-likelihood')
saveplot(gcf, [wd '/GLM_loglikelihood_compare_crossval.eps'])

%Plot coherence...
for l = 1:length(lambdas)
    lambda = lambdas(l);
    display(['Lambda: ' num2str(lambda)])
    jackknifecoherence_crossval(wd, fn_out, nfolds, lambda)
end