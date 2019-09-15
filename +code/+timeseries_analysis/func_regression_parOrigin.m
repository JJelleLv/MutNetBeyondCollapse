function [mdlOrigin,regParOrigin,regParPvalsOrigin,regParConfOrigin,regRsqOrigin,regResidOrigin,regPvalOrigin] = func_regression_parOrigin(regXdata,regYdata)

%% regression - linear through origin
mdlOrigin=fitlm(regXdata,regYdata,'y~x1-1');

regParOrigin=mdlOrigin.Coefficients.Estimate;
regParPvalsOrigin=mdlOrigin.Coefficients.pValue; %% pValue model parameters
regParConfOrigin=coefCI(mdlOrigin); %% confidence interval

regRsqOrigin=mdlOrigin.Rsquared.Adjusted;
regPvalOrigin=mdlOrigin.coefTest; %% geeft een soort overall p-waarde (zelfde als p-waarde van 'b' wanneer Y=a+bX)

regResidOrigin=mdlOrigin.Residuals.Raw;