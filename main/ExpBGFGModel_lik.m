function [llr,ll_fore,ll_back] = ExpBGFGModel_lik(model,im)

d2 = (im - model.fg_mu).^2 ./ model.fg_sigma.^2 / 2;
ll_fore = -d2 - model.fg_log_Z;

d2 = (im - model.bg_mu).^2 ./ model.bg_sigma.^2 / 2;
ll_back = -d2 - model.bg_log_Z;

llr = ll_fore - ll_back;