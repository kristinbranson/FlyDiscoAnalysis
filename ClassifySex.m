function [sex,diagnostics] = ClassifySex(area,mu_area,var_area,ptrans,state2sex)

x_trans_fun = @(areacurr,areaprev,scurr) normpdf(areacurr,mu_area(scurr),var_area(scurr));
log_ptrans = log(ptrans);
log_s_trans_fun = @(scurr,sprev) log_ptrans(sprev,scurr);

[sexstate,normscore] = first_order_viterbi(area(:),2,x_trans_fun,log_s_trans_fun);
sex = state2sex(sexstate);

if nargout > 1,
  diagnostics = struct;
  diagnostics.normhmmscore = normscore;
  diagnostics.nswaps = nnz(sexstate(1:end-1)~=sexstate(2:end));
  diagnostics.meanabsdev = nanmean(abs( area - mu_area(sexstate)' ));  
end