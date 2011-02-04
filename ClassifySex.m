function sex = ClassifySex(area,mu_area,var_area,ptrans,state2sex)

x_trans_fun = @(areacurr,areaprev,scurr) normpdf(areacurr,mu_area(scurr),var_area);
log_ptrans = log(ptrans);
log_s_trans_fun = @(scurr,sprev) log_ptrans(sprev,scurr);

sexstate = first_order_viterbi(area(:),2,x_trans_fun,log_s_trans_fun);
sex = state2sex(sexstate);
