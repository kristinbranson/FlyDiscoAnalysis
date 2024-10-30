function [trx, diagnostics, mu_area, var_area, state2sex, ll] = ...
  add_sex_field_to_trx_computed_using_machine_vision(trx, X, version, analysis_protocol, sexclassifier_params, sexclassifierin, dosave, sexclassifieroutmatfile)

%% learn a 2-state HMM for area in an unsupervised manner

fprintf('gender = "b", learning 2-state HMM...\n');

% initialize parameters
nstates = 2;
if isfield(sexclassifierin,'ptrans')
  fprintf(1,'''ptrans'' field present in SC settings...\n');

  ptrans = sexclassifierin.ptrans;
  psame = 1-ptrans;
else
  fprintf(1,'''ptrans'' field not present in SC settings. Using ''psame'' field...\n');
 
  psame = sexclassifierin.psame;
  ptrans = 1-psame;
end
fprintf(1,' ... setting [psame ptrans] = [%.03g %.03g]\n',psame,ptrans);
ptransmat = ones(nstates)*ptrans;
ptransmat(eye(nstates)==1) = psame;

if isfield(sexclassifier_params,'frac_female'),
  prior = [1-sexclassifier_params.frac_female,sexclassifier_params.frac_female];
else
  prior = ones(1,nstates)/nstates;
end
state2sex = cell(1,nstates);

% em for hmm
[mu_area,var_area,ll,mu_area_km,var_area_km]=hmm_multiseq_1d(X,nstates,ptransmat,...
  sexclassifier_params.niters_em,sexclassifier_params.tol_em,...
  [],prior);
if mu_area(1) > mu_area(2),
  mu_area = mu_area(end:-1:1);
  var_area = var_area(end:-1:1);
end
state2sex{argmax(mu_area)} = 'F';
state2sex{argmin(mu_area)} = 'M';

fprintf(1,'mu_area(1) mu_area_km(1) mu_area(2) mu_area_km(2): %0.3f %0.3f %0.3f %0.3f\n',...
  mu_area(1),mu_area_km(1),mu_area(2),mu_area_km(2));
fprintf(1,'vr_area(1) vr_area_km(1) vr_area(2) vr_area_km(2): %0.3f %0.3f %0.3f %0.3f\n',...
  var_area(1),var_area_km(1),var_area(2),var_area_km(2));
  
  
%% save classifier

filterorder = sexclassifierin.areasmooth_filterorder;
maxfreq = sexclassifierin.areasmooth_maxfreq;
maxerrx = sexclassifierin.areasmooth_maxerrx;

if dosave,
  try
    if exist(sexclassifieroutmatfile,'file'),
      delete(sexclassifieroutmatfile);
    end
  catch ,
    % Ignore any errors that occur
  end
  try
    save(sexclassifieroutmatfile,'mu_area','var_area','ptrans','prior','ll',...
      'nstates','state2sex','maxerrx','maxfreq','filterorder','version','analysis_protocol');
  catch ME,
    warning('FlyDiscoClassifySex:save',...
      'Could not save to file %s: %s',sexclassifieroutmatfile,getReport(ME));
    fprintf('Could not save to file %s: %s\n',sexclassifieroutmatfile,getReport(ME));
  end
end

%% classify sex

for fly = 1:numel(trx),
  
  % Viterbi to classify per-frame
  [trx(fly).sex,diagnostics(fly)] = ...
    ClassifySex(X{fly}',mu_area,var_area,ptransmat,state2sex);  %#ok<AGROW>
  
end
