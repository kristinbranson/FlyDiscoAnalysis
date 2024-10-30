function [trx, diagnostics, mu_area, var_area, state2sex, ll] = add_predetermined_sex_field_to_trx(trx, metadata, X)

fprintf('Metadata-specified gender is ''%s'', so not doing sex classification, just setting sex to ''%s'' for all flies\n', metadata.gender, upper(metadata.gender));

% set sex to metadata.gender for all flies
% also set diagnostics that we can
mean_area_all = nanmean(cat(1,X{:})); %#ok<NANMEAN> 
var_area_all = nanstd(cat(1,X{:}),1); %#ok<NANSTD> 
for fly = 1:numel(trx),
  trx(fly).sex = repmat({upper(metadata.gender)},[1,trx(fly).nframes]);
  diagnostics_curr = struct();
  diagnostics_curr.normhmmscore = nan;
  diagnostics_curr.nswaps = 0;
  diagnostics_curr.meanabsdev = nanmean(abs(X{fly}-mean_area_all)); %#ok<NANMEAN> 
  diagnostics(fly) = diagnostics_curr; %#ok<AGROW>
end

mu_area = nan(1,2);
var_area = nan(1,2);
if strcmpi(metadata.gender,'f'),
  mu_area(2) = mean_area_all;
  var_area(2) = var_area_all;
elseif strcmpi(metadata.gender,'m'),
  mu_area(1) = mean_area_all;
  var_area(1) = var_area_all;
end
state2sex = {'M','F'};
ll = [];
