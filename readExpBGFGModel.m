function model = readExpBGFGModel(matfilename)

if ~exist(matfilename,'file'),
  error('ExpBGFGModel %s does not exist',matfilename);
end

model = load(matfilename);
fns = {'fg_mu','fg_sigma','bg_mu','bg_sigma','obj_detection_dist_frac_fg',...
  'obj_detection_dist_frac_bg','always_bg_mask_frac','always_bg_mask'};
missing = ~isfield(model,fns);
if any(missing),
  error(['Missing fields of ExpBGFGModel:',sprintf(' %s',fns{missing})]);
end