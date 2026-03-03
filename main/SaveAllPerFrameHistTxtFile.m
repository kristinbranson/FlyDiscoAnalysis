function SaveAllPerFrameHistTxtFile(histtxtsavename,histperfly,histperexp)

if exist(histtxtsavename,'file'),
  try
    delete(histtxtsavename);
  end
end

% open file
fid = fopen(histtxtsavename,'w');
if fid < 0,
  error('Could not open file %s for writing',histtxtsavename);
end

fns = fieldnames(histperfly);
nfns = numel(fns);
fields = cell(1,nfns);
flyconditions = cell(1,nfns);
frameconditions = cell(1,nfns);
for i = 1:nfns,
  m = regexp(fns{i},'^(?<field>.*)_fly(?<flyconditions>.*)_frame(?<frameconditions>.*)$','names','once');
  if isempty(m),
    error('Could not parse %s into field, fly, frame conditions',fns{i});
  end
  fields{i} = m.field;
  flyconditions{i} = m.flyconditions;
  frameconditions{i} = m.frameconditions;
end

allfields = unique(fields);
nfields = numel(allfields);

for i = 1:nfields,

  % all fns with this field
  field = allfields{i};
  idx = find(strcmp(fields,field));
  nidx = numel(idx);

  %% per-fly hist
    
  % which flies are analyzed per condition
  isdata = cell(1,nidx);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    isdata{jj} = ~isnan(histperfly.(fn).Z);
  end
  
  % conditions concatenated together for this field
  fprintf(fid,'%s_conditions',field);
  for jj = 1:nidx,
    j = idx(jj);
    fprintf(fid,',fly%s_frame%s',flyconditions{j},frameconditions{j});
  end

  %% per-fly hist

  % number of flies analyzed per condition
  fprintf(fid,'\n%s_nflies_analyzed',field);  
  for jj = 1:nidx,
    fprintf(fid,',%d',nnz(isdata{jj}));
  end
  
  % flies analyzed per condition
  fprintf(fid,'\n%s_flies_analyzed',field);
  for jj = 1:nidx,
    fprintf(fid,',%d',find(isdata{jj}));
  end
  
  % linear fractions -- this will be nbins + 2 as it also include
  % fracless, fracmore 
  fprintf(fid,'\n%s_frac_linear_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    frac = cat(1,histperfly.(fn).fracless_linear(isdata{jj}),...
      histperfly.(fn).frac_linear(:,isdata{jj}),...
      histperfly.(fn).fracmore_linear(isdata{jj}));
    fprintf(fid,',%f',frac);
  end
  
  % log fractions
  fprintf(fid,'\n%s_frac_log_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    frac = cat(1,histperfly.(fn).fracless_log(isdata{jj}),...
      histperfly.(fn).frac_log(:,isdata{jj}),...
      histperfly.(fn).fracmore_log(isdata{jj}));
    fprintf(fid,',%f',frac);
  end
  
  fprintf(fid,'\n%s_Z_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%d',histperfly.(fn).Z);
  end
  
  fprintf(fid,'\n%s_fracframesanalyzed_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',histperfly.(fn).fracframesanalyzed);
  end
  fprintf(fid,'\n');

  %% per-exp hist

  fprintf(fid,'%s_mean_frac_linear_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    frac = cat(2,histperexp.(fn).meanfracless_linear,...
      histperexp.(fn).meanfrac_linear,histperexp.(fn).meanfracmore_linear);
    fprintf(fid,',%f',frac);
  end
  
  fprintf(fid,'\n%s_std_frac_linear_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    frac = cat(2,histperexp.(fn).stdfracless_linear,...
      histperexp.(fn).stdfrac_linear,histperexp.(fn).stdfracmore_linear);
    fprintf(fid,',%f',frac);
  end

  fprintf(fid,'\n%s_mean_frac_log_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    frac = cat(2,histperexp.(fn).meanfracless_log,...
      histperexp.(fn).meanfrac_log,histperexp.(fn).meanfracmore_log);
    fprintf(fid,',%f',frac);
  end
  
  fprintf(fid,'\n%s_std_frac_log_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    frac = cat(2,histperexp.(fn).stdfracless_log,...
      histperexp.(fn).stdfrac_log,histperexp.(fn).stdfracmore_log);
    fprintf(fid,',%f',frac);
  end
  
  fprintf(fid,'\n%s_Z_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%d',histperexp.(fn).Z);
  end
  
  fprintf(fid,'\n%s_meanZ_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',histperexp.(fn).meanZ);
  end
  fprintf(fid,'\n');
  
end

fclose(fid);