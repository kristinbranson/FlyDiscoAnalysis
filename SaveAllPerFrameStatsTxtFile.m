function SaveAllPerFrameStatsTxtFile(statstxtsavename,statsperfly,statsperexp)
  
% open file
fid = fopen(statstxtsavename,'w');
if fid < 0,
  error('Could not open file %s for writing',statstxtsavename);
end

fns = fieldnames(statsperfly);
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

  %% per-fly stats
    
  % which flies are analyzed per condition
  isdata = cell(1,nidx);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    isdata{jj} = ~isnan(statsperfly.(fn).Z);
  end
  
  % conditions concatenated together for this field
  fprintf(fid,'%s_conditions',field);
  for jj = 1:nidx,
    j = idx(jj);
    fprintf(fid,',fly%s_frame%s',flyconditions{j},frameconditions{j});
  end
  
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
  
  % mean per fly
  fprintf(fid,'\n%s_mean_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperfly.(fn).mean(isdata{jj}));
  end
  
  % std per fly
  fprintf(fid,'\n%s_std_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperfly.(fn).std(isdata{jj}));
  end
  
  % number of frames per fly
  fprintf(fid,'\n%s_Z_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%d',statsperfly.(fn).Z(isdata{jj}));
  end

  % fraction of trajectory per fly
  fprintf(fid,'\n%s_fracframesanalyzed_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperfly.(fn).fracframesanalyzed(isdata{jj}));
  end
  
  % percentiles per fly
  fprintf(fid,'\n%s_prctiles_perfly',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperfly.(fn).prctiles(:,isdata{jj}));
  end
  fprintf(fid,'\n');
  
  %% per-exp stats

  fprintf(fid,'%s_meanmean_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperexp.(fn).meanmean);
  end
  
  fprintf(fid,'\n%s_stdmean_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperexp.(fn).stdmean);
  end
  
  fprintf(fid,'\n%s_meanstd_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperexp.(fn).meanstd);
  end
  
  fprintf(fid,'\n%s_Z_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperexp.(fn).Z);
  end
  
  fprintf(fid,'\n%s_meanZ_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperexp.(fn).meanZ);
  end
  
  fprintf(fid,'\n%s_meanprctiles_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperexp.(fn).meanprctiles);
  end
  
  fprintf(fid,'\n%s_stdprctiles_perexp',field);
  for jj = 1:nidx,
    j = idx(jj);
    fn = fns{j};
    fprintf(fid,',%f',statsperexp.(fn).stdprctiles);
  end
  fprintf(fid,'\n');

end

fclose(fid);