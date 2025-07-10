function [success,msg] = SaveStatsCSVFile(stats,i,csvfilename,varargin)

[statfns] = myparse(varargin,'statfns',{});

success = false;
msg = '';

if isempty(statfns),
  statfns = fieldnames(stats.means);
end

fid = fopen(csvfilename,'w');
if fid < 1,
  msg = sprintf('Could not open file %s for writing',csvfilename);
  return;
end

ispvalue = isfield(stats,'pvalue_bigger') && isfield(stats,'pvalue_smaller');
iszscore = isfield(stats,'zscores');

fprintf(fid,'field,value');
if iszscore,
  fprintf(fid,',zscore');
end
if ispvalue,
  fprintf(fid,',pvalue_bigger,pvalue_smaller');
end
fprintf(fid,'\n');

for j = 1:numel(statfns),
  
  fn = statfns{j};
  fprintf(fid,'%s,%e',fn,stats.means.(fn)(i));

  if iszscore,
    fprintf(fid,',%e',stats.zscores.(fn)(i));
  end
  
  if ispvalue,
    fprintf(fid,',%e,%e',stats.pvalue_bigger.(fn)(i),stats.pvalue_smaller.(fn)(i));
  end

  fprintf(fid,'\n');
  
end

fclose(fid);
success = true;