function SavePerFrameHistTxtFile(histtxtsavename,fn,histperfly,histperexp)

% open file
fid = fopen(histtxtsavename,'a');
if fid < 0,
  error('Could not open file %s for appending',histtxtsavename);
end

% per-fly hist
isdata = ~isnan(histperfly.Z);
fprintf(fid,'%s_nflies_analyzed,%d\n',fn,nnz(isdata));
fprintf(fid,'%s_flies_analyzed',fn);
fprintf(fid,',%d',find(isdata));
fprintf(fid,'\n%s_frac_linear_perfly',fn);
frac = cat(1,histperfly.fracless_linear(isdata),...
  histperfly.frac_linear(:,isdata),...
  histperfly.fracmore_linear(isdata));
fprintf(fid,',%f',frac);
fprintf(fid,'\n%s_frac_log_perfly',fn);
frac = cat(1,histperfly.fracless_log(isdata),...
  histperfly.frac_log(:,isdata),...
  histperfly.fracmore_log(isdata));
fprintf(fid,',%f',frac);
fprintf(fid,'\n%s_Z_perfly',fn);
fprintf(fid,',%d',histperfly.Z);
fprintf(fid,'\n%s_fracframesanalyzed_perfly',fn);
fprintf(fid,',%f',histperfly.fracframesanalyzed);
fprintf(fid,'\n');

% per-exp hist
frac = cat(2,histperexp.meanfracless_linear,...
  histperexp.meanfrac_linear,histperexp.meanfracmore_linear);
fprintf(fid,'%s_mean_frac_linear_perexp',fn);
fprintf(fid,',%f',frac);
frac = cat(2,histperexp.stdfracless_linear,...
  histperexp.stdfrac_linear,histperexp.stdfracmore_linear);
fprintf(fid,'\n%s_std_frac_linear_perexp',fn);
fprintf(fid,',%f',frac);

frac = cat(2,histperexp.meanfracless_log,...
  histperexp.meanfrac_log,histperexp.meanfracmore_log);
fprintf(fid,'%s_mean_frac_log_perexp',fn);
fprintf(fid,',%f',frac);
frac = cat(2,histperexp.stdfracless_log,...
  histperexp.stdfrac_log,histperexp.stdfracmore_log);
fprintf(fid,'\n%s_std_frac_log_perexp',fn);
fprintf(fid,',%f',frac);

fprintf(fid,'\n%s_Z_perexp,%d\n',fn,histperexp.Z);
fprintf(fid,'%s_meanZ_perexp,%f\n',fn,histperexp.meanZ);

fclose(fid);