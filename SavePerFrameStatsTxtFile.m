function SavePerFrameStatsTxtFile(statstxtsavename,fn,statsperfly,statsperexp)
  
% open file
fid = fopen(statstxtsavename,'a');
if fid < 0,
  error('Could not open file %s for appending',statstxtsavename);
end

%nprctiles = numel(stats_params.prctiles_compute);

% per-fly stats
isdata = ~isnan(statsperfly.Z);
fprintf(fid,'%s_nflies_analyzed,%d\n',fn,nnz(isdata));
fprintf(fid,'%s_flies_analyzed',fn);
fprintf(fid,',%d',find(isdata));
fprintf(fid,'\n%s_mean_perfly',fn);
fprintf(fid,',%f',statsperfly.mean(isdata));
fprintf(fid,'\n%s_std_perfly',fn);
fprintf(fid,',%f',statsperfly.std(isdata));
fprintf(fid,'\n%s_Z_perfly',fn);
fprintf(fid,',%d',statsperfly.Z(isdata));
fprintf(fid,'\n%s_fracframesanalyzed_perfly',fn);
fprintf(fid,',%f',statsperfly.fracframesanalyzed(isdata));
fprintf(fid,'\n%s_prctiles_perfly',fn);
fprintf(fid,',%f',statsperfly.prctiles(:,isdata));
% for j = 1:nprctiles,
%    if round(stats_params.prctiles_compute(j)) == stats_params.prctiles_compute(j),
%     s = sprintf('%03d',stats_params.prctiles_compute(j));
%   else
%     s = sprintf('%03g',stats_params.prctiles_compute(j));
%   end
%   fprintf(fid,'\n%s_prctile_%s_perfly',fn,s);
%   fprintf(fid,',%f',statsperfly.prctiles(j,isdata));
% end
fprintf(fid,'\n');

% per-exp stats
fprintf(fid,'%s_meanmean_perexp,%f\n',fn,statsperexp.meanmean);
fprintf(fid,'%s_stdmean_perexp,%f\n',fn,statsperexp.stdmean);
fprintf(fid,'%s_meanstd_perexp,%f\n',fn,statsperexp.meanstd);
fprintf(fid,'%s_Z_perexp,%d\n',fn,statsperexp.Z);
fprintf(fid,'%s_meanZ_perexp,%f\n',fn,statsperexp.meanZ);
fprintf(fid,'%s_meanprctiles_perexp',fn);
fprintf(fid,',%f',statsperexp.meanprctiles);
fprintf(fid,'\n%s_stdprctiles_perexp',fn);
fprintf(fid,',%f',statsperexp.stdprctiles);
% for j = 1:nprctiles,
%   if round(stats_params.prctiles_compute(j)) == stats_params.prctiles_compute(j),
%     s = sprintf('%03d',stats_params.prctiles_compute(j));
%   else
%     s = sprintf('%03g',stats_params.prctiles_compute(j));
%   end
%   fprintf(fid,'%s_meanprctile_%s_perexp,%f\n',fn,s,statsperexp.meanprctiles(j));
%   fprintf(fid,'%s_stdprctile_%s_perexp,%f\n',fn,s,statsperexp.stdprctiles(j));
% end


fclose(fid);