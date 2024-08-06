function jobid = parse_bsub_jobid(result)

jobid = regexp(result,'Job <(\d+)>','once','tokens');
assert(~isempty(jobid));
jobid = str2double(jobid{1});