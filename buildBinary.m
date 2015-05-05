function buildBinary(projname)

projfile = [projname '.prj'];
assert(exist(projfile,'file')==2,'Can''t find file ''%s''.',projfile);

% take snapshot
REPO_SNAPSHOT_SCRIPT = 'repo_snapshot.sh';
ssfile = which(REPO_SNAPSHOT_SCRIPT);
assert(~isempty(ssfile),'Can''t find ''%s''.',REPO_SNAPSHOT_SCRIPT);
cmd = sprintf('%s > repo_snapshot_%s.txt',ssfile,datestr(now,'yyyymmddTHHMMSS'));
system(cmd);

% build
deploytool('-build',projfile); % need to figure out how to record build log
% dryname = sprintf('dry_bld_%s_%s.txt',projname,datestr(now,'yyyymmddTHHMMSS'));
% diary(dryname);
% mkdir(projname);
% mkdir(projname,'src');
% mcc -o FlyBowlAutomaticChecks_Complete -W main:FlyBowlAutomaticChecks_Complete -T link:exe -d /groups/flyprojects/home/leea30/flybowl/FBA/FlyBowlAnalysis/FlyBowlAutomaticChecks_Complete/src -w enable:specified_file_mismatch -w enable:repeated_file -w enable:switch_ignored -w enable:missing_lib_sentinel -w enable:demo_license -v /groups/flyprojects/home/leea30/flybowl/FBA/FlyBowlAnalysis/FlyBowlAutomaticChecks_Complete.m -a /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling/cleanup_ctrax_data.m -a /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling/load_tracks.m -a /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe/Macguffin.m -a /groups/branson/home/bransonk/tracking/code/JCtrax/misc/myparse.m -a /groups/branson/home/bransonk/tracking/code/JCtrax/misc/myparse_nocheck.m -a /groups/flyprojects/home/leea30/flybowl/FBA/FlyBowlAnalysis/ReadMetadataFile.m -a /groups/flyprojects/home/leea30/flybowl/FBA/FlyBowlAnalysis/ReadParams.m;
% diary off;

% post
% move build outputs, snapshot, diary into build dirs
% change permissions for binaries


