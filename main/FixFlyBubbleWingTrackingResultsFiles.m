function [errmsg,sucmsg] = FixFlyBubbleWingTrackingResultsFiles(expdir)

errmsg = [];
sucmsg= [];

% hardcode names maybe bad idea
%wing tracking trx filename
wingtrxfilestr = 'wingtracking_results.mat';

%new name for wing tracking perframe file 
wingperframefilestr = 'wingtracking_perframedata.mat';

%
trxfilestr = 'registered_trx.mat';

% check if trxfilestr exists
if ~exist(fullfile(expdir, trxfilestr),'file')
        errmsg = sprintf('%s, registered trx does not exist \n',expdir)
        return
end 

% check if wingtracking has run
if ~exist(fullfile(expdir, wingtrxfilestr),'file')
        errmsg = sprintf('%s, wing tracking results do not exist \n',expdir)
        return
end 

% check that wingtrxfilestr doesn't contain trx
wingresultsVariables = who('-file',fullfile(expdir, wingtrxfilestr));
if ~ismember('trx',wingresultsVariables)
    % mv wingtrxfilestr to wingperframefilestr
    movefile(fullfile(expdir,wingtrxfilestr),fullfile(expdir,wingperframefilestr));
else
    sucmsg = sprintf('%s, wingresults contains trx\n',expdir)
    return
end


% check trxfilestr for wing tracking fields
load(fullfile(expdir, trxfilestr),'trx')
if isfield(trx,'wing_anglel')
    % make relative softlink wingtrxfilestr to trxfilestr
    pwdprev = pwd;
    cd(expdir)
    cmd = sprintf('ln -s ./%s ./%s', trxfilestr,wingtrxfilestr);
    system(cmd);
    cd(pwdprev);
    sucmsg = sprintf('%s, created wing tracking trx\n',expdir)
else
    errmsg = sprintf('%s, registered_trx does not contain wing tracking data\n',expdir)
    return
end 

