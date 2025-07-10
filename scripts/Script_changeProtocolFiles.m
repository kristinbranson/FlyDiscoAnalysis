%% update protocol files to have ms duration instead of secs
% 20210330, 20210421 updated
modpath;
%% find all the experiment directories collected on or before 20210330
% pull expdirs and load all metdata
rootdatadir = '/groups/branson/bransonlab/flydisco_data';
metadatafile = 'Metadata.xml';

%inputs to getExperimentDirsFlyDisco: 'metadatafile','Metadata.xml','screen_type','*','line_name','*','date','*');
[expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'date','2021*');
idx = [];
% find experiments earlier than 20210331
for i = 1:numel(expdirstruct)
    date_curr = expdirstruct(i).date;
    date_curr = datenum(date_curr,'yyyymmdd');
    idx(i) = date_curr < datenum('20210331','yyyymmdd');
end
idx = logical(idx);

explist = {expdirstruct(idx).file_system_path};

%% make a back file and then edit protocol.mat

protocolfilename = 'protocol.mat';
backupprotocol = 'protocol.mat.bak';
rootdatadir = '/groups/branson/bransonlab/flydisco_data';

% explist = {'VNC_JRC_SS31881_RigB_20210325T144340'};

for i = 1:numel(explist)
    expdir = explist{i};   
    protocolfile = fullfile(rootdatadir,expdir,protocolfilename);
    protcolfilebackup = fullfile(rootdatadir,expdir,backupprotocol);
    if exist(protocolfile,'file')
        cmd = sprintf('cp %s %s \n',protocolfile,protcolfilebackup);
        system(cmd);
        load(protocolfile)
        for j = 1:numel(protocol.duration)
            protocol.duration(j) = protocol.duration(j)*1000;
        end
        save(protocolfile,'protocol');
        sprintf('editted %s \n',protocolfile)
    end
end