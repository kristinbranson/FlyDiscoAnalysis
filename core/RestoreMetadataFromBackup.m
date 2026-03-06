function [success] = RestoreMetadataFromBackup(metadatafile,metadatabackupfile)
success = false;

[a,b,c] = fileparts(metadatafile);
TODELmetadatafile = fullfile(a,['TODEL_',b,c]);

if exist(TODELmetadatafile,'file')
    delete(TODELmetadatafile)
end

% if existing metadata backup to tmp file
if exist(metadatafile,'file')
    copyfile(metadatafile,TODELmetadatafile);
    if ~exist(TODELmetadatafile,'file')
        success = false;
    end
end

if exist(metadatabackupfile,'file')
    copyfile(metadatabackupfile, metadatafile,'f')
    if exist(metadatafile,'file')
        success = true;
    end
end
