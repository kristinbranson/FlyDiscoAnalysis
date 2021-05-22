% running ReSaveMetadata

logfiledir = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/MetadataFixes';
logfilename = 'expdirs_wk5_metachanges_logofauto.csv';

logfile = fullfile(logfiledir,logfilename);

%% per experiment
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210506_testingMetadataMods/VNC_JRC_SS52985_RigA_20210503T152107';
metadatafilename = 'Metadata.xml';
resavefilename = fullfile(expdir,metadatafilename);
metadata = ReadMetadataFile(resavefilename);

% save changes to logfile 
fid = fopen(logfile,'a');
 
%% manual changes create struct

changestruct = struct;
changestruct.experimenter = 'ALiceRobie';
changestruct.hours_sorted = 1400;

%% replace metadata fields with changes
fprintf(fid,'%s, ',expdir)
metadatafieldsTOBECHANGED = fieldnames(changestruct);
metadatafields = fieldnames(metadata);
for i = 1:numel(metadatafieldsTOBECHANGED)
    if isnumeric(changestruct.(metadatafieldsTOBECHANGED{i})) && isnumeric(metadata.(metadatafieldsTOBECHANGED{i}))
    fprintf(fid,'%s,%d,%d,',(metadatafieldsTOBECHANGED{i}),metadata.(metadatafieldsTOBECHANGED{i}),changestruct.(metadatafieldsTOBECHANGED{i}));
    metadata.(metadatafieldsTOBECHANGED{i}) = changestruct.(metadatafieldsTOBECHANGED{i}); 
    else 
        fprintf(fid,'%s,%s,%s,',(metadatafieldsTOBECHANGED{i}),metadata.(metadatafieldsTOBECHANGED{i}),changestruct.(metadatafieldsTOBECHANGED{i}));
    metadata.(metadatafieldsTOBECHANGED{i}) = changestruct.(metadatafieldsTOBECHANGED{i}); 
    end
end
fprintf(fid,'\n');
fclose(fid);

%% resave metadata and backup


[success] = ResaveMetadata(metadata,resavefilename);

