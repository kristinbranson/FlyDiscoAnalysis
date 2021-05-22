% running ReSaveMetadata

logfiledir = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/MetadataFixes';
logfilename = 'expdirs_wk5_metachanges_logofautoTEMP.csv';

logfile = fullfile(logfiledir,logfilename);

% save changes to logfile
fid2 = fopen(logfile,'a');

%% explist 
rootdatadir = '/groups/branson/bransonlab/flydisco_data';
explist  = {'VNC_EXT_VGLUT-GAL4_RigD_20210420T142706', ...
'VNC_JRC_SS50975_RigD_20210420T130327', ...
'VNC_JRC_SS51901_RigD_20210420T131534', ...
'VNC_JRC_SS51902_RigD_20210420T153114', ...
'VNC_JRC_SS59198_RigD_20210420T133940', ...
'VNC_JRC_SS60231_RigD_20210420T135100', ...
'VNC_JRC_SS60234_RigD_20210420T140255', ...
'VNC_JRC_SS60240_RigD_20210420T141455', ...
'VNC_JRC_SS60246_RigD_20210420T154237', ...
'VNC_JRC_SS60632_RigD_20210420T144027', ...
'VNC_JRC_SS61066_RigD_20210420T145159', ...
'VNC_JRC_SS64190_RigD_20210420T150349', ...
'VNC_JRC_SS68333_RigD_20210420T151651', ...
'VNC_YNA_K_162984_RigD_20210420T132759'};
%testing - accidently did in reall data direct DATA NOT PROTECTED PROPERLY
% explist = {'VNC_JHS_K_85321_RigA_20210322T164447'};

    

for j = 1:numel(explist)
    % per experiment
    expdir = fullfile(rootdatadir,explist{j})
    metadatafilename = 'Metadata.xml';
    resavefilename = fullfile(expdir,metadatafilename);
    metadata = ReadMetadataFile(resavefilename);
    

    % manual changes create struct    
    changestruct = struct;
    changestruct.experimenter = 'managanc';

    
    % replace metadata fields with changes
    fprintf(fid2,'%s, ',expdir);
    metadatafieldsTOBECHANGED = fieldnames(changestruct);
    metadatafields = fieldnames(metadata);
    for i = 1:numel(metadatafieldsTOBECHANGED)
        if isnumeric(changestruct.(metadatafieldsTOBECHANGED{i})) && isnumeric(metadata.(metadatafieldsTOBECHANGED{i}))
            fprintf(fid2,'%s,%d,%d,',(metadatafieldsTOBECHANGED{i}),metadata.(metadatafieldsTOBECHANGED{i}),changestruct.(metadatafieldsTOBECHANGED{i}));
            metadata.(metadatafieldsTOBECHANGED{i}) = changestruct.(metadatafieldsTOBECHANGED{i});
        else
            fprintf(fid2,'%s,%s,%s,',(metadatafieldsTOBECHANGED{i}),metadata.(metadatafieldsTOBECHANGED{i}),changestruct.(metadatafieldsTOBECHANGED{i}));
            metadata.(metadatafieldsTOBECHANGED{i}) = changestruct.(metadatafieldsTOBECHANGED{i});
        end
    end
    fprintf(fid2,'\n');
    
    % resave metadata and backup   
    [success] = ResaveMetadata(metadata,resavefilename);
    
end

fclose(fid2);
