function [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,varargin)
%rootdatadir - dir to be searched
%datastructname - out data struct
% limit search with screen_type, line or date (can truncate with *)
% assuming standard names
sxclfilestr = 'sexclassifier_diagnostics.txt';
autochcksinfilestr = 'automatic_checks_incoming_results.txt';
autochckscompfilestr = 'automatic_checks_complete_results.txt';
manchcksfilestr = 'manual_fail.txt';


moviefilestr = 'movie.ufmf';
trxfilestr = 'registered_trx.mat';
% TODO 
% add option to not return expdirs without Metadata file
% add option to add autochecks incoming and completed
% add option to check for aborted complete or failed files. 
%%% screen_type from experiment name NOT accurate. change to expname, add
%%% real screen_type filter
[metadatafile, expdirname,line_name,date,nflies,autocheckin,movielength,TrajNum,autocheckcomplete,manualcheck] = myparse(varargin,'metadatafile','Metadata.xml','expdirname','*','line_name','*', ...
    'date','*','nflies',true,'autocheckin',true,'movielength',true,'TrajNum',true,'autocheckcomplete',true,'manualcheck',false);


searchname = sprintf('%s%s%s',date,line_name,expdirname);
%remove extra wildcards
searchname = strrep(searchname,'***','*');
searchname = strrep(searchname,'**','*');

input = fullfile(rootdatadir,searchname);

tmp = dir(input);

expdirstruct = [];
for i = 1:numel(tmp)
    if ismember(tmp(i).name,{'.','..'})
        continue;
    end
    expdir = tmp(i).name;
    metadatafilestr = fullfile(rootdatadir,expdir,metadatafile);
    if exist(metadatafilestr,'file')
        try
    tmpM = ReadMetadataFile(metadatafilestr);
        catch
        fprintf('%s failed ReadMetadataFile\n',expdir)       
        end
    tmpM.file_system_path = fullfile(rootdatadir,expdir);  
    % add date
    daTe = tmpM.exp_datetime(1:8);
    tmpM.date = daTe;    
    tmpM.NoMetadata = false;
    else
        tmpM.file_system_path = expdir;   
        tmpM.NoMetadata = true;
    end
    % add nflies
    if nflies
    sxclfile = fullfile(rootdatadir,expdir,sxclfilestr);
    if exist(sxclfile,'file')
        sex = ReadParams(sxclfile);
        tmpM.nflies = sex.median_nflies;  
    else 
        tmpM.nflies = nan; 
    end
    end
    
    % add autochecks incoming
    if autocheckin
        try
            autochcksinfile = fullfile(rootdatadir,expdir,autochcksinfilestr);
            if exist(autochcksinfile,'file')
                autochcksin = ReadParams(autochcksinfile);
                tmpM.automated_pf_incoming = autochcksin.automated_pf;
                if tmpM.automated_pf_incoming == 'F',
                tmpM.automated_pf_incoming_category = autochcksin.automated_pf_category;
                end
                
            else
                tmpM.automated_pf_incoming = 'NaN';
            end
        catch
            tmpM.automated_pf_incoming = 'could not load file';
        end
    end
    
    % add autochecks complete
    if autocheckcomplete
        try
            autochckscompfile = fullfile(rootdatadir,expdir,autochckscompfilestr);
            if exist(autochckscompfile,'file')
                autochckscomplete = ReadParams(autochckscompfile);
                tmpM.automated_pf = autochckscomplete.automated_pf;
                if tmpM.automated_pf == 'F',
                tmpM.automated_pf_category = autochckscomplete.automated_pf_category;
                tmpM.notes_curation = autochckscomplete.notes_curation;
                end
                
            else
                tmpM.automated_pf = 'NaN';
            end
        catch
            tmpM.automated_pf = 'could not load file';
        end
    end
    % add manual failes
    if manualcheck
        try
            manchcksfile = fullfile(rootdatadir,expdir,manchcksfilestr);
            if exist(manchcksfile,'file')
                tmpM.manual_fail = 'F';
                manchcks = textread(manchcksfile,'%s');               
                tmpM.manual_fail_category = strjoin(manchcks);
            else
                tmpM.manual_fail = 'U';
            end
        catch
            tmpM.manual_fail = 'error check for file';           
        end
    end

    
    
%     % add FlyDiscoAnalysisStatus
%     analysiscompletefile = fullfile(rootdatadir,expdir,analysiscompletefilestr);
%     if FlyDiscoAnalysisStatus
%     if exist(analysiscompletefile,'file')
%         tmpM.FlyDiscoAnalysisStatus = '1';
%     else 
%         tmpM.FlyDiscoAnalysisStatus = '0';
%     end
%     end
    % add nframes
    if movielength
        moviefile = fullfile(rootdatadir,expdir,moviefilestr);
        if movielength
            if exist(moviefile,'file')
                headerinfo = ufmf_read_header(moviefile);
                nframes = headerinfo.nframes;
                tmpM.nframes = nframes;
            else
                tmpM.nframes = nan;
            end
        end
    end
    % add trajectory num
    if TrajNum
        trxfile = fullfile(rootdatadir,expdir,trxfilestr);
        if exist(trxfile,'file')
            load(trxfile,'trx');
            tmpM.trajnum = numel(trx);
            tmpM.pxpermm = trx(1).pxpermm;
        else
            tmpM.trajnum = nan;
            tmpM.pxpermm = nan;
        end
    end
    
    % compile struct 
    expdirstruct = structappend(expdirstruct,tmpM);  
end

    

