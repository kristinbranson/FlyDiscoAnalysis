function [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,varargin)
%rootdatadir - dir to be searched
%datastructname - out data struct
% limit search with screen_type, line or date (can truncate with *)
% assuming standard names
sxclfilestr = 'sexclassifier_diagnostics.txt';
autochcksinfilestr = 'automatic_checks_incoming_results.txt';
analysiscompletefilestr = 'ANALYSIS-COMPLETED-SUCCESSFULLY';
moviefilestr = 'movie.ufmf';
% TODO 
% add option to not return expdirs without Metadata file
% add option to add autochecks incoming and completed
% add option to check for aborted complete or failed files. 
[metadatafile, screen_type,line_name,date,nflies,autocheckin,FlyDiscoAnalysisStatus,movielength] = myparse(varargin,'metadatafile','Metadata.xml','screen_type','*','line_name','*', ...
    'date','*','nflies',false,'autocheckin',false,'FlyDiscoAnalysisStatus', false,'movielength',false);


searchname = sprintf('%s%s%s',screen_type,line_name,date);
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
    tmpM = ReadMetadataFile(metadatafilestr);
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
    
    % add FlyDiscoAnalysisStatus
    analysiscompletefile = fullfile(rootdatadir,expdir,analysiscompletefilestr);
    if FlyDiscoAnalysisStatus
    if exist(analysiscompletefile,'file')
        tmpM.FlyDiscoAnalysisStatus = '1';
    else 
        tmpM.FlyDiscoAnalysisStatus = '0';
    end
    end
    % add nframes
    moviefile = fullfile(rootdatadir,expdir,moviefilestr);
    if movielength
        if exist(moviefile,'file')
            [~,nframes] = get_readframe_fcn(moviefile);
            tmpM.nframes = nframes;
        else
            tmpM.nframes = nan;
        end
    end
    
    % compile struct 
    expdirstruct = structappend(expdirstruct,tmpM);  
end

    

