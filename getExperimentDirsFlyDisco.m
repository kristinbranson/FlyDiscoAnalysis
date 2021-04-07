function [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,varargin)
%rootdatadir - dir to be searched
%datastructname - out data struct
% limit search with screen_type, line or date (can truncate with *)


% TODO 
% add option to not return expdirs without Metadata file
% add option to add autochecks incoming and completed
% add option to check for aborted complete or failed files. 
[metadatafile, screen_type,line_name,date] = myparse(varargin,'metadatafile','Metadata.xml','screen_type','*','line_name','*','date','*');

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
    tmpM.file_system_path = expdir;  
    % add date
    daTe = tmpM.exp_datetime(1:8);
    tmpM.date = daTe;    
    tmpM.NoMetadata = false;
    else
        tmpM.file_system_path = expdir;   
        tmpM.NoMetadata = true;
    end
    
    expdirstruct = structappend(expdirstruct,tmpM);   
    
end

    

