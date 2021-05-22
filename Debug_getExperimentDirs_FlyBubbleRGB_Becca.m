%Deubg get all experiment directories 
%% setup path
addpath /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/JAABA/filehandling
addpath /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/JAABA/misc
%% pull expdirs and load all metdata
rootdatadir = '/groups/branson/bransonlab/flydisco_data';

metadatafile = 'Metadata.xml';

% %inputs to getExperimentDirsFlyDisco: 'metadatafile','Metadata.xml','screen_type','*','line_name','*','date','*');
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'screen_type','VNC*','date','20210405*');

% updated inputs 'metadatafile','Metadata.xml','screen_type','*','line_name','*', ...
%     'date','*','nflies',false,'autocheckin',false,'FlyDiscoAnalysisStatus', false);
[expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'screen_type','VNC*','date','2021042*',...
    'nflies',true,'autocheckin',true,'FlyDiscoAnalysisStatus', true);

savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_2021042_metadatachanges.csv';

%% make tsv with ALL metadata fields

expdirtable = struct2table(expdirstruct);
writetable(expdirtable,savefile,'Delimiter','tab');

%% make csv file for all experiments (use for pulling data for metadata changes) 
fid = fopen([savefile,'.tsv'],'w');

fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','linename','experimentor','notes_tech','notes_behav','flag_redo');

for i = 1:numel(expdirstruct)
  [~,expname] = fileparts(expdirstruct(i).file_system_path);
  datestr = expdirstruct(i).date;
  datetimestr = expdirstruct(i).exp_datetime;
  linename = expdirstruct(i).line;
  experimenter = expdirstruct(i).experimenter;
  notestech = expdirstruct(i).notes_technical;
  notesbeh = expdirstruct(i).notes_behavioral;
  redoflag = expdirstruct(i).flag_redo;


%   notes = experiments_eddison(i).notes_curation;
  fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,linename,experimenter,notestech,notesbeh,string(redoflag));
end

fclose(fid);
