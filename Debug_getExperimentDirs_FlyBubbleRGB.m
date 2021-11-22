%Deubg get all experiment directories 
%% setup path
% addpath /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/JAABA/filehandling
% addpath /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/JAABA/misc
modpath
% metadatafile = 'Metadata.xml';
% %% pull expdirs and load all metdata
rootdatadir = '/groups/branson/bransonlab/flydisco_data';

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_wk18wk19wk20_metadatachanges';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC','date','202110*','autocheckin',true,'autocheckcomplete',true);

% rootdatadir = '/groups/branson/bransonlab/from_tier2/fly_bubble/bubble_data';
% savefile = '/groups/branson/home/robiea/Projects_data/FlyBubble/metadata_ANY_bubble_data_20211008';

savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC_20211104';
[expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC','autocheckin',true,'autocheckcomplete',true);


% %inputs to getExperimentDirsFlyDisco: 'metadatafile','Metadata.xml','expdirname','*','line_name','*', ...
%     'date','*','nflies',false,'autocheckin',false,'movielength',false,'TrajNum',false,'autocheckcomplete',false);

% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,varagin);



% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC','autocheckin',true,'autocheckcomplete',true,'TrajNum',true);
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','autocheckin',true);


%% select on expdirs with certain led protocol
idx = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
expdirstruct= expdirstruct(idx);


%% make tsv with ALL metadata fields
savefilecsv = [savefile,'.csv'];
expdirtable = struct2table(expdirstruct);
writetable(expdirtable,savefilecsv,'Delimiter','tab');

%% save mat file
metadata = expdirstruct;
save(savefile,'metadata')

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
%% make csv file for all experiments (use for pulling data for metadata changes) with autochecks
fid = fopen([savefile,'.tsv'],'w');

fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','linename','experimentor','notes_tech','notes_behav','automated_pf');

for i = 1:numel(expdirstruct)  
  [~,expname] = fileparts(expdirstruct(i).file_system_path);
  datestr = expdirstruct(i).date;
  datetimestr = expdirstruct(i).exp_datetime;
  linename = expdirstruct(i).line;
  experimenter = expdirstruct(i).experimenter;
  notestech = expdirstruct(i).notes_technical;
  notesbeh = expdirstruct(i).notes_behavioral;
  autopf = expdirstruct(i).automated_pf;


%   notes = experiments_eddison(i).notes_curation;
  fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,linename,experimenter,notestech,notesbeh,autopf);
end

fclose(fid);
