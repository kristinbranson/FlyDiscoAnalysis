%Deubg get all experiment directories 
%% setup path
% addpath /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/JAABA/filehandling
% addpath /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/JAABA/misc
modpath
% metadatafile = 'Metadata.xml';
% %% pull expdirs and load all metdata
rootdatadir = '/groups/branson/bransonlab/flydisco_data';

% %inputs to getExperimentDirsFlyDisco: 'metadatafile','Metadata.xml','expdirname','*','line_name','*', ...
%     'date','*','nflies',false,'autocheckin',false,'movielength',false,'TrajNum',false,'autocheckcomplete',false);

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_allflydisco_20230831';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','movielength',false,'TrajNum',false,'autocheckin',false,'autocheckcomplete',true);


% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_allflydisco_20230831';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','movielength',false,'TrajNum',false,'autocheckin',false,'autocheckcomplete',true,'manualcheck',true);

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_LPC1';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','LPC1','date','_2023*','movielength',false,'TrajNum',false,'autocheckin',true,'autocheckcomplete',true,'manualcheck',true);


% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC_week1';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC','date','202103*','autocheckin',true,'autocheckcomplete',true);


% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_AmpRec';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','AmpRec','date','202208*','autocheckin',true,'autocheckcomplete',true);

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNCall_20240313';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC*','autocheckin',true,'autocheckcomplete',true,'manualcheck',true);


% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC2_2023';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC2*','date','2023*','autocheckin',true,'autocheckcomplete',true,'manualcheck',true,'TrajNum',true);


% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC3_julytemperatureproblem';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC3*','date','202407*','autocheckin',true,'autocheckcomplete',true,'manualcheck',false,'TrajNum',false);


% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC3_2024_traj';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC3*','date','2024*','autocheckin',true,'autocheckcomplete',true,'manualcheck',true,'TrajNum',true);

savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_2024_allVNC';
[expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC*','autocheckin',true,'autocheckcomplete',true,'manualcheck',true,'TrajNum',false);


% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC_20240711';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC*','autocheckin',true,'autocheckcomplete',true,'manualcheck',true,'TrajNum',false);

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC2_20240711';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC2*','autocheckin',true,'autocheckcomplete',true,'manualcheck',false,'TrajNum',false);

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC2202305_20240822';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml'car,'expdirname','VNC2*','date','202305*','autocheckin',true,'autocheckcomplete',true,'manualcheck',false,'TrajNum',false);

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC2202407_20240823';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC3*','date','202407*','autocheckin',true,'autocheckcomplete',true,'manualcheck',false,'TrajNum',false);

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_allVNC_tofindmanualF_20240824';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC*','autocheckin',true,'autocheckcomplete',true,'manualcheck',true,'TrajNum',false);


% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC2_20220915';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC*','autocheckin',true,'autocheckcomplete',true);

% rootdatadir = '/groups/branson/bransonlab/from_tier2/fly_bubble/bubble_data';
% savefile = '/groups/branson/home/robiea/Projects_data/FlyBubble/metadata_ANY_bubble_data_20211008';

% savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/';
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC','autocheckin',true,'autocheckcomplete',true);



% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,varagin);



% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC','autocheckin',true,'autocheckcomplete',true,'TrajNum',true);
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','autocheckin',true);
% save raw pull data
nowdatetime = datestr(now,'yyyymmddTHHMMSS');
save([savefile,'_',nowdatetime,'.mat'],'expdirstruct');
%% rerun to add autochecks


%% select on expdirs with certain led protocol
% idx = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
idx = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC2');
expdirstruct = expdirstruct(idx);

%% both screentypes for VNC
% expdirstruct = metadata;
idx1 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
idx2 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC2');
% idx3 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_led5secVNC');
% idxall = idx1 + idx2 + idx3;
idxall = idx1 +idx2;
expdirstruct2 = expdirstruct(logical(idxall));

%% find experiments more recent than target date
startdate = datenum('20210328T000000','yyyymmddTHHMMSS');
dates = datenum({expdirstruct.exp_datetime},'yyyymmddTHHMMSS');
idxd = find(dates > startdate);
expdirstruct = expdirstruct(idxd);
%% find experiments older than target date
startdate = datenum('20210328T000000','yyyymmddTHHMMSS');
dates = datenum({expdirstruct.exp_datetime},'yyyymmddTHHMMSS');
idxd = find(dates < startdate);
expdirstruct = expdirstruct(idxd);
%% make tsv with ALL metadata fields
savefilecsv = [savefile,'.csv'];
expdirtable = struct2table(expdirstruct);
writetable(expdirtable,savefilecsv,'Delimiter','tab');

%% save mat file
metadata = expdirstruct2;
save([savefile,'_allVNC.mat'],'metadata')

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
fid = fopen([savefile,'_metadata','.tsv'],'w');

fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','linename','experimentor','gender','notes_tech','notes_behav','automated_pf');

for i = 1:numel(expdirstruct)  
  [~,expname] = fileparts(expdirstruct(i).file_system_path);
  datestr = expdirstruct(i).date;
  datetimestr = expdirstruct(i).exp_datetime;
  linename = expdirstruct(i).line;
  experimenter = expdirstruct(i).experimenter;
  gender = expdirstruct(i).gender;
  notestech = expdirstruct(i).notes_technical;
  notesbeh = expdirstruct(i).notes_behavioral;
  autopf = expdirstruct(i).automated_pf;


%   notes = experiments_eddison(i).notes_curation;
  fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,linename,experimenter,gender,notestech,notesbeh,autopf);
end

fclose(fid);

%% make csv file for all experiments (use for pulling data for metadata changes) with autochecks and manual pf
fid = fopen([savefile,'_metadata','.tsv'],'w');

fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','linename','experimentor','gender','notes_tech','notes_behav','automated_pf','manual_f');

for i = 1:numel(expdirstruct)  
  [~,expname] = fileparts(expdirstruct(i).file_system_path);
  datestr = expdirstruct(i).date;
  datetimestr = expdirstruct(i).exp_datetime;
  linename = expdirstruct(i).line;
  experimenter = expdirstruct(i).experimenter;
  gender = expdirstruct(i).gender;
  notestech = expdirstruct(i).notes_technical;
  notesbeh = expdirstruct(i).notes_behavioral;
  autopf = expdirstruct(i).automated_pf;
  manf = expdirstruct(i).manual_fail;

%   notes = experiments_eddison(i).notes_curation;
  fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,linename,experimenter,gender,notestech,notesbeh,autopf,manf);
end

fclose(fid);
%% make experiment list for VNC and VNC2 with only autochecks = p  
expdirstructOG = expdirstruct;

% screen_type VNC2
idx = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC2');
expdirstruct = expdirstruct(idx);

pidx = find([expdirstruct.automated_pf] == 'P');

explist = {expdirstruct(pidx).file_system_path};
save([savefile,'_VNC2passAuto.mat'],'explist')

% screen_type VNC
expdirstruct = expdirstructOG;
idx = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
expdirstruct = expdirstruct(idx);

pidx = find([expdirstruct.automated_pf] == 'P');

explist = {expdirstruct(pidx).file_system_path};
save([savefile,'_VNCpassAuto.mat'],'explist')

% all list
expdirstruct = expdirstructOG;
idx1 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
idx2 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC2');
idxall = idx1 + idx2;
expdirstruct = expdirstruct(logical(idxall));
pidx = find([expdirstruct.automated_pf] == 'P');

explist = {expdirstruct(pidx).file_system_path};
save([savefile,'_VNCallpassAuto.mat'],'explist');

savefiletext = [savefile,'_VNCallpass.txt'];
fid = fopen(savefiletext,'w');
for i = 1:numel(explist)
    fprintf(fid,'%s\n',explist{i});
end
fclose(fid);

%% check for manual fails in expdirs with autochecks = P
idx = true(1,numel(explist));
for i = 1:numel(explist)
    expdir = explist{i};
    if exist(fullfile(expdir,'manual_fail.txt'),'file')
        idx(i) = false;
    end
end

explist = explist(idx);

savefiletext = [savefile,'_VNCall_passAandM.mat'];
save(savefiletext,'explist');

%% make metadata struct with only auto = P and manual = P 
% 20230130
% start with just VNC and VNC2 datastruct
% auto_pf = P
% pidx = find([expdirstruct2.automated_pf] == 'P');
% pidx = find([expdirstruct.automated_pf] == 'P');
pidx = strcmp({expdirstruct.automated_pf},'P');
expdirstruct3 = expdirstruct(pidx);
% no manual_fail.txt file
explist = {expdirstruct3.file_system_path};
idx = true(1,numel(explist));
for i = 1:numel(explist)
    expdir = explist{i};
    if exist(fullfile(expdir,'manual_fail.txt'),'file')
        idx(i) = false;
    end
end

metadata = expdirstruct3(idx);
savefiletext = [savefile,'_VNC_passAandM_metadata.mat'];
save(savefiletext,'metadata');

%% filter VNC2 and VNC3 out of VNC data
savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC_20240711';
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC_20240711_20240711T175926.mat
idx = strcmp({metadata.screen_type},'non_olympiad_dickson_VNC');
expdirstruct = expdirstruct(idx);
%input into cell above 

%% fix problem with unique and expdirstruct

out = unique(cellfun(@num2str,{expdirstruct.screen_type},'uni',0))

%% pull just the HHMMSS from experiment date time
timeofday =[];
datetime = {expdirstruct.exp_datetime};
for i = 1:numel(expdirstruct)
timeofday(i,:) = datetime{i}(end-5:end);
expdirstruct(i).timeofday = timeofday(i,:);
end


nowdatetime = datestr(now,'yyyymmddTHHMMSS');
save([savefile,'_addedtimeofday_',nowdatetime,'.mat'],'expdirstruct');
%% make csv file for all experiments (use for pulling data for metadata changes) with time of day 
fid = fopen([savefile,'_timeofday','.tsv'],'w');

fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','timeofday','linename','experimentor','gender','notes_tech','notes_behav','automated_pf');

for i = 1:numel(expdirstruct)  
  [~,expname] = fileparts(expdirstruct(i).file_system_path);
  datestr = expdirstruct(i).date;
  datetimestr = expdirstruct(i).exp_datetime;
  timeofdaystr = expdirstruct(i).timeofday;
  linename = expdirstruct(i).line;
  experimenter = expdirstruct(i).experimenter;
  gender = expdirstruct(i).gender;
  notestech = expdirstruct(i).notes_technical;
  notesbeh = expdirstruct(i).notes_behavioral;
  autopf = expdirstruct(i).automated_pf;


%   notes = experiments_eddison(i).notes_curation;
  fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,timeofdaystr,linename,experimenter,gender,notestech,notesbeh,autopf);
end

fclose(fid);
%% make csv file for all experiments (use for pulling data for metadata changes) with time of day 
fid = fopen([savefile,'_autofailcat','.tsv'],'w');

% fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','linename','automated_pf','automated_pf_category','manual_fail','manual_fail_category');
fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t \n','expname','date','datetime','linename','automated_pf','automated_pf_category');

for i = 1:numel(expdirstruct)
    [~,expname] = fileparts(expdirstruct(i).file_system_path);
    datestr = expdirstruct(i).date;
    datetimestr = expdirstruct(i).exp_datetime;
    linename = expdirstruct(i).line;
    autopf = expdirstruct(i).automated_pf;
    autopfcat = expdirstruct(i).automated_pf_category;
    manpf = expdirstruct(i).manual_fail;
%     manpfcat = expdirstruct(i).manual_fail_category;
    %   notes = experiments_eddison(i).notes_curation;
%     fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,linename,autopf,autopfcat,manpf,manpfcat);
        fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t \n',expname, datestr,datetimestr,linename,autopf,autopfcat);
end

fclose(fid);

%% list of experiments with manual fail
idx = strcmp({expdirstruct.manual_fail},'F');
idx2 = find(idx);
savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/moveto_nearlineVNCscreen_notinuse.txt';
fid = fopen(savefile,'w');
for i = 1:numel(idx2)
    fprintf(fid,'%s\n',expdirstruct(idx2(i)).file_system_path);
end
fclose(fid);