%Deubg get all experiment directories 
%% setup path
modpath
% metadatafile = 'Metadata.xml';
% %% pull expdirs and load all metdata
rootdatadir = '/groups/rubin/data0/rubinlab/flydisco_data/schrettercBub';

%Setting up inputs...
% %inputs to getExperimentDirsFlyDisco: 'metadatafile','Metadata.xml','expdirname','*','line_name','*', ...
%     'date','*','nflies',false,'autocheckin',false,'movielength',false,'TrajNum',false,'autocheckcomplete',false);

% savefile = ''; Example for checking
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','movielength',false,'TrajNum',false,'autocheckin',false,'autocheckcomplete',true);
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','movielength',false,'TrajNum',false,'autocheckin',false,'autocheckcomplete',true,'manualcheck',true);

savefile = '/groups/rubin/home/schretterc/Documents/Test_15mmBubblesRGB/expdirs_multiBubFails';
%savefile = '/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/VisionProject_TuTuAsANDIB112/TuTuAs/Female/aIPgTrpAChrim/TuTuA_2SS2';
% [expdirstruct] = getExperimentDirsFlyDiscoKS(rootdatadir,'metadatafile','metaData.xml','line_name','_norpASSempty*','date','20230321*');
[expdirstruct] = getExperimentDirsFlyDiscoKS(rootdatadir,'metadatafile','metaData.xml','date','20240903*');

% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,varagin);

% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','expdirname','VNC','autocheckin',true,'autocheckcomplete',true,'TrajNum',true);
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'metadatafile','Metadata.xml','autocheckin',true);
% save raw pull data
nowdatetime = datestr(now,'yyyymmddTHHMMSS');
save([savefile,'_',nowdatetime,'.mat'],'expdirstruct');

%% make a list for the ones that failed autochecks complete and are a certain screen type
%only uncomment below if more than 1 screen type used
% expdirstructOG = expdirstruct;

% screen_type 1
idx = strcmp({expdirstruct.screen_type},'multibub2');
expdirstruct = expdirstruct(idx);

%change below to F if want failed experiments
pidx = find([expdirstruct.automated_pf] == 'P'); 

explist = {expdirstruct(pidx).file_system_path};
save([savefile,'_multiBub2_20241105List.mat'],'explist')

% % screen_type 2
% expdirstruct = expdirstructOG;
% idx = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
% expdirstruct = expdirstruct(idx);
% 
% pidx = find([expdirstruct.automated_pf] == 'P');
% 
% explist = {expdirstruct(pidx).file_system_path};
% save([savefile,'_VNCpassAuto.mat'],'explist')

% all list
% expdirstruct = expdirstructOG;
% idx1 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
% % idx2 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC2');
% % idxall = idx1 + idx2;
% expdirstruct = expdirstruct(logical(idxall));
% pidx = find([expdirstruct.automated_pf] == 'P');
% 
% explist = {expdirstruct(pidx).file_system_path};
% save([savefile,'_VNCallpassAuto.mat'],'explist');

savefiletext = [savefile,'_multiBub2_20241105.txt'];
fid = fopen(savefiletext,'w');
for i = 1:numel(explist)
    fprintf(fid,'%s\n',explist{i});
end
fclose(fid);

%% select on expdirs with certain led protocol
% idx = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
idx = strcmp({expdirstruct.screen_type},'mmultibub2');
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
pidx = find([expdirstruct2.automated_pf] == 'P');
expdirstruct3 = expdirstruct2(pidx);
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
savefiletext = [savefile,'_VNCall_passAandM_metadata.mat'];
save(savefiletext,'metadata');


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
fid = fopen([savefile,'.tsv'],'w');

fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','linename','automated_pf','automated_pf_category','manual_fail','manual_fail_category');

for i = 1:numel(expdirstruct)
    [~,expname] = fileparts(expdirstruct(i).file_system_path);
    datestr = expdirstruct(i).date;
    datetimestr = expdirstruct(i).exp_datetime;
    linename = expdirstruct(i).line;
    autopf = expdirstruct(i).automated_pf;
    autopfcat = expdirstruct(i).automated_pf_category;
    manpf = expdirstruct(i).manual_fail;
    manpfcat = expdirstruct(i).manual_fail_category;
    %   notes = experiments_eddison(i).notes_curation;
    fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,linename,autopf,autopfcat,manpf,manpfcat);
end

fclose(fid);