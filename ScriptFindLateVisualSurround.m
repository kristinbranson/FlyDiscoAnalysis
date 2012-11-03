%% load in the data

load /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/gal4screen/FractimeDataAllGAL4s20121030.mat data;

%% parameters

manual_curator = 'bransonk';
timestamp = datestr(now,'yyyymmddTHHMMSS');

%% read in all the video diagnostics

tmp = load(fullfile(data(1).file_system_path,'video_diagnostics.mat'));
fns = fieldnames(tmp);

for i = 1:numel(data),
  filename = fullfile(data(i).file_system_path,'video_diagnostics.mat');
  if ~exist(filename,'file'),
    fprintf('File %s does not exist\n',filename);
    for j = 1:numel(fns),
      fnin = fns{j};
      fnout = ['video_diagnostics_',fns{j}];
      data(i).(fnout) = nan;
    end
    continue;
  end
  tmp = load(fullfile(data(i).file_system_path,'video_diagnostics.mat'));
  for j = 1:numel(fns),
    fnin = fns{j};
    fnout = ['video_diagnostics_',fns{j}];
    data(i).(fnout) = tmp.(fnin);
  end
  if mod(i,100) == 0,
    fprintf('%d\n',i);
  end
end

idxmissing = find(isnan([data.video_diagnostics_fracFramesWithNBoxes07000]));
fprintf('The following experiments have no video diagnostics:\n');
fprintf('%s\n',data(idxmissing).file_system_path);

%% initialize

isproblem = nan(1,numel(v));
notes = cell(1,numel(v));

%% look at experiments sorted by fracFramesWithNBoxes07000

[v,order] = sort([data.video_diagnostics_fracFramesWithNBoxes07000],2,'descend');
firsti = find(~isnan(v),1);

[isproblem,notes] = ExamineVideoDiagnostics(data,order(firsti:end),'isproblem',isproblem,'notes',notes);

%% look at experiments sorted by fracFramesWithNBoxes01000

[v,order] = sort([data.video_diagnostics_fracFramesWithNBoxes01000],2,'descend');
firsti = find(~isnan(v),1);

[isproblem,notes] = ExamineVideoDiagnostics(data,order(firsti:end),'isproblem',isproblem,'notes',notes);


%% look at experiments sorted by ufmf_diagnostics_summary_maxBandWidth

[v,order] = sort([data.ufmf_diagnostics_summary_maxBandWidth],2,'descend');
firsti = find(~isnan(v),1);

[isproblem,notes] = ExamineVideoDiagnostics(data,order(firsti:end),'isproblem',isproblem,'notes',notes);

%% look at experiments sorted by ctrax_diagnostics_nlarge_ignored

[v,order] = sort([data.ctrax_diagnostics_nlarge_ignored],2,'descend');
firsti = find(~isnan(v),1);

[isproblem,notes] = ExamineVideoDiagnostics(data,order(firsti:end),'isproblem',isproblem,'notes',notes);

%% look at all experiments with technical notes about the surround

isnotes = ~cellfun(@isempty,regexpi({data.notes_technical},'visual|surround')) & isproblem~=1;
idxnotes = find(isnotes);

[isproblem,notes] = ExamineVideoDiagnostics(data,idxnotes,'isproblem',isproblem,'notes',notes);

%% look at all other experiments from sets with problems

[sets,~,setidx] = unique({data.set});

issameset = false(1,numel(data));
for i = find(isproblem==1),
  issameset = issameset | (setidx(i) == setidx & isproblem~=1);
end
idxsameset = find(issameset);

[isproblem,notes] = ExamineVideoDiagnostics(data,idxsameset,'isproblem',isproblem,'notes',notes);

%% add in visual surround problem for all empty notes for problems

idx = isproblem==1 & cellfun(@(x) isempty(x) || iscell(x)&&isempty(x{1}),notes);
[notes{idx}] = deal('Visual surround problem');

%% output results to a file

savefilename = sprintf('VideoDiagnosticsCuration%s.tsv',timestamp);

fid = fopen(savefilename,'w');
if fid < 1,
  error('Could not open file %s for writing. Make sure it is not open in another program.',savefilename);
end
fprintf(fid,'experiment\tmanual_pf\tmanual_curator\tmanual_curation_date\tnotes_curation\n');
for i = find(isproblem==1),
  if isempty(data(i).notes_curation),
    notes_curation_curr = {};
  elseif ~iscell(data(i).notes_curation),
    notes_curation_curr = {data(i).notes_curation};
  else
    note_curation_curr = data(i).notes_curation;
  end
  removeidx = ismember(notes_curation_curr,{'NULL','None'});
  notes_curation_curr(removeidx) = []; %#ok<SAGROW>
  
  notes_curation_curr{end+1} = notes{i}; %#ok<SAGROW>
  if iscell(notes_curation_curr),
    notes_curation_curr = sprintf('%s\\n',notes_curation_curr{:});
    % remove extra \n
    notes_curation_curr = notes_curation_curr(1:end-2);
  else
    notes_curation_curr = data(i).notes_curation;
  end
  % only \n?
  if numel(notes_curation_curr) >= 2 && strcmp(notes_curation_curr(end-1:end),'\n'),
    notes_curation_curr = notes_curation_curr(1:end-2);
  end
  % \n"?
  if numel(notes_curation_curr) >= 3 && strcmp(notes_curation_curr(end-2:end),'\n"'),
    notes_curation_curr = [notes_curation_curr(1:end-3),'"'];
  end
  fprintf(fid,'%s\t%s\t',data(i).experiment_name,'F');
  % print curation time, curator if experiment variables
  fprintf(fid,'%s\t%s\t',manual_curator,timestamp);
  % print notes
  fprintf(fid,'%s\n',notes_curation_curr);
end
fclose(fid);

save(sprintf('VideoDiagnosticsCuration%s.mat',timestamp),'notes','isproblem');