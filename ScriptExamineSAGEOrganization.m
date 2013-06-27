% score field names
tmp = SAGE.Lab('olympiad');
tmp = tmp.assay('bowl');
bowlAssay_data = tmp.dataSet('score');
scorefields = bowlAssay_data.fields;
fprintf('Score field names:\n');
fprintf('%s\n',scorefields.name);

% data field names
bowlAssay_data = tmp.dataSet('data');
datafields = bowlAssay_data.fields;
fprintf('\nData field names:\n');
fprintf('%s\n',datafields.name);

if ~isempty(setdiff({datafields.name},{scorefields.name})) || ...
    ~isempty(setdiff({scorefields.name},{datafields.name})),
  fprintf('Score and data field names do not match\n');
end

% score data names
tmpscore = SAGEGetBowlData('dataset','score','experiment_name','FlyBowl_GMR_86E11_AE_01_TrpA_Rig1Plate10BowlC_20110323T132841');
scoredata = setdiff(fieldnames(tmpscore),{scorefields.name});
fprintf('\nScore data names:\n');
fprintf('%s\n',scoredata{:})

% data only data names
tmpdata = SAGEGetBowlData('dataset','data','experiment_name','FlyBowl_GMR_86E11_AE_01_TrpA_Rig1Plate10BowlC_20110323T132841');
fprintf('\nData only data names:\n');
datadata = setdiff(fieldnames(tmpdata),union({scorefields.name},scoredata));
fprintf('%s\n',datadata{:});

% compare all bias diagnostics with SAGE data
bd = ReadParams('/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/GMR_86E11_AE_01_TrpA_Rig1Plate10BowlC_20110323T132841/bias_diagnostics.txt');
biasdiagnsoticscores = regexp(scoredata,'^bias_diagnostics_(.*)$','once','tokens');
biasdiagnsoticscores(cellfun(@isempty,biasdiagnsoticscores)) = [];
biasdiagnsoticscores = cellfun(@(x) x{1},biasdiagnsoticscores,'UniformOutput',false);

biasdiagnsoticdata = regexp(datadata,'^bias_diagnostics_(.*)$','once','tokens');
biasdiagnsoticdata(cellfun(@isempty,biasdiagnsoticdata)) = [];
biasdiagnsoticdata = cellfun(@(x) x{1},biasdiagnsoticdata,'UniformOutput',false);

if ~isempty(setdiff(union(biasdiagnsoticscores,biasdiagnsoticdata),fieldnames(bd))) && ...
    ~isempty(setdiff(fieldnames(bd),union(biasdiagnsoticscores,biasdiagnsoticdata))),
  fprintf('Difference in bias diangostic fields in SAGE vs file\n');
end