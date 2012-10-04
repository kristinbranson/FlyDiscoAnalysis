expdir_new = '/groups/branson/bransonlab/tracking_data/olympiad/HackHitData/GMR_MB031B_TrpA_Rig2Plate17BowlA_20120610T132136';
expdir_old = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/GMR_MB031B_TrpA_Rig2Plate17BowlA_20120610T132136';

perframedir = fullfile(expdir_new,'perframe');
perframefns = dir(fullfile(perframedir,'*.mat'));
perframefns = regexprep({perframefns.name},'\.mat$','');
isscoreorlabel = ~cellfun(@isempty,regexp(perframefns,'([sS]core)|([lL]abel)','once'));
perframefns(isscoreorlabel) = [];
perframefns = setdiff(perframefns,{'sex'});

for j = 1:numel(perframefns),
  fn = perframefns{j};
  dataold = load(fullfile(expdir_old,'perframe',[fn,'.mat']));
  datanew = load(fullfile(expdir_new,'perframe',[fn,'.mat']));
  
  for i = 1:numel(dataold.data),
    err = max(abs(dataold.data{i}-datanew.data{i}));
    if err > 0,
      fprintf('%s(%d): %f\n',fn,i,err);
    end
  end
end