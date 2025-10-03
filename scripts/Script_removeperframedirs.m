explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/explist_socialbubblereprocessing.txt','%s');
rootdir = '/groups/branson/bransonlab/flybubble_social';

for i = 1:numel(explist)
    expdir = fullfile(rootdir,explist{i});
    rmdir(fullfile(expdir,'perframe'),'s')
end
