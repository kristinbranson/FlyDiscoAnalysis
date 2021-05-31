addpath JAABA/misc;
addpath JAABA/filehandling;
%expdir = '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS46706_RigA_20210512T132556';
expdir = 'testdata/VNC_JRC_SS46706_RigA_20210512T132556';
analysis_protocol = '20210531_flybubble_LED';
settingsdir = 'settings';
FlyDiscoComputePerFrameStats(expdir,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);


%%

moviefile = fullfile(expdir,'movie.ufmf');
[readframe,~,fid,headerinfo] = get_readframe_fcn(moviefile);
td = load(fullfile(expdir,'wingtracking_results.mat'));

%%
f = 9649;

figure(1234);
clf;
imagesc(readframe(f));
axis image;
colormap gray;
hold on;

hflies = gobjects(numel(td.trx));
hwings = gobjects(numel(td.trx));
for fly = 1:numel(td.trx),
  hflies(fly) = plot(nan,nan,'-');
  hwings(fly) = plot(nan,nan,'-','Color',get(hflies(fly),'Color'));
  pos = struct;
  pos.x = td.trx(fly).x(f + td.trx(fly).off);
  pos.y = td.trx(fly).y(f + td.trx(fly).off);
  pos.theta = td.trx(fly).theta(f + td.trx(fly).off);
  pos.a = td.trx(fly).a(f + td.trx(fly).off);
  pos.b = td.trx(fly).b(f + td.trx(fly).off);
  pos.xwingl = td.trx(fly).xwingl(f + td.trx(fly).off);
  pos.ywingl = td.trx(fly).ywingl(f + td.trx(fly).off);
  pos.xwingr = td.trx(fly).xwingr(f + td.trx(fly).off);
  pos.ywingr = td.trx(fly).ywingr(f + td.trx(fly).off);
  updatewingedfly(hflies(fly),hwings(fly),pos);
  plot(pos.xwingl,pos.ywingl,'.','Color',get(hflies(fly),'Color'));
end


%%

fclose(fid);