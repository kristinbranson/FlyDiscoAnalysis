modpath;
lblfile = '/groups/branson/home/bransonlab/APTlabelfiles/20220826_Grone_60K_aligned.lbl';
explistfile = '/groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/expdirs_VNC2_20220921_VNCallpass.txt';
settingsdir = 'settings';
analysis_protocol = '20220913_flybubble_LED_VNC2';
rootoutdir = '/groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/testdata';
sshhost = 'bransonlab@login1.int.janelia.org';

expdirlist = importdata(explistfile);

nexps = numel(expdirlist);

%%

success = false(nexps,1);
msgs = cell(nexps,1);
jobid = nan(nexps,1);

for expi = 21:nexps,

  [success(expi),msgs{expi},jobid(expi)] = FlyDiscoAPTTrack(expdirlist{expi},'settingsdir',settingsdir,...
    'analysis_protocol',analysis_protocol,...
    'sshhost',sshhost,...
    'dryrun',false,...
    'verbose',2,...
    'dowait',false,...
    'dooverwrite',false);

end

%% check status

jobresult = nan(nexps,1);
lsf_status = cell(nexps,1);
idx = 21:nexps;
[jobresultcurr,lsf_statuscurr] = get_bsub_job_status(jobid(idx),sshhost);
idx1 = ~isnan(jobresultcurr);
jobresult(idx(idx1)) = jobresultcurr(idx1);
lsf_status(idx(idx1)) = lsf_statuscurr(idx1);

%% plot status

expdir = '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigA_20220419T094935';
movfile = fullfile(expdir,'movie.ufmf');
trkfile = '/groups/branson/home/bransonk/tracking/code/APT/deepnet/VNC2_EXT_VGLUT-GAL4_RigA_20220419T094935.trk';

trk = TrkFile.load(trkfile);
[readframe] = get_readframe_fcn(movfile);

%%

roottestoutdir = 'aptchecktrkresults20220923';

hfigsample = figure(123);
nsamples = 10;
set(hfigsample,'Units','pixels','Position',[100,100,1144,1144]);
hfigxy = figure(124);
set(hfigxy,'Units','pixels','Position',[100,100,1800,700]);
maxvel = 10;
trkfileexists = false(nexps,1);
for expi = 1:nexps,
%   if isnan(jobid(expi)),
%     continue;
%   end
  expdir = expdirlist{expi};
  [~,expname] = fileparts(expdir);
  testoutdir = fullfile(roottestoutdir,expname);
  movfile = fullfile(expdir,'movie.ufmf');
  trkfile = fullfile(expdir,'apt.trk');
  trkfileexists(expi) = exist(trkfile,'file') > 0;
  parttrkfile = [trkfile,'.part'];

  if exist(trkfile,'file'),
    trk = TrkFile.load(trkfile);
    fprintf('%d Loading full trk file %s\n',expi,trkfile);
  elseif exist(parttrkfile,'file'),
    trk = TrkFile.load(parttrkfile);
    fprintf('%d Loading part trk file %s\n',expi,parttrkfile);
  else
    continue;
  end
  if ~exist(testoutdir,'dir'),
    mkdir(testoutdir);
  end
  [readframe,~,fid] = get_readframe_fcn(movfile);

  Ts = trk.getEndFrame();
  maxT = max(Ts);
  samplets = unique(round(linspace(1,maxT-1,nsamples)));

  for samplei = 1:numel(samplets),
    t = samplets(samplei);

    im = readframe(t);
    [~,xy] = trk.getPTrkFrame(t,'collapse',true);

    clf(hfigsample);
    hax = axes('Position',[0,0,1,1],'Parent',hfigsample);
    imagesc(im,'Parent',hax,[0,255]); 
    axis(hax,'image','off');
    hold(hax,'on');
    plot(hax,squeeze(xy(:,1,:))',squeeze(xy(:,2,:))','.');
    colormap(hax,'gray');
    npts = size(xy,1);
    colors = jet(npts);
    set(hax,'ColorOrder',colors);
    if samplei == 1,
      truesize(hfigsample);
    end
    drawnow;

    outfilenames = SaveFigLotsOfWays(hfigsample,fullfile(testoutdir,sprintf('SampleFrame%05d',t)),{'png'});

  end

  clf(hfigxy);
  hax = axes('Position',[.05,.1,.89,.87],'Parent',hfigxy);

  colors = jet(trk.ntracklets)*.7;
  for tgt = 1:trk.ntracklets,
    [xy,~,fr] = trk.getPTrkTgt(tgt);
    meanvel = squeeze(mean(sqrt(sum((xy(:,:,2:end)-xy(:,:,1:end-1)).^2,2)),1,'omitnan'));
    plot(hax,fr(1:end-1)',tgt-1+meanvel/maxvel,'-','Color',colors(tgt,:));
    hold(hax,'on');
  end
  set(hax,'YLim',[0,trk.ntracklets],'XLim',[0,maxT]);
  xlabel(hax,'Frame');
  ylabel(hax,sprintf('Mean landmark speed (px/fr) / %d, target',maxvel));
  box(hax,'off');
  drawnow;
  outfilenames = SaveFigLotsOfWays(hfigxy,fullfile(testoutdir,'LandmarkSpeed'),{'pdf','png'});

  fclose(fid);

end