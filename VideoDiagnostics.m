function video_diagnostics = VideoDiagnostics(expdir,varargin)

% compute average frame rate from first and last timestamps and number of frames

video_diagnostics = struct;

%% parse parameters
[analysis_protocol,settingsdir,datalocparamsfilestr,nboxes_lims,nax_ufmf,hfig,radius_ufmf,fig_pos,DEBUG] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'nboxes_lims',[1000,2000,3000,4000,5000,6000,7000],...
  'nax_ufmf',[3,5],...
  'hfig',1,...
  'radius_ufmf',100,...
  'fig_pos',[1,1,1548,926],...
  'DEBUG',false);

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% load timestamps from ctrax_results

moviefile = fullfile(expdir,dataloc_params.moviefilestr);
isvideo = exist(moviefile,'file');
if ~isvideo,
  video_diagnostics.mean_frame_rate = -1;
  for i = 1:numel(nboxes_lims),
    fn = sprintf('fracFramesWithNBoxes%05d',nboxes_lims(i));
    video_diagnostics.(fn) = -1;
  end
else
  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
  
  % get average frame rate
  [~,t0] = readframe(1);
  [~,t1] = readframe(nframes);
  video_diagnostics.mean_frame_rate = (nframes-1)/(t1-t0);
  
  % get nboxes per frame
  nboxes = ufmf_read_nboxes(headerinfo,1:nframes);
  
  % how many frames have more than nboxes_lims boxes
  for i = 1:numel(nboxes_lims),
    fn = sprintf('fracFramesWithNBoxes%05d',nboxes_lims(i));
    video_diagnostics.(fn) = nnz(nboxes >= nboxes_lims(i)) / nframes;
  end

  % make an image of the frames with the most boxes
  figure(hfig);
  clf;
  set(hfig,'Position',fig_pos);
  hax = createsubplots(nax_ufmf(1),nax_ufmf(2),.01);
  hax = reshape(hax,nax_ufmf([2,1]))';
  for i = 1:prod(nax_ufmf),
    [nboxes_curr,f] = max(nboxes);
    if nboxes_curr < 0,
      delete(hax(i:end));
      break;
    end
    [im,~,~,~,mu] = readframe(f);
    im_rgb = uint8(repmat(im,[1,1,3]));
    tmp = im;
    tmp(im~=mu) = 0;
    im_rgb(:,:,2) = tmp;
    image(im_rgb,'Parent',hax(i));
    hold(hax(i),'on');
    text(1,1,sprintf('%05d',f),'Color','g','Parent',hax(i),'HorizontalAlignment','left','VerticalAlignment','top','FontSize',18);
    axis(hax(i),'image','off');
    nboxes(max(1,f-radius_ufmf):min(nframes,f+radius_ufmf)) = -1;
  end
  
  fclose(fid);
  
end

%% save image

if isvideo && ~DEBUG,
  savename = fullfile(expdir,dataloc_params.videodiagnosticsimagefilestr);
  try
    if exist(savename,'file'),
      delete(savename);
    end
    save2png(savename,hfig);
  catch ME,
    warning('Could not write video diagnostics image to file %s:\n%s',savename,getReport(ME));
  end

end

%% save to mat file
if ~DEBUG,
  videodiagnosticsmatfilename = fullfile(expdir,dataloc_params.videodiagnosticsmatfilestr);
  try
    save(videodiagnosticsmatfilename,'-struct','video_diagnostics');
  catch ME,
    warning('Could not save video diagnostics to mat file %s:\n%s',videodiagnosticsmatfilename,getReport(ME));
  end
end

%% write to text file

videodiagnosticsfilename = fullfile(expdir,dataloc_params.videodiagnosticsfilestr);
if DEBUG,
  fid = 1;
else
  fid = fopen(videodiagnosticsfilename,'w');
end
if fid < 0,
  warning('Could not open video diagnostics file %s for writing',videodiagnosticsfilename);
else
  fns = fieldnames(video_diagnostics);
  for i = 1:numel(fns),
    fprintf(fid,'%s',fns{i});
    fprintf(fid,',%f',video_diagnostics.(fns{i}));
    fprintf(fid,'\n');
  end
  if ~DEBUG,
    fclose(fid);
  end
end