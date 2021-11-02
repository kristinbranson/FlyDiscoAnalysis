% [succeeded,aviname] = make_apt_result_movie('param',value,...)
% Options:
% 'moviename': name of raw movie to annotate
% 'aptname': name of mat file containing trajectories, or struct of
% trajectories
% 'aviname': name of movie to output to
% 'colors': colors to plot each fly
% 'nzoomr': number of rows of zoom fly boxes
% 'nzoomc': number of columns of zoom fly boxes
% 'boxradius': radius of the zoom fly box in pixels
% 'taillength': number of frames of trajectory to plot behind each fly
% 'zoomflies': flies to zoom in on
% 'fps': frames per second of output movie
% 'maxnframes': number of frames to output
% 'firstframe': first frame to output
% 'compression': compressor to use when outputting (pc only, will be set to
% 'none' for linux). 
% 'figpos': position of figure
% if any parameters are not given, the user will be prompted for these
function [succeeded,aviname,figpos,height,width] = make_apt_result_movie(varargin)

succeeded = false;
defaults.boxradius = 1.5;
defaults.taillength = 100;
defaults.fps = 20;
defaults.zoomflies = [];
defaults.nzoomr = 5;
defaults.nzoomc = 3;
defaults.compression = 'None';
defaults.headlandmark = 1;
defaults.taillandmark = 7;
defaults.apttailalpha = .5;
defaults.colormult = .7;
if isdisplay,
  defaults.lmkmarkersize = 6;
else
  defaults.lmkmarkersize = 10;
end
allowedcompressions = {'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE', 'None','Uncompressed AVI','Motion JPEG AVI'};
useVideoWriter = exist('VideoWriter','file');
[moviename,aptname,trxname,aviname,colors,zoomflies,nzoomr,nzoomc,boxradius,...
  taillength,fps,maxnframes,firstframes,compression,figpos,movietitle,...
  useVideoWriter,headlandmark,taillandmark,...
  avifileTempDataFile,titletext,showtimestamps,dynamicflyselection,...
  doshowsex,doplotapttail,doplotcentertail,doplotwings,doflipud,dofliplr] = ...
  myparse(varargin,'moviename','','aptname','','trxname','','aviname','','colors',[],'zoomflies',[],'nzoomr',nan,'nzoomc',nan,...
  'boxradius',nan,'taillength',nan,'fps',nan,'maxnframes',nan,'firstframes',[],'compression','',...
  'figpos',[],'movietitle','','useVideoWriter',useVideoWriter,...
  'headlandmark',[],'taillandmark',[],...
  'avifileTempDataFile','',...
  'titletext',true,...
  'showtimestamps',false, ...
  'dynamicflyselection',true,...
  'doshowsex',true,...
  'doplotapttail',true,...
  'doplotcentertail',false,...
  'doplotwings',true,...
  'flipud',false,'fliplr',false);

if ~ischar(compression),
  compression = '';
end
if ~isempty(compression) && ~any(strcmpi(compression,allowedcompressions)),
  fprintf('Unknown compressor %s\n',compression);
  compression = '';
end

if ~ischar(moviename) || isempty(moviename) || ~exist(moviename,'file'),
  fprintf('Choose raw movie to annotate\n');
  helpmsg = 'Choose raw movie to annotate';
  [movienameonly,moviepath] = uigetfilehelp({'*.fmf';'*.sbfmf';'*.avi'},'Choose raw movie to annotate','','helpmsg',helpmsg);
  if ~ischar(movienameonly),
    return;
  end
  moviename = [moviepath,movienameonly];
else
  [moviepath,movienameonly] = split_path_and_filename(moviename);
end
[readframe1,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
  
if doflipud && dofliplr,
  readframe = @(x) fliplr(flipud(readframe1(x)));  %#ok<FLUDLR>
elseif doflipud,
  readframe = @(x) flipud(readframe1(x));
elseif dofliplr
  readframe = @(x) fliplr(readframe1(x));
else
  readframe = readframe1;
end
  
if fid < 0,
  uiwait(msgbox(sprintf('Could not read in movie %s',moviename)));
  return;
end

if ~ischar(aptname) && isequal(class(aptname),'TrkFile'),
  apttrk = aptname;
else
  if ~ischar(aptname) || isempty(aptname) || ~exist(aptname,'file'),
    fprintf('Choose mat file containing flies'' APT trajectories corresponding to movie %s\n',moviename);
    [~,ext] = splitext(moviename);
    aptname = [moviename(1:end-length(ext)+1),'trk'];
    helpmsg = sprintf('Choose APT trk file to annotate the movie %s with',moviename);
    [trknameonly,trkpath] = uigetfilehelp('*.trk',sprintf('Choose APT trk file for %s',movienameonly),aptname,'helpmsg',helpmsg);
    if ~ischar(trknameonly),
      return;
    end
    aptname = [trkpath,trknameonly];
  end
  if ~exist(aptname,'file'),
    loadsucceeded = false;
    return;
  end
  try
    apttrk = TrkFile.load(aptname);
  catch 
    loadsucceeded = false;
    return;
  end
end
ntgts = apttrk.ntlts;

if ~ischar(trxname) || isempty(trxname) || ~exist(trxname,'file'),
  fprintf('Choose mat file containing flies'' Ctrax/FlyTracker trajectories corresponding to movie %s\n',moviename);
  [~,ext] = splitext(moviename);
  trxname = [moviename(1:end-length(ext)+1),'mat'];
  helpmsg = sprintf('Choose Ctrax/FlyTracker trx file to annotate the movie %s with. Hit cancel to not use a trx file',moviename);
  [trxnameonly,trxpath] = uigetfilehelp('*.mat',sprintf('Choose Ctrax/FlyTracker trx file for %s',movienameonly),trxname,'helpmsg',helpmsg);
  if ischar(trxnameonly),
    trxname = [trxpath,trxnameonly];
  else
    trxname = '';
  end
end
istimestamps = false;
if ~isempty(trxname),
  istrx = true;
  istimestamps = true;
  [trx,trxname,loadsucceeded,timestamps] = load_tracks(trxname);
  if ~loadsucceeded,
    return;
  end
else
  istrx = false;
  if isfield(headerinfo,'timestamps'),
    timestamps = headerinfo.timestamps;
    istimestamps = true;
  end
end

isbodylandmarks = ~isempty(headlandmark) && ~isempty(taillandmark);
if isbodylandmarks,
  bodylandmarks = [headlandmark,taillandmark];
else
  bodylandmarks = 1:apttrk.npts;
end

% output avi
haveaviname = false;
if ischar(aviname) && ~isempty(aviname)
  [~,ext] = splitext(aviname);
  if strcmpi(ext,'.avi'),
    haveaviname = true;
  end
end
if ~haveaviname,
  fprintf('Choose avi file to output annotated version of %s\n',moviename);
  [base,~] = splitext(movienameonly);
  aviname = [moviepath,'ctraxresults_',base,'.avi'];
  helpmsg = {};
  helpmsg{1} = 'Choose avi file to write annotated movie to.';
  helpmsg{2} = sprintf('Raw movie input: %s',moviename);
  helpmsg{3} = sprintf('APT trk file name input: %s',aptname);
  [avinameonly,avipath] = uiputfilehelp('*.avi',sprintf('Choose output avi for %s',movienameonly),aviname,'helpmsg',helpmsg);
  if ~ischar(avinameonly),
    return;
  end
  aviname = [avipath,avinameonly];
end

% make sure sex is set
doshowsex = doshowsex && istrx && isfield(trx,'sex') && ~any(cellfun(@isempty,{trx.sex}));
if doshowsex,
  sexes = {};
  for i = 1:numel(trx),
    if ~iscell(trx(i).sex),
      trx(i).sex = num2cell(trx(i).sex);
    end
    sexes = union(sexes,unique(trx(i).sex));
  end
  sexes = upper(sexes);
  sexes = sort(sexes);
  sexmarkers = {'none','*','x','o','+','d','s','p','h'};
  if numel(sexes) > numel(sexmarkers),
    sexmarkers = repmat(sexmarkers,[1,ceil(numel(sexes)/numel(sexmarkers))]);
  end
  sexmarkers = sexmarkers(1:numel(sexes));
end

doplotwings = doplotwings && istrx && all(isfield(trx,{'xwingl','ywingl','xwingr','ywingr'}));
doplotwings_perfly = repmat(doplotwings,[1,ntgts]);
if doplotwings && isfield(trx,'wingtype'),
  for i = 1:numel(trx),
    doplotwings_perfly(i) = all(strcmpi(trx(i).wingtype,'full'));
  end
end

nzoom = length(zoomflies);
if nzoom > 0,
  if max(zoomflies) > ntgts || min(zoomflies) < 1,
    uiwait(msgbox('Illegal values for zoomflies'));
    return;
  end
end

prompts = {};
defaultanswers = {};
if isnan(nzoomr),
  prompts{end+1} = 'Number of rows of zoomed fly boxes';
  defaultanswers{end+1} = num2str(defaults.nzoomr);
end
if isnan(nzoomc) && nzoom == 0,
  prompts{end+1} = 'Number of columns of zoomed fly boxes';
  defaultanswers{end+1} = num2str(defaults.nzoomc);
end
if isnan(boxradius),
  prompts{end+1} = 'Radius of zoomed fly box (in pixels)';
  
  
  maxd = [];
  for fly = 1:ntgts,
    ff = apttrk.startframes(fly);
    ef = apttrk.endframes(fly);
    nf = ef-ff+1;
    xy = apttrk.getPTrkTgt(fly);
    xy = xy(:,:,ff:ef);
    if istrx,
      centerx = trx(fly).x;
      centery = trx(fly).y;
    else
      centerx = mean(xy(bodylandmarks,1,ff:ef),1);
      centery = mean(xy(bodylandmarks,2,ff:ef),1);
    end
    
    d = sqrt((xy(:,1,:)-reshape(centerx,[1,1,nf])).^2 + (xy(:,2,:)-reshape(centery,[1,1,nf])).^2);
    maxd = [maxd;squeeze(max(d,[],1))]; %#ok<AGROW>
    
  end
  defaultanswers{end+1} = num2str(median(maxd)*defaults.boxradius);
end
if isnan(taillength),
  prompts{end+1} = 'Length of plotted tail trajectory (in frames)';
  defaultanswers{end+1} = num2str(defaults.taillength);
end
if isnan(fps) && itrx && isfield(trx,'fps'),
  fps = trx(1).fps;
end
if isnan(fps),
  prompts{end+1} = 'Output movie frames per second';
  defaultanswers{end+1} = num2str(defaults.fps);
end
if isnan(maxnframes),
  prompts{end+1} = 'Max number of frames to output';
  defaultanswers{end+1} = num2str(nframes);
end
if isempty(firstframes),
  prompts{end+1} = 'First frame to output';
  defaultanswers{end+1} = num2str(1);
end

compressionprompt = ['Compressor (must be one of ',...
  sprintf('%s, ',allowedcompressions{1:end-1}),allowedcompressions{end},')'];
if ~ispc && ~useVideoWriter,
  compression = 'None';
elseif isempty(compression),
  prompts{end+1} = compressionprompt;
  defaultanswers{end+1} = defaults.compression;
end
if ~isbodylandmarks && ~istrx,
  headlandmarkprompt = sprintf('Head landmark index (1 to %d)',apttrk.npts);
  taillandmarkprompt = sprintf('Tail landmark index (1 to %d)',apttrk.npts);
  prompts{end+1} = headlandmarkprompt;
  defaultanswers{end+1} = defaults.headlandmark;  
  prompts{end+1} = taillandmarkprompt;
  defaultanswers{end+1} = defaults.taillandmark;  
end


if ~isempty(prompts),
  while true,
    answers = inputdlg(prompts,'make_apt_result_movie parameters',1,defaultanswers);
    if isempty(answers),
      return;
    end
    failed = false;
    for i = 1:length(answers),
      if strcmp(prompts{i},compressionprompt)
        j = strmatch(answers{i},allowedcompressions);  %#ok<MATCH2>
        if isempty(j),
          uiwait(msgbox(sprintf('Illegal compressor: %s',answers{i})));
          failed = true;
          break;
        end
        compression = allowedcompressions{j};
        answers{i} = allowedcompressions{j};
        defaultanswers{i} = compression; %#ok<AGROW>
        continue;
      end
      answers{i} = str2double(answers{i});
      if isempty(answers{i}) || answers{i} < 0,
        uiwait(msgbox('All answers must be positive numbers'));
        failed = true;
        break;
      end
      switch prompts{i},
        case 'Number of rows of zoomed fly boxes',
          nzoomr = ceil(answers{i});
          
          % check that nzoomr is not too big
          if nzoom > 0,
            if ~isnan(nzoomc) && nzoomr > ceil(nzoom / nzoomc),
              fprintf('nzoomr = %d > ceil(nzoom = %d / nzoomc = %d), decreasing nzoomr\n',nzoomr,nzoom,nzoomc);
              nzoomr = ceil(nzoom/nzoomc);
            elseif nzoomr > nzoom,
              fprintf('nzoomr = %d > nflies to plot = %d, decreasing nzoomr\n',nzoomr,nzoom);
              nzoomr = nzoom;
            end
          else
            if ~isnan(nzoomc) && nzoomr > ceil(ntgts / nzoomc),
              fprintf('nzoomr = %d > ceil(nflies = %d / nzoomc = %d), decreasing nzoomr\n',nzoomr,ntgts,nzoomc);
              nzoomr = ceil(ntgts/nzoomc);
            elseif nzoomr > ntgts,
              fprintf('nzoomr = %d > nflies = %d, decreasing nzoomr\n',nzoomr,ntgts);
              nzoomr = ntgts;
            end
          end
          answers{i} = nzoomr;
        case 'Number of columns of zoomed fly boxes',
          nzoomc = ceil(answers{i});
          % check that nzoomr is not too big
          if nzoom > 0,
            if ~isnan(nzoomr) && nzoomc > ceil(nzoom / nzoomr),
              fprintf('nzoomc = %d > ceil(nzoom = %d / nzoomr = %d), decreasing nzoomc\n',nzoomc,nzoom,nzoomr);
              nzoomc = ceil(nzoom/nzoomr);
            elseif nzoomc > nzoom,
              fprintf('nzoomc = %d > nflies to plot = %d, decreasing nzoomc\n',nzoomc,nzoom);
              nzoomc = nzoom;
            end
          else
            if ~isnan(nzoomr) && nzoomc > ceil(ntgts / nzoomr),
              fprintf('nzoomc = %d > ceil(nflies = %d / nzoomr = %d), decreasing nzoomc\n',nzoomc,ntgts,nzoomr);
              nzoomc = ceil(ntgts/nzoomr);
            elseif nzoomc > ntgts,
              fprintf('nzoomc = %d > nflies = %d, decreasing nzoomc\n',nzoomc,ntgts);
              nzoomc = ntgts;
            end
          end
          answers{i} = nzoomc;
        case 'Radius of zoomed fly box (in pixels)',
          boxradius = ceil(answers{i});
          answers{i} = boxradius;
        case 'Length of plotted tail trajectory (in frames)',
          taillength = ceil(answers{i});
          answers{i} = taillength;
        case 'Output movie frames per second',
          fps = answers{i};
          answers{i} = fps;
        case 'Max number of frames to output',
          maxnframes = min(ceil(answers{i}),nframes);
          answers{i} = maxnframes;
        case 'First frame to output',
          firstframes = max(1,ceil(answers{i}));
          answers{i} = firstframes;
        case headlandmarkprompt,
          headlandmark = min(apttrk.npts,max(1,round(answers{i})));
          if answers{i} ~= headlandmark,
            fprintf('Illegal value %f for head landmark index, rounding & bounding to %d\n',answers{i},headlandmark);
          end
        case taillandmarkprompt,
          taillandmark = min(apttrk.npts,max(1,round(answers{i})));
          if taillandmark == headlandmark,
            if headlandmark < apttrk.npts,
              taillandmark = headlandmark + 1;
            else
              taillandmark = headlandmark - 1;
            end
          end
          if answers{i} ~= headlandmark,
            fprintf('Illegal value %f for tail landmark index, rounding & bounding to %d\n',answers{i},headlandmark);
          end
          
      end
      defaultanswers{i} = num2str(answers{i}); %#ok<AGROW>
    end
    if failed,
      continue;
    end
    break;
  end
end

if ~isbodylandmarks && ~istrx,
  bodylandmarks = [headlandmark,taillandmark];
end

if nzoom > 0,
  nzoomc = ceil(nzoom/nzoomr);
end
endframes = min(nframes,firstframes+maxnframes-1);
im = readframe(firstframes(1));
[nr,nc,ncolors] = size(im);

% choose some random flies to zoom in on
nzoom = nzoomr*nzoomc;

nframesoverlap = zeros(1,ntgts);
for i = 1:numel(firstframes),
  nframesoverlap = nframesoverlap + ...
      max(0,min(endframes(i),apttrk.endframes)-max(firstframes(i),apttrk.startframes) + 1);
end
fliesmaybeplot = find(nframesoverlap > 0);

if isempty(zoomflies),
  if length(fliesmaybeplot) < nzoom,
    zoomflies = [fliesmaybeplot,nan(1,nzoom-length(fliesmaybeplot))];
    fprintf('Not enough flies to plot\n');
  else
    
    if dynamicflyselection,
      
      % choose flies to plot in each frame
      fliesplotperframe = nan(nzoom,nframes);
      
      % flies currently chosen
      fliesplotcurr = nan(nzoom,1);
      allfliesplotted = [];
      
      % loop through all frames
      
      for i = 1:numel(firstframes),
        for f = firstframes(i):endframes(i),
          
          % flies that are currently alive
          isalive = false(1,ntgts);
          for fly = 1:ntgts,
            isalive(fly) = apttrk.startframes(fly) <= f && apttrk.endframes(fly) >= f;
          end
          fliesalive = find(isalive);
          
          % remove newly dead flies
          openzoomboxes = isnan(fliesplotcurr) | ~ismember(fliesplotcurr,fliesalive);
          fliesplotcurr(openzoomboxes) = nan;
          
          % how many flies do we need to choose
          nflieschoose = nnz(openzoomboxes);
          
          % flies we can choose to add
          fliesleft = setdiff(fliesalive,fliesplotcurr);
          
          nfliesshort = max(0,nflieschoose-numel(fliesleft));
          if nfliesshort > 0,
            newfliesplot = fliesleft;
          else
            [~,order] = sort(-apttrk.endframes(fliesleft));
            newfliesplot = fliesleft(order(1:nflieschoose));
          end
          allfliesplotted = [allfliesplotted,newfliesplot]; %#ok<AGROW>
          fliesplotcurr(openzoomboxes) = [newfliesplot,nan(1,nfliesshort)];
          fliesplotperframe(:,f) = fliesplotcurr;
          
        end
      end

      zoomflies = fliesplotperframe(:,firstframes(1));
      
    else
      fliesmaybeplot = fliesmaybeplot(randperm(length(fliesmaybeplot)));
      [~,flieswithmostframes] = sort(-nframesoverlap(fliesmaybeplot));
      zoomflies = sort(fliesmaybeplot(flieswithmostframes(1:nzoom)));
    end
  end
elseif nzoom > length(zoomflies),
  zoomflies = [zoomflies,nan(1,nzoom-length(zoomflies))];
end
zoomflies = reshape(zoomflies,[nzoomr,nzoomc]);
rowszoom = floor(nr/nzoomr);


% colors of the flies
if isempty(colors),
  if dynamicflyselection && exist('allfliesplotted','var'),
    zoomfliesreal = allfliesplotted;
  else
    zoomfliesreal = zoomflies(~isnan(zoomflies));
  end
  ncolorstmp = max(ntgts,64);
  colors0 = jet(ncolorstmp).*defaults.colormult;
  colors0 = colors0(round(linspace(1,ncolorstmp,ntgts)),:);
  fliesnotplot = setdiff(1:ntgts,fliesmaybeplot);
  fliesnotzoom = setdiff(fliesmaybeplot,zoomfliesreal(:)');
  colors = nan(ntgts,3);
  coloridx = round(linspace(1,size(colors0,1),numel(zoomfliesreal)));
  colors(zoomfliesreal(:),:) = colors0(coloridx,:);
  colors0(coloridx,:) = [];
  if ~isempty(fliesnotzoom),
    coloridx = round(linspace(1,size(colors0,1),numel(fliesnotzoom)));
    colors(fliesnotzoom,:) = colors0(coloridx,:);
    colors0(coloridx,:) = [];
  end
  if ~isempty(fliesnotplot),
    colors(fliesnotplot,:) = colors0;
  end
elseif size(colors,1) ~= ntgts,
  colors = colors(modrange(0:ntgts-1,size(colors,1))+1,:);
end

% if ishandle(1),
%   close(1);
% end
fig_name = 'the make_apt_result_movie() figure' ;
fig = findobj(groot, 'Type', 'figure', 'Name', fig_name) ;
if isempty(fig) ,
    fig = figure('Name', fig_name) ;
end
clf(fig);
hax = axes(fig) ;
hold(hax,'on');
%hax = gca;
set(hax,'position',[0,0,1,1]);
axis(hax,'off');

isdisplay1 = true;
%isdisplay1 = ispc || ~strcmpi(get(fig,'XDisplay'),'nodisplay');

% corners of zoom boxes in plotted image coords
x0 = nc+(0:nzoomc-1)*rowszoom+1;
y0 = (0:nzoomr-1)*rowszoom+1;
x1 = x0 + rowszoom - 1;
y1 = y0 + rowszoom - 1;

% relative frame offset
nframesoff = apttrk.startframes - 1;

% pre-allocate
himzoom = zeros(nzoomr,nzoomc);
htailapt = zeros(1,ntgts);
htailcenter = zeros(1,ntgts);
hlandmarks = zeros(1,ntgts);
hwing = zeros(1,ntgts);
hsexmarker = zeros(1,ntgts);
scalefactor = rowszoom / (2*boxradius+1);
hzoom = zeros(nzoomr,nzoomc);
hzoomwing = zeros(nzoomr,nzoomc);
htextzoom = zeros(nzoomr,nzoomc);

if ~istrx,
  trx = apt2trx(apttrk,bodylandmarks);
end

frame_count_per_fprintf = 100 ;
for segi = 1:numel(firstframes),
  firstframe = firstframes(segi);
  endframe = endframes(segi);
  fprintf('Adding frames from segment %d, frames %d-%d\n', segi, firstframe, endframe) ;

  for frame = firstframe:endframe,
    if frame==firstframe ,
        tic_id = tic() ;
    else
      if mod(frame - firstframe,frame_count_per_fprintf) == 0,
          elapsed_time = toc(tic_id) ;  % seconds
          frame_pace = elapsed_time/frame_count_per_fprintf ;
          fprintf('Just wrote frame %d, write rate = %f s/fr\n',frame,frame_pace);
          print_matlab_memory_usage() ;
          tic_id = tic() ;
      end
    end
    
    % relative frame
    idx = frame - nframesoff;
    
    isalive = frame >= apttrk.startframes & frame <= apttrk.endframes;
    
    % draw the unzoomed image
    im = uint8(readframe(frame));
    if ncolors == 1,
      im = repmat(im,[1,1,3]);
    end
    if frame == firstframes(1),
      him = image([1,nc],[1,nr],im);
      axis image;
      axis([.5,x1(end)+.5,.5,y1(end)+.5]);
      axis off;
    else
      set(him,'cdata',im);
    end
    
    % draw frame number text box
    framestr = sprintf('Frame %d,frame',frame);
    if istimestamps,
      framestr = [framestr,sprintf(',t = %.2f s',timestamps(frame)-timestamps(1))];
    end
    if ~isempty(movietitle),
      framestr = {framestr,movietitle}; %#ok<AGROW>
    end
    % text doesn't show up in no display mode
    if titletext && isdisplay1,
      if frame == firstframes(1),
        htext = text(.5,.5,framestr,'Parent',hax,'BackgroundColor','w','Color','g','VerticalAlignment','bottom','interpreter','none');
      else
        set(htext,'String',framestr);
      end
    end
    % draw frame number text box
    if istimestamps,
      timestampstr = sprintf('%.2f s',timestamps(frame)-timestamps(1));
    else
      timestampstr = sprintf('Fr %d',frame);
    end
%     if ~isempty(movietitle),
%       timestampstr = {timestampstr,movietitle}; %#ok<AGROW>
%     end
    % text doesn't show up in no display mode
    if showtimestamps && isdisplay1,
      if frame == firstframes(1),
        htext = text(.5,hax.YLim(2)-(hax.YLim(2)/15),timestampstr,'Parent',hax,'BackgroundColor','k','Color','w','VerticalAlignment','bottom','interpreter','none');
      else
        set(htext,'String',timestampstr);
      end
    end
    
    % draw the zoomed image
    if dynamicflyselection && exist('fliesplotperframe','var'),
      zoomflies = reshape(fliesplotperframe(:,frame),[nzoomr,nzoomc]);
    end
    [tfhaspred,xy,tfocc] = apttrk.getPTrkFrame(frame);
    for i = 1:nzoomr,
      for j = 1:nzoomc,
        fly = zoomflies(i,j);
        
        % Is fly visible?
        if isnan(fly) || ~isalive(fly) ,
          is_fly_visible = false ;
        else
          is_fly_visible = tfhaspred(fly) && any(all(isfinite(xy(:,:,fly)),2),1);
        end
               
        if is_fly_visible ,
          % grab a box around (x,y)
          x = round(trx(fly).x(idx(fly)));
          y = round(trx(fly).y(idx(fly)));

          boxradx1 = min(boxradius,x-1);
          boxradx2 = min(boxradius,size(im,2)-x);
          boxrady1 = min(boxradius,y-1);
          boxrady2 = min(boxradius,size(im,1)-y);
          box = uint8(zeros(2*boxradius+1));
          box(boxradius+1-boxrady1:boxradius+1+boxrady2, boxradius+1-boxradx1:boxradius+1+boxradx2) = ...
            im(y-boxrady1:y+boxrady2,x-boxradx1:x+boxradx2);
          if frame == firstframes(1),
            himzoom(i,j) = image([x0(j),x1(j)],[y0(i),y1(i)],repmat(box,[1,1,3]));
          else
            set(himzoom(i,j),'cdata',repmat(box,[1,1,3]));
          end          
        else
          % If fly is not visible, show a placeholder
          if frame == firstframes(1),
            himzoom(i,j) = image([x0(j),x1(j)],[y0(i),y1(i)],repmat(uint8(123),[boxradius*2+1,boxradius*2+1,3]));
          else
            set(himzoom(i,j),'cdata',repmat(uint8(123),[boxradius*2+1,boxradius*2+1,3]));
          end
          continue
        end
        
        
      end
    end
    
    % plot the zoomed out position
    if frame == firstframes(1),
      for fly = 1:ntgts,
        
        if doplotcentertail,
          htailcenter(fly) = plot(nan,nan,'-','color',colors(fly,:));
        end
        if doplotapttail,
          htailapt(fly) = patch(nan,nan,[0,0,0],'FaceColor','none','EdgeColor',colors(fly,:),'EdgeAlpha',defaults.apttailalpha);
          %htailapt(fly) = plot(nan,nan,'-','color',colors(fly,:),'LineWidth',defaults.apttailwidth);
        end
        if doshowsex,
          hsexmarker(fly) = plot(nan,nan,'.','color',colors(fly,:),'markerfacecolor',colors(fly,:));
        end
        if doplotwings && doplotwings_perfly(fly),
          hwing(fly) = plot(nan,nan,'.-','color',colors(fly,:));
        end
        hlandmarks(fly) = plot(nan,nan,'.','color',colors(fly,:)*defaults.colormult,'MarkerSize',defaults.lmkmarkersize);
      end
    end
        
    for fly = 1:ntgts,
      if isalive(fly),
        i0 = max(1,idx(fly)-taillength);
        if doplotcentertail,
          set(htailcenter(fly),'xdata',trx(fly).x(i0:idx(fly)),...
            'ydata',trx(fly).y(i0:idx(fly)));
        end
        if doplotapttail,
          frame0 = max(1,frame-taillength);
          taillengthcurr = frame-frame0+1;
          [haspredtail,xytail] = apttrk.getPTrkFT(frame0:frame,fly);
          xytail(:,:,~haspredtail) = nan;
          xytail = permute(xytail,[3,1,2]);
          xytail = cat(1,xytail,nan([1,apttrk.npts,2]));
          xytail = reshape(xytail,[(taillengthcurr+1)*apttrk.npts,2]);
          set(htailapt(fly),'XData',xytail(:,1),'YData',xytail(:,2));
        end
        if doshowsex,
          sexi = find(strcmpi(trx(fly).sex{idx(fly)},sexes),1);
          sexmarker = sexmarkers{sexi};
          set(hsexmarker(fly),'xdata',trx(fly).x(idx(fly)),...
            'ydata',trx(fly).y(idx(fly)),...
            'color',colors(fly,:),...
            'marker',sexmarker,...
            'markerfacecolor',colors(fly,:));
        end
        set(hlandmarks(fly),'XData',xy(:,1,fly),'YData',xy(:,2,fly));
        if doplotwings && doplotwings_perfly(fly),
          xwing = [trx(fly).xwingl(idx(fly)),trx(fly).x(idx(fly)),trx(fly).xwingr(idx(fly))];
          ywing = [trx(fly).ywingl(idx(fly)),trx(fly).y(idx(fly)),trx(fly).ywingr(idx(fly))];
          set(hwing(fly),'XData',xwing,'YData',ywing);
        end
      else
        if doplotcentertail,
          set(htailcenter(fly),'xdata',[],'ydata',[]);
        end
        if doplotapttail,
          set(htailapt(fly),'xdata',[],'ydata',[]);
        end
        set(hlandmarks(fly),'xdata',[],'ydata',[]);
        if doshowsex,
          set(hsexmarker(fly),'xdata',[],'ydata',[]);
        end
        if doplotwings && doplotwings_perfly(fly),
          set(hwing(fly),'XData',[],'YData',[]);
        end
      end
    end
    
    % plot the zoomed views
    for i = 1:nzoomr,
      for j = 1:nzoomc,
        fly = zoomflies(i,j);
        if ~isnan(fly) && isalive(fly),
          x = trx(fly).x(idx(fly));
          y = trx(fly).y(idx(fly));
          offx = -round(x) + boxradius + .5;
          offy = -round(y) + boxradius + .5;
          x = x + offx; %boxradius + (x - round(x))+.5;
          y = y + offy; %boxradius + (y - round(y))+.5;
          x = x * scalefactor;
          y = y * scalefactor;
          x = x + x0(j) - 1;
          y = y + y0(i) - 1;
          a = trx(fly).a(idx(fly))*scalefactor;
          b = trx(fly).b(idx(fly))*scalefactor;
          theta = trx(fly).theta(idx(fly));
          
          lmkx = xy(:,1,fly)+offx;
          lmky = xy(:,2,fly)+offy;
          lmkx = lmkx * scalefactor;
          lmky = lmky * scalefactor;
          lmkx = lmkx + x0(j) - 1;
          lmky = lmky + y0(i) - 1;
          
          if doshowsex,
            s = sprintf('%d, %s',fly,trx(fly).sex{idx(fly)});
          else
            s = sprintf('%d',fly);
          end
          if doplotwings && doplotwings_perfly(fly),
            xwingl = trx(fly).xwingl(idx(fly)) - round(trx(fly).x(idx(fly))) + boxradius + .5;
            ywingl = trx(fly).ywingl(idx(fly)) - round(trx(fly).y(idx(fly))) + boxradius + .5;
            xwingl = xwingl * scalefactor;
            ywingl = ywingl * scalefactor;
            xwingl = xwingl + x0(j) - 1;
            ywingl = ywingl + y0(i) - 1;
            xwingr = trx(fly).xwingr(idx(fly)) - round(trx(fly).x(idx(fly))) + boxradius + .5;
            ywingr = trx(fly).ywingr(idx(fly)) - round(trx(fly).y(idx(fly))) + boxradius + .5;
            xwingr = xwingr * scalefactor;
            ywingr = ywingr * scalefactor;
            xwingr = xwingr + x0(j) - 1;
            ywingr = ywingr + y0(i) - 1;
            xwing = [xwingl,x,xwingr];
            ywing = [ywingl,y,ywingr];
          end

          if frame == firstframes(1),
            hzoom(i,j) = plot(lmkx,lmky,'.','color',colors(fly,:),'MarkerSize',defaults.lmkmarkersize);%drawflyo(x,y,theta,a,b);
            if doplotwings && doplotwings_perfly(fly),
              hzoomwing(i,j) = plot(xwing,ywing,'.-','color',colors(fly,:));
            else
              hzoomwing(i,j) = plot(nan,nan,'.-','color',colors(fly,:));
            end
            if isdisplay1,
              htextzoom(i,j) = text((x0(j)+x1(j))/2,.95*y0(i)+.05*y1(i),s,...
                'color',colors(fly,:),'horizontalalignment','center',...
                'verticalalignment','bottom','fontweight','bold');
            else
              if doshowsex,
                sexi = find(strcmpi(trx(fly).sex{idx(fly)},sexes),1);
                sexmarker = sexmarkers{sexi};
                htextzoom(i,j) = plot(x,y,'.',...
                  'color',colors(fly,:),'marker',sexmarker,'markerfacecolor',colors(fly,:));
              end
            end
          else
            set(hzoom(i,j),'XData',lmkx,'YData',lmky);
            if doplotwings && doplotwings_perfly(fly),
              set(hzoomwing(i,j),'XData',xwing,'YData',ywing,'Color',colors(fly,:));
            else
              set(hzoomwing(i,j),'XData',nan,'YData',nan,'Color',colors(fly,:));
            end
            if isdisplay1,
              set(htextzoom(i,j),'string',s,'color',colors(fly,:));
            else
              if doshowsex,
                sexi = find(strcmpi(trx(fly).sex{idx(fly)},sexes),1);
                sexmarker = sexmarkers{sexi};
                set(htextzoom(i,j),'xdata',x,...
                  'ydata',y,...
                  'color',colors(fly,:),...
                  'marker',sexmarker,...
                  'markerfacecolor',colors(fly,:));
              end
            end
          end
          set(hzoom(i,j),'color',colors(fly,:));
        else
          if frame == firstframes(1),
            hzoom(i,j) = plot(nan,nan,'-');
            hzoomwing(i,j) = plot(nan,nan,'.-');
            if isdisplay1,
              htextzoom(i,j) = text((x0(j)+x1(j))/2,.95*y0(i)+.05*y1(i),'',...
                'horizontalalignment','center',...
                'verticalalignment','bottom','fontweight','bold');
            else
              if doshowsex,
                htextzoom(i,j) = plot(nan,nan,'.','marker','none');
              end
            end
          else
            set(hzoom(i,j),'xdata',[],'ydata',[]);
            if doplotwings,
              set(hzoomwing(i,j),'XData',[],'YData',[]);
            end
            if isdisplay1,
              set(htextzoom(i,j),'string','');
            else
              if doshowsex,
                set(htextzoom(i,j),'xdata',[],'ydata',[]);
              end
            end
          end
        end
      end
    end
    
    if frame == firstframes(1),
      if ~isempty(figpos),
        set(fig,'Position',figpos);
      else
        input('Resize figure to the desired size, hit enter when done.');
        figpos = get(fig,'Position');
      end
      %set(fig,'visible','off');
      if useVideoWriter,
        if strcmpi(compression,'None') || strcmpi(compression,'Uncompressed AVI'),
          profile = 'Uncompressed AVI';
        else
          profile = 'Motion JPEG AVI';
        end
        aviobj = VideoWriter(aviname,profile); %#ok<TNMLP>
        set(aviobj,'FrameRate',fps);
        if ~strcmpi(profile,'Uncompressed AVI'),
          set(aviobj,'Quality',100);
        end
        open(aviobj);
      else
        if isempty(avifileTempDataFile),
          aviobj = avifile(aviname,'fps',fps,'quality',100,'compression',compression); 
        else
          aviobj = myavifile(aviname,'fps',fps,'quality',100,'compression',compression,...
            'TempDataFile',avifileTempDataFile); 
          fprintf('Temporary data file for avi writing: %s\n',aviobj.TempDataFile);
        end
      end
    end
    
    fr = getframe(hax);
    %fr = struct() ;
    %fr.cdata = uint8(randi(256, [428 600 3])-1) ;
    %fr.colormap = [] ;
    if frame == firstframes(1),
      height = size(fr.cdata,1);
      width = size(fr.cdata,2);
%       fr = getframe_invisible(hax);
%       [height,width,~] = size(fr);
      fprintf('Size of frame is %d x %d\n',height,width);
%       gfdata = getframe_initialize(hax);
%       [fr,height,width] = getframe_invisible_nocheck(gfdata,[height,width],false,false);
% 
%       height = ceil(height/4)*4;
%       width = ceil(width/4)*4;
%       fr = getframe_invisible(hax,[height,width]);
    else
      height1 = size(fr.cdata,1);
      width1 = size(fr.cdata,2);
      if height1 < height,
        dheight1 = floor((height-height1)/2);
        dheight2 = (height-height1)-dheight1;
        fr.cdata = fr.cdata(dheight1:end-dheight2,:,:);
      elseif height1 > height,
        dheight1 = floor((height1-height)/2);
        dheight2 = (height1-height)-dheight1;
        fr.cdata = cat(1,zeros([dheight1,width1,3],class(fr.cdata)),fr.cdata,zeros([dheight2,width1,3],class(fr.cdata)));
      end
      if width1 < width,
        dwidth1 = floor((width-width1)/2);
        dwidth2 = (width-width1)-dwidth1;
        fr.cdata = fr.cdata(:,dwidth1:end-dwidth2,:);
      elseif width1 > width,
        dwidth1 = floor((width1-width)/2);
        dwidth2 = (width1-width)-dwidth1;
        fr.cdata = cat(2,zeros([height,dwidth1,3],class(fr.cdata)),fr.cdata,zeros([height,dwidth2,3],class(fr.cdata)));
      end
      %fr = getframe_invisible_nocheck(gfdata,[height,width],false);
    end
    if useVideoWriter,
      writeVideo(aviobj,fr);
    else
      aviobj = addframe(aviobj,fr);
    end
    %set(fig,'Position',figpos);
  end
  
end
  
fprintf('Finishing AVI...\n');

if useVideoWriter,
  close(aviobj);
else
  aviobj = close(aviobj); %#ok<NASGU>
end
if fid > 0,
  fclose(fid);
end

fprintf('Cleanup...\n');

%getframe_cleanup(gfdata);

succeeded = true;