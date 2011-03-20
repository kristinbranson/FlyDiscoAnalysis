function trx = FlyBowlRegisterTrx(expdir,varargin)

[analysis_protocol,settingsdir,registrationparamsfilestr,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'registrationparamsfilestr','registration_params.txt',...
  'datalocparamsfilestr','dataloc_params.txt');

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
      
%% read in registration params

% name of parameters file
registrationparamsfile = fullfile(settingsdir,analysis_protocol,registrationparamsfilestr);
if ~exist(registrationparamsfile,'file'),
  error('Registration params file %s does not exist',registrationparamsfile);
end
% read
registration_params = ReadParams(registrationparamsfile);

%% detect registration marks

% name of annotation file
annfile = fullfile(expdir,dataloc_params.annfilestr);

% name of movie file
moviefile = fullfile(expdir,dataloc_params.moviefilestr);

% template filename should be relative to settings directory
if isfield(registration_params,'bowlMarkerType') && ...
    ~ismember(registration_params.bowlMarkerType,{'gradient'}),
  registration_params.bowlMarkerType = fullfile(settingsdir,analysis_protocol,registration_params.bowlMarkerType);
end

registration_params = struct2paramscell(registration_params);
% file to save image to
if isfield(dataloc_params,'registrationimagefilestr'),
  registration_params(end+1:end+2) = {'imsavename',fullfile(expdir,dataloc_params.registrationimagefilestr)};
end
% detect
try
  registration_data = detectRegistrationMarks(registration_params{:},'annName',annfile,'movieName',moviefile);
catch ME,
  error(['Error detecting registration marks:\n',getReport(ME)]);
end


%% apply registration

% name of input trx mat file
ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);

% name of movie
moviefile = fullfile(expdir,dataloc_params.moviefilestr);

% load trajectories
[trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname',annfile); %#ok<NASGU>
if ~succeeded,
  error('Could not load trajectories from file %s',ctraxfile);
end

% apply to trajectories
for fly = 1:length(trx),
  
  % apply transformation to 4 extremal points on the ellipse
  xnose0 = trx(fly).x + 2*trx(fly).a.*cos(trx(fly).theta);
  ynose0 = trx(fly).y + 2*trx(fly).a.*sin(trx(fly).theta);
  xtail0 = trx(fly).x - 2*trx(fly).a.*cos(trx(fly).theta);
  ytail0 = trx(fly).y - 2*trx(fly).a.*sin(trx(fly).theta);
  xleft0 = trx(fly).x + 2*trx(fly).b.*cos(trx(fly).theta-pi/2);
  yleft0 = trx(fly).y + 2*trx(fly).b.*sin(trx(fly).theta-pi/2);
  xright0 = trx(fly).x + 2*trx(fly).b.*cos(trx(fly).theta+pi/2);
  yright0 = trx(fly).y + 2*trx(fly).b.*sin(trx(fly).theta+pi/2);
  [xnose1,ynose1] = registration_data.registerfn(xnose0,ynose0);
  [xtail1,ytail1] = registration_data.registerfn(xtail0,ytail0);
  [xleft1,yleft1] = registration_data.registerfn(xleft0,yleft0);
  [xright1,yright1] = registration_data.registerfn(xright0,yright0);
  % compute the center as the mean of these four points
  x1 = (xnose1+xtail1+xleft1+xright1)/4;
  y1 = (ynose1+ytail1+yleft1+yright1)/4;
  % compute the major axis from the nose to tail distance
  a1 = sqrt( (xnose1-xtail1).^2 + (ynose1-ytail1).^2 ) / 4;
  % compute the minor axis length as the left to right distance
  b1 = sqrt( (xleft1-xright1).^2 + (yleft1-yright1).^2 ) / 4;
  % compute the orientation from the nose and tail points only
  theta1 = atan2(ynose1-ytail1,xnose1-xtail1);
  % store the registerd positions
  trx(fly).x_mm = x1;
  trx(fly).y_mm = y1;
  trx(fly).a_mm = a1;
  trx(fly).b_mm = b1;
  trx(fly).theta_mm = theta1;
  
  % add dt
  % TODO: fix timestamps after fix errors!
  if isfield(trx,'timestamps'),
    trx(fly).dt = diff(trx(fly).timestamps);
  else
    trx(fly).dt = repmat(1/trx(fly).fps,[1,trx(fly).nframes-1]);
  end

  trx(fly).fps = 1/mean(trx(fly).dt);
  trx(fly).pxpermm = 1 / registration_data.scale;
  
end

%% save registered trx to file

% name of output trx mat file
trxfile = fullfile(expdir,dataloc_params.trxfilestr);

save(trxfile,'trx','timestamps');

%% save params to mat file

registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
tmp = rmfield(registration_data,'registerfn'); %#ok<NASGU>
save(registrationmatfile,'-struct','tmp');

%% save params to text file
registrationtxtfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
fid = fopen(registrationtxtfile,'w');
fnssave = {'offX','offY','offTheta','scale','bowlMarkerTheta','featureStrengths',...
  'circleCenterX','circleCenterY','circleRadius'};
for i = 1:numel(fnssave),
  fn = fnssave{i};
  fprintf(fid,'%s,%f\n',fn,registration_data.(fn));
end
fclose(fid);
