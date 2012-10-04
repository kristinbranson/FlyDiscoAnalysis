function FlyBowlExtraDiagnostics(expdir,varargin)

[dotemperature,dobkgd,dobias,dovideo,leftovers] = ...
  myparse_nocheck(varargin,'dotemperature',true,...
  'dobkgd',true,...
  'dobias',true,...
  'dovideo',true);
if ischar(dotemperature),
  dotemperature = str2double(dotemperature)~=0;
end
if ischar(dobkgd),
  dobkgd = str2double(dotemperature)~=0;
end
if ischar(dobias),
  dobias = str2double(dotemperature)~=0;
end
if ischar(dovideo),
  dovideo = str2double(dotemperature)~=0;
end

%% compute temperature diagnostics
if dotemperature,
  TemperatureDiagnostics(expdir,leftovers{:});
end

%% compute background model diagnostics
if dobkgd,
  BkgdModelDiagnostics(expdir,leftovers{:});
end

%% compute bias diagnostics
if dobias,
  BowlBiasDiagnostics(expdir,leftovers{:});
end

%% compute average frame rate from first and last timestamps and number of frames
if dovideo,
  VideoDiagnostics(expdir,leftovers{:});
end

%% 

%%

if isdeployed,
  close all;
end