function FlyBowlExtraDiagnostics(expdir,varargin)

%% compute temperature diagnostics
TemperatureDiagnostics(expdir,varargin{:});

%% compute background model diagnostics
BkgdModelDiagnostics(expdir,varargin{:});

%% compute bias diagnostics
BowlBiasDiagnostics(expdir,varargin{:});

%% compute average frame rate from first and last timestamps and number of frames
VideoDiagnostics(expdir,varargin{:});

%% 

%%

if isdeployed,
  close all;
end