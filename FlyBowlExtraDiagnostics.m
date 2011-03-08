function FlyBowlExtraDiagnostics(expdir,varargin)

%% compute temperature diagnostics
TemperatureDiagnostics(expdir,varargin{:});

%% compute background model diagnostics
BkgdModelDiagnostics(expdir,varargin{:});

%% compute bias diagnostics
BowlBiasDiagnostics(expdir,varargin{:});