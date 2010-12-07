classdef CtraxDiagnostics < CtraxStatsBase

  properties
            
  end
  
  properties (SetAccess = private, GetAccess = public)
    
  end
  
  properties (Hidden = true)
    
  end
  
  properties (Constant)
      
  end
  
  methods
    
    function obj = CtraxDiagnostics(varargin)
      % obj = CtraxDiagnostics('param',value,...)
      % 
      % Constructor:
      % For each optional parameter 'param', the member variable 'param' is
      % set to value. 
      % Experiment directory is set to expdir, relevant files are opened
      % and processed. 
      
      obj@CtraxStatsBase(varargin{:});
      
    end  
    
    heatmap = CenterPositionHeatmap(obj,varargin)
    heatmap = CenterPositionBinEntries(obj,varargin)
    heatmap = CenterPositionHeatmapPolar(obj,varargin)
    heatmap = CenterPositionBinEntriesPolar(obj,varargin)
    
    heatmap = BowlAngleBias(obj,varargin)
    
  end
  
end