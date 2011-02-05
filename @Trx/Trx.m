classdef Trx < handle

  % Trx class
  %
  % The main data consists of any number of experiments (each corresponding
  % to one video of multiple flies). 
  % ROOTREADDIR: There is a root directory ROOTREADDIR that is the root
  % directory for all data we only need read access to -- currently the
  % video and metadata file (someday, when Ctrax is part of the pipeline,
  % the trajectories output by Ctrax will be included here). 
  % ROOTWRITEDIR: There is a root directory ROOTWRITEDIR that contains all
  % data for which we need write access; data that can be created by the
  % FlyBowlAnalysis code. This is currently the registered trajectories, sex 
  % classification, per-frame derived behavior statistics. 
  % EXPDIR_BASES: Within each of these root directories, there will be a
  % subdirectory of the same name for each experiment (name stored in
  % EXPDIR_BASES). 
  % CTRAXFILESTR: For each experiment subdirectory of ROOTWRITEDIR
  % (currently), there will be a mat file containing the trajectories
  % output by Ctrax. The name of this file within the experiment
  % subdirectory is defined by CTRAXFILESTR. 
  % REGISTER: The trajectories can be registered (pixels converted to
  % milimeters, coordinate system translated and rotated so that they align
  % between experiments) with the function REGISTER. 
  % TRXFILESTR: The registered trajectories are stored in the file
  % TRXFILESTR within the experiment directory in ROOTWRITEDIR. 
  % DETECTREGISTRATIONPARAMS: The parameters for REGISTER are stored in the
  % struct DETECTREGISTRATIONPARAMS. 
  % PERFRAMEDIR: The derived, per-frame time series for a given experiment
  % are stored within the subdirectory PERFRAMEDIR of each experiment
  % subdirectory in ROOTWRITEDIR. Within PERFRAMEDIR, there is a mat file
  % for each type of per-frame derived measurement, e.g. speed, change in
  % orientation. These mat files each contain the time series for all flies
  % in the given experiment, and are created or loaded on demand. 
  % PERFRAMEHISTORY: For each per-frame derived measure, the last access
  % time for any experiment is stored in the NFNS x 2 PERFRAMEHISTORY. 
  % MAXDATACACHED: When the total number of doubles stored in data exceeds
  % MAXDATACACHED, the least recently used per-frame derived measurement is
  % cleared from the memory cache. 

  
  properties (SetAccess = protected, GetAccess = public)

    % number of flies
    nflies = 0;
    
    % number of experiment directories
    nexpdirs = 0;
    
    % number of flies per movie
    nfliespermovie = [];
    
    %% parameters of data locations
    
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
    analysis_protocol = 'current';
    
    datalocparamsfilestr = 'dataloc_params.txt';
    
    dataloc_params = [];

    % names of the experiment directories, relative to root data directory
    expdir_bases = {};
    
    % full path to experiment directories in rootreaddir
    expdir_reads = {};
    
    % full path to experiment directories in rootwritedir
    expdir_writes = {};
    
    % names of mat files containing registered, processed trajectories,
    % with sex classified
    trxfiles = {};

    % movie locations
    movienames = {};
    
    %% landmark parameters

    landmark_params = [];
    
    %% per-frame parameters
    perframe_params = [];
    % default field of view for computing angle subtended
    %fov = pi;
    
    % smoothing orientation
    %thetafil = [1 4 6 4 1]/16;
    
    % units of per-frame properties
    units = struct;

    
    %% sex classifier parameters
    sexclassifier_params = [];
        
    %% data caching
    
    % history of per-frame properties loaded & their last access time
    perframehistory = cell(0,2);
    
    % maximum number of doubles in data
    maxdatacached = 2^27;
    
    % number of doubles currently stored
    ndatacached = 0;
    
    % per experiment
    ndatacachedperexp = [];
    
    % data cached for each experiment
    datacached = {};
    
    %% video info
    
    % number of frames in the video
    movie_nframes = [];
    
    % image sizes
    nrs = [];
    ncs = [];
    ncolors = [];
    
    % movie size in mm
    width_mms = [];
    height_mms = [];
    
    pxpermm = [];
    fps = [];
    
    %% trajectory frames
    
    firstframes = [];
    endframes = [];
    nframes = [];

    %% indexing stuff
    
    exp2flies = {};
    fly2exp = [];
    
  end
  
  properties (Hidden = true)
    
  end
  
  properties (Constant)
    
  end
  
  methods
    
    function obj = Trx(varargin)

      % all arguments are parameters
      for i = 1:2:nargin-1,
        obj.(varargin{i}) = varargin{i+1};
      end
      
      obj.ReadAllParams();
      
    end
        
    % deconstructor
    function delete(obj)

      obj.nflies = 0;
      obj.expdir_bases = {};
      obj.expdir_reads = {};
      obj.expdir_writes = {};
      obj.perframehistory = cell(0,2);
      obj.ndatacached = 0;
      obj.datacached = {};
      
    end
    
    function flyidx = getFlyIdx(obj,n,fly)
      flyidx = obj.exp2flies{n}(fly);
    end
    
    function [n,fly] = getExpFly(obj,flyidx)
      n = obj.fly2exp(flyidx);
      fly = nan(size(flyidx));
      for i = 1:numel(flyidx),
        flycurr = find(obj.exp2flies{n} == flyidx(i),1);
        if isempty(flycurr), 
          error('Sanity check: flyidx %d mapped to exp %d, but this is not in exp2flies{%d}',flyidx(i),n,n);
        end
        fly(i) = flycurr;
      end
    end

%     function n = numel(obj)
%       
%       n = obj.nflies;
%       
%     end
% 
%     function varargout = size(obj)
% 
%       n = obj.nflies;
%       if nargout <= 1,
%         varargout = {n};
%       else
%         varargout = num2cell([1,n,ones(1,nargout-2)]);
%       end
%       
%     end

    function [data,units] = ComputePerFrameData(obj,fn,n)
      
      funname = sprintf('compute_%s',fn);
      [data,units] = feval(funname,obj,n);
      filename = obj.GetPerFrameFile(fn,n);
      save(filename,'data','units');
      
    end

    %
    
    % function declarations

    % read parameters
    ReadAllParams(obj,varargin)
    ReadLandmarkParams(obj)
    ReadPerFrameParams(obj)
    ReadSexClassifierParams(obj)
        
    % AddExpDir(obj,expdir,vidinfo)
    % add a new experiment directory
    AddExpDir(obj,expdir,vidinfo)
    
    % RemoveExpDir(obj,expdir)
    % remove an experiment directory
    RemoveExpDir(obj,expdir)
        
    % StoreTrajectories(obj,n,traj)
    % store trajectories traj for experiment n
    StoreTrajectories(obj,n,traj)
    
    FreeDataCache(obj,ndataadd)
    
    [varargout] = subsref(obj,s)
    
    [varargout] = GetPerFrameData(obj,fn,varargin)
    
    SetPerFrameData(obj,fn,x,varargin)
    
    UpdatePerFrameAccessTime(obj,fn)
    
    x = LoadPerFrameData(obj,fn,n)
    
    CleanPerFrameData(obj,fn,varargin)
    
  end
  
  methods(Static)
    
    fns = PerFrameFieldNames()
    
  end
  
end