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
  
  properties
    
  end
  
  properties (SetAccess = private, GetAccess = public)
    
    % number of trajectories
    ntrx = 0;
    
    % number of experiment directories
    nexpdirs = 0;
    
    % names of the experiment directories, relative to root data directory
    expdir_bases = {};
    
    % full path to experiment directories in rootreaddir
    expdir_reads = {};
    
    % full path to experiment directories in rootwritedir
    expdir_writes = {};
    
    % cached data
    data = {};
    
    % parameters
    
    % directory within which experiment movies are contained
    rootreaddir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    
    % directory within which experiment tracking & analysis data are
    % contained
    rootwritedir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118';
    
    % derived, per-frame measurement directory, relative to experiment directory
    perframedir = 'perframe';
    
    % base name of mat file
    ctraxfilestr = 'ctrax_results.mat';
    
    % after registering
    trxfilestr = 'registered_trx.mat';
    
    % history of per-frame properties loaded & their last access time
    perframehistory = cell(0,2);
    
    % maximum number of doubles in data
    maxdatacache = 2^27;
    
    % number of doubles currently stored
    ndatacached = 0;
    
    % data cached for each experiment
    datacached = {};
    
    % parameters for registration
    detectregistrationparams = struct;
    
    % where the bowl marker should be, based on bowl
    %bowl2MarkerAngle = 3*pi/4;
        
    % default arena center and radius
    %arena_center_mm = [0,0];
    %arena_radius_mm = 127/2;
    
    % default field of view for computing angle subtended
    fov = pi;
    
    % smoothing orientation
    thetafil = [1 4 6 4 1]/16;
    
  end
  
  properties (Hidden = true)
    
  end
  
  properties (Constant)
    
  end
  
  methods
    
    function obj = Trx(varargin)

      % all arguments are parameters
      for i = 1:2:nargin-1,
        set(obj,varargin{i},varargin{i+1});
      end
      
    end
    
    % deconstructor
    function delete(obj)

      obj.ntrx = 0;
      obj.expdir_bases = {};
      obj.expdir_reads = {};
      obj.expdir_writes = {};
      obj.data = {};
      obj.perframehistory = cell(0,2);
      obj.ndatacached = 0;
      obj.datacached = {};
      
    end
        
  end
  
end