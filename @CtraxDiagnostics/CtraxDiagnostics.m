classdef CtraxDiagnostics
  properties
    
    % current experiment directory
    expdir = '';
    
    % base name of movie
    moviefilestr = 'movie.ufmf';
    
    % base name of ann file
    annfilestr = 'movie.ufmf.ann';
    
    % base name of mat file
    ctraxfilestr = 'ctrax_results.mat';
    
    % after calibrating
    trxfilestr = 'calibrated_trx.mat';
    
    % if NOWRITEACCESS, we aren't able to write to the experiment directories
    % temporarily, data will be in resultsdir
    NOWRITEACCESS = true;
    resultsdir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118';
    
    % current file names
    moviefile = '';
    annfile = '';
    ctraxfile = '';
    trxfile = '';
    
    % fid of movie 
    moviefid = -1;
    
    % whether an experiment directory is currently open
    isopen = false;
    
    % names of images read from annotation file
    annfile_images = {'background_median','background_mean',...
      'background_mad','background_std',...
      'hfnorm','fracframesisback',...
      'background_center','background_dev'};

  end
  
  methods
    
    function obj = CtraxDiagnostics(expdir,varargin)
      % optional parameters
      if nargin >= 2,
        for i = 1:2:length(varargin)-1,
          obj.(varargin{i}) = varargin{i+1};
        end
      end
      % open input expdir
      obj.isopen = false;
      if nargin >= 1,
        obj.SetExpDir(expdir);
      end
    end
    
  end
  
end