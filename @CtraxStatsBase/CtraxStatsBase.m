classdef CtraxStatsBase < handle
  properties
            
    % base name of movie
    moviefilestr = 'movie.ufmf';
    
    % base name of ann file
    annfilestr = 'movie.ufmf.ann';
    
    % base name of mat file
    ctraxfilestr = 'ctrax_results.mat';
    
    % after registering
    trxfilestr = 'registered_trx.mat';
    
    % landmark-based measurements
    landmarksfilestr = 'derived_landmarks_trx.mat';
    % closest fly-based measurements
    closestflyfilestr = 'derived_closestfly_trx.mat';
    % speed-based measurements
    speedfilestr = 'derived_speed_trx.mat';
    
    % root directory containing all experiment directories
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    
    % name of file to write registration results to
    registrationfilestr = 'registration_data.mat';

    % parameters for registration
    detectregistrationparams = struct;
    
    % if NOWRITEACCESS, we aren't able to write to the experiment directories
    % temporarily, data will be in resultsdir
    NOWRITEACCESS = true;
    resultsdir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118';
    
    % default arena center and radius
    arena_center_mm = [0,0];
    arena_radius_mm = 127/2;
    
    % default field of view for computing angle subtended
    fov = pi;
    
    % smoothing orientation
    thetafil = [1 4 6 4 1]/16;
    
    % HistogramMeasurements parameters
    histogrammeasurements_nbins = 100;
    histogrammeasurements_doplot = true;
    histogrammeasurements_averaging = 'perexp_perfly';
    histogrammeasurements_jackknife = '';
    histogrammeasurements_axesstyleparams = {'TickDir','out'};
    histogrammeasurements_histstyleparams = {'linestyle','-','linewidth',4,'marker','.','color',zeros(1,3)};
    histogrammeasurements_perflyhiststyleparams = {'linestyle','-','linewidth',.5,'marker','none'};
    histogrammeasurements_perexphiststyleparams = {'linestyle','-','linewidth',1,'marker','none'};
    histogrammeasurements_errstyleparams = {'facecolor',.8*ones(1,3),'edgecolor','none'};
    histogrammeasurements_legendstyleparams = {};
    histogrammeasurements_ploterrorbars = 'stderr';
    histogrammeasurements_plotperfly = false;
    histogrammeasurements_plotperexp = false;
    histogrammeasurements_docomputestd = true;
    histogrammeasurements_docomputestderr = true;
    histogrammeasurements_figpos = [1,1,1200,600];
    histogrammeasurements_lim_prctile = [1,99];
    
  end
  
  properties (SetAccess = private, GetAccess = public)

    % number of expdirs open
    nexpdirs = 0;
    
    % open experiment directories
    expdirs = {}; % full path for reading movie
    write_expdirs = {}; % full path for reading trx, writing
    expdir_bases = {}; % path within rootdatadir
    
    % current file names
    moviefiles = {};
    annfiles = {};
    ctraxfiles = {};
    trxfiles = {};
    landmarksfiles = {};
    closestflyfiles = {};
    speedfiles = {};
    registrationfiles = {};
    
    % function handle for reading frames from current movie
    readframes = {};
    
    % number of frames in current movies
    nframes = [];
    
    % file pointer for current movie
    moviefids = [];
    
    % header info read from current movie
    headerinfos = {};

    % movie frame size
    nrs = [];
    ncs = [];
    ncolors = [];
    width_mms = [];
    height_mms = [];

    % annotation info
    anns = {};
    
    % current trajectories
    trx = [];
    % units of fields of struct
    units = struct('x',parseunits('px'),...
      'y',parseunits('px'),...
      'theta',parseunits('rad'),...
      'a',parseunits('px'),...
      'b',parseunits('px'),...
      'theta_mm',parseunits('rad'),...
      'x_mm',parseunits('mm'),...
      'y_mm',parseunits('mm'),...
      'a_mm',parseunits('mm'),...
      'b_mm',parseunits('mm'),...
      'timestamps',parseunits('s'),...
      'dt',parseunits('s'));
    
    % number of flies in each movie
    nfliespermovie = [];
    % total number of flies
    nflies = 0;
    
    % which flies are in a given movie
    movie2flies = {};
    
    % which movie each fly corresponds to 
    fly2movie = [];
    
    % which extra measurements have been computed
    didComputeLandmarkMeasurements = false;
    didComputeClosestFlyMeasurements = false;
    didComputeSpeedMeasurements = false;
  end
  
  properties (Hidden = true)
    
  end
  
  properties (Constant)
    
    % names of images read from annotation file
    annfile_images = {'background_median','background_mean',...
      'background_mad','background_std',...
      'hfnorm','fracframesisback',...
      'background_center','background_dev'};
    

  
  end
  
  methods
    
    % obj = CtraxDiagnostics('param',value,...)
    %
    % Constructor:
    % For each optional parameter 'param', the member variable 'param' is
    % set to value.
    % Experiment directory is set to expdir, relevant files are opened
    % and processed.
    function obj = CtraxStatsBase(varargin)
      
      % optional parameters
      for i = 1:2:length(varargin)-1,
        obj.(varargin{i}) = varargin{i+1};
      end
    end
    
    % deconstructor
    function delete(obj)
      
      obj.RemoveAllExpDirs()
      
    end
    
    % obj.AddExpDir(expdir)
    % Adds the trajectories and all other data associated with experiment 
    % expdir to the data represented by obj. The movie associated with
    % expdir is also opened for reading. If the trajectories for expdir
    % have not yet been registered, then registration is performed and the
    % results are saved to the trxfile. If expdir is already represented,
    % it is removed then added again.  
    AddExpDir(obj,expdir)
    
    % obj.RemoveExpDir(expdir)
    % Removes the trajectories and all other data associated with experiment 
    % expdir to the data represented by obj. The movie file associated with
    % expdir is also closed. 
    RemoveExpDir(obj,expdir)
    
    % obj.RemoveAllExpDirs()
    % Removes all the data stored in obj, closes the relevant file handles.
    RemoveAllExpDirs(obj)
    
    % obj.DeleteRegistrationFiles([expdirs])
    % Deletes the files created during the registration step -- the
    % registerd trajectories (<expdir>/<trxfilestr>) and the registration
    % parameters (<expdir>/<trxfilestr>). If expdirs is input, then the
    % registration files for the specified expdirs will be deleted. If not
    % specified, then registration files for all loaded experiments will be
    % deleted.
    DeleteRegistrationFiles(obj,expdirs)
    
    % obj.ComputeLandmarkMeasurements()
    % Compute per-frame measurements of distances, angles to landmarks in
    % the arena, as well as speeds related to these coordinate systems for
    % all currently loaded trajectories. If new experiment directories are
    % added after ComputeLandmarkMeasurements is called, landmark-derived
    % measurements will be computed for the new trajectories when the
    % experiment is added. Currently, the only landmark is the circular
    % arena wall, which is assumed to have an origin arena_center_mm and
    % radius arena_radius_mm in the registered coordinate system. The
    % landmark-based derived measurements will be stored to
    % <expdir>/<landmarksfilestr>.
    ComputeLandmarkMeasurements(obj)

    % obj.DeleteLandmarkMeasurementFiles([expdirs])
    % Deletes the files created during the ComputeLandmarkMeasurements step
    % -- the file containing the landmark-based derived measurements
    % <expdir>/<landmarksfilestr>. If expdirs is input, then the
    % landmark files for the specified expdirs will be deleted. If not
    % specified, then landmark files for all loaded experiments will be
    % deleted.
    DeleteLandmarkMeasurementFiles(obj,expdirs)
    
    % obj.ComputeClosestFlyMeasurements()
    % Compute per-frame measurements of distances, angles to the nearest fly
    % to each fly, as well as speeds related to these coordinate systems
    % for all currently loaded trajectories. If new experiment directories
    % are added after ComputeClosestFlyMeasurements is called, closest
    % fly-derived measurements will be computed for the new trajectories
    % when the experiment is added. The closest fly-based derived
    % measurements will be stored to <expdir>/<closestflyfilestr>. 
    ComputeClosestFlyMeasurements(obj)
    
    % obj.DeleteClosestFlyMeasurementFiles([expdirs])
    % Deletes the files created during the ComputeClosestFlyMeasurements step
    % -- the file containing the closest fly-based derived measurements
    % <expdir>/<closestflyfilestr>. If expdirs is input, then the
    % closest-fly files for the specified expdirs will be deleted. If not
    % specified, then closest fly files for all loaded experiments will be
    % deleted.
    DeleteClosestFlyMeasurementFiles(obj,expdirs)

    % obj.ComputeSpeedMeasurements()
    % Compute per-frame measurements of speeds and accelerations for all
    % currently loaded trajectories. If new experiment directories are
    % added after ComputeSpeedMeasurements is called, speeds will
    % be computed for the new trajectories when the experiment is added.
    % The speed measurements will be stored to <expdir>/<speedfilestr>. 
    ComputeSpeedMeasurements(obj)
    
    % obj.DeleteSpeedMeasurementFiles([expdirs])
    % Deletes the files created during the ComputeSpeedMeasurements step --
    % the file containing the speed measurements <expdir>/<speedfilestr>.
    % If expdirs is input, then the speed files for the specified expdirs
    % will be deleted. If not specified, then speed files for all loaded
    % experiments will be deleted.
    DeleteSpeedMeasurementFiles(obj,expdirs)
    
  end
  
  methods (Static)
    
  end
  
end