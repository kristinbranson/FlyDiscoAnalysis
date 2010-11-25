classdef CtraxDiagnostics
  properties
    
    % current experiment directory
    expdir = '';
    
    % base name of movie
    moviefilestr = 'movie.ufmf';
    
    % base name of ann file
    annfilestr = 'movie.ufmf.ann';
    
    % base name of mat file
    matfilestr = 'movie.mat';
    
    % if NOWRITEACCESS, we aren't able to write to the experiment directories
    % temporarily, data will be in resultsdir
    NOWRITEACCESS = true;
    resultsdir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118';
    
    % current file names
    moviefile = '';
    annfile = '';
    matfile = '';
    
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
    function obj = CtraxDiagnostics(expdir)
      obj.isopen = false;
      if nargin >= 1,
        obj.SetExpDir(expdir);
      end
    end
    function delete(obj)
      % delete(obj)
      % clear obj
      % Closes the SerialPort if open before deallocating. 
      if obj.IsOpen,
        fprintf('Closing port %s.\n',obj.Port);
        fclose(obj.SerialPort);
        obj.IsOpen = false;
        obj.SerialPort = [];
      end
    end
    function set.Port(obj,value)
      if obj.IsOpen,
        error('Cannot set Port while SerialPort is open');
      end
      obj.Port = value;
    end
    function set.BaudRate(obj,value)
      obj.BaudRate = value;
      if obj.IsOpen,
        set(obj.SerialPort,'BaudRate',value);
      end
    end
    function set.DataBits(obj,value)
      obj.DataBits = value;
      if obj.IsOpen,
        set(obj.SerialPort,'DataBits',value);
      end
    end
    function set.Parity(obj,value)
      obj.Parity = value;
      if obj.IsOpen,
        set(obj.SerialPort,'Parity',value);
      end
    end
    function set.StopBits(obj,value)
      obj.StopBits = value;
      if obj.IsOpen,
        set(obj.SerialPort,'StopBits',value);
      end
    end
    function set.FlowControl(obj,value)
      obj.FlowControl = value;
      if obj.IsOpen,
        set(obj.SerialPort,'FlowControl',value);
      end
    end
    function set.Terminator(obj,value)
      obj.Terminator = value;
      if obj.IsOpen,
        set(obj.SerialPort,'Terminator',value);
      end
    end
    function set.ReadAsyncMode(obj,value)
      obj.ReadAsyncMode = value;
      if obj.IsOpen,
        set(obj.SerialPort,'ReadAsyncMode',value);
      end
    end
    [success,errormsg] = open(obj)
    [success,errormsg] = close(obj)
    [temp,humid,success,errormsg] = read(obj,nReadings)
    [success,errormsg] = flush(obj)
  end
  methods(Static)
    portsAvailable = getAvailableSerialPorts()
  end
end