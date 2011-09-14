classdef BowlExperimentDetailHandler < OlyDat.ExperimentDetailHandler
%BoxExperimentDetailHandler Concrete ExperimentDetailHandler for Box assay
%   For use with OlyDat.Browser. At the moment, the experiment detail for
%   the box consists of two plots, a temperature trace and a vibration
%   trace.

    properties (Hidden)
    end
    
    properties (Dependent)
      fDetailHasOpened
    end
    
    properties (Public)
    end
    
    methods
        function tf = get.fDetailHasOpened(obj)
            tf = false;
        end
    end
    
    methods
        
        function open(obj,data)
          if ~isempty(data)
            obj.refresh(data);
          end
        end
                
        function refresh(obj,data)
          if isempty(data),
            return;
          end
          fprintf('Experiment: %s\n',data.experiment_name(numel('FlyBowl_')+1:end));
          if ~isfield(data,'file_system_path'),
            return;
          end
          fprintf('Path: %s\n',data.file_system_path);
          if ~exist(fullfile(data.file_system_path,'video_diagnostics.png'),'file'),
            VideoDiagnostics(data.file_system_path,'debug',true);
          end
          moviename = fullfile(data.file_system_path,'movie.ufmf');
          if exist(moviename,'file'),
            showufmf('UFMFName',moviename);
          end
          if ispc,
            winopen(data.file_system_path);
          else
            web(data.file_system_path,'-browser');
          end
        end
        
        function close(obj)
        end
        
    end
    
    
end