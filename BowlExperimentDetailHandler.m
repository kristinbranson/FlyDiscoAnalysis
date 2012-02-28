classdef BowlExperimentDetailHandler < OlyDat.ExperimentDetailHandler
%BoxExperimentDetailHandler Concrete ExperimentDetailHandler for Box assay
%   For use with OlyDat.Browser. At the moment, the experiment detail for
%   the box consists of two plots, a temperature trace and a vibration
%   trace.

    properties (Hidden)
      videodiagnostics_params = {}; 
      compute_extra_diagnostics = false;
    end
    
    properties (Dependent)
      fDetailHasOpened
    end
    
    methods
        function tf = get.fDetailHasOpened(obj)
            tf = false;
        end
    end
    
    methods
      
      function SetVideoDiagnosticsParams(obj,varargin)
        obj.videodiagnostics_params = varargin;
      end
      
      function SetComputeExtraDiagnostics(obj,val)
        obj.compute_extra_diagnostics = val;
      end
        
        function open(obj,data,varargin)
          
          if isempty(data)
            return;
          end
          
          if ~isfield(data,'file_system_path'),
            return;
          end
          fprintf('Path: %s\n',data.file_system_path);
          hwait = nan;
          if exist(data.file_system_path,'dir'),
            hwait = mywaitbar(0,hwait,'Opening experiment directory...');
            if ispc,
              winopen(data.file_system_path);
            else
              web(data.file_system_path,'-browser');
            end
            drawnow;
          else
            uiwait(warndlg(sprintf('Experiment %s does not exist',data.file_system_path)));
          end

          if obj.compute_extra_diagnostics && ~exist(fullfile(data.file_system_path,'video_diagnostics.png'),'file'),
            hwait = mywaitbar(.1,hwait,'Computing video diagnostics...');
            VideoDiagnostics(data.file_system_path,'debug',true,obj.videodiagnostics_params{:});
            drawnow;
          end
          moviename = fullfile(data.file_system_path,'movie.ufmf');
          if obj.compute_extra_diagnostics && exist(moviename,'file'),
            hwait = mywaitbar(.9,hwait,'Showing UFMF...');
            showufmf('UFMFName',moviename);
            drawnow;
          end
          if ishandle(hwait), delete(hwait); end
          
        end
                
        function refresh(obj,data)

          if isempty(data)
            return;
          end
          fprintf('Experiment: %s\n',data.experiment_name(numel('FlyBowl_')+1:end));

        end
        
        function close(obj)
        end
        
    end
    
    
end