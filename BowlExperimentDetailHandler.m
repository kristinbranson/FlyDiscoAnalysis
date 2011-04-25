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
          if ~isempty(data) && isfield(data,'file_system_path'),
            if ispc,
              winopen(data.file_system_path);
            else
              web(data.file_system_path,'-browser');
            end
          end
        end
        
        function close(obj)

        end
        
    end
    
    
end