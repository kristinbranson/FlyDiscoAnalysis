classdef SingleStatTubeAveragedPlot_FlyBowl < SingleStatTubeAveragedPlot
% virtual class: extends the SingleStatTubeAveragedPlot in ways common for
% all FlyBowl plotting
    
    properties
        FailureFlags = {'manual_pf',{'F'},'automated_pf',{'F'},'flag_redo',1};
    end
    
    methods
      function tfFail = detectFailure(obj,data)
        
        tfFail = false(1,numel(data));
        for i = 1:2:numel(obj.FailureFlags)-1,
          fn = obj.FailureFlags{i};
          val = obj.FailureFlags{i+1};
          if isnumeric(val),
            tfFail = tfFail | ismember([data.(fn)],val);
          else
            tfFail = tfFail | ismember({data.(fn)},val);
          end
        end
        
      end
    end
    
end