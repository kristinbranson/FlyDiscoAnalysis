classdef Test < handle

  properties
    
    x = [];
    
  end
  
  methods
    
    function obj = Test(x)
      obj.x = x;
    end
    
    function n = numel(obj)
      n = numel(obj.x);
    end
    
    function varargout = size(obj)
      
      varargout = cell(1,max(1,nargout));
      [varargout{:}] = size(obj.x);
    
    end
    
    function varargout = subsref(obj,s)

      if (s(1).type == '.') && ismember(s(1).subs,methods(obj)),
        [varargout{1:nargout}] = builtin('subsref',obj,s);
        return;
      end
      
      varargout = num2cell(repmat(123,[1,max(1,nargout)]));
    end
    
    function res = test(obj)
      res = obj.x + 234;
    end
    
    function res = testtest(obj)
      res = testtesttest(obj);
    end
    
  end
  
end
