function  bioToggleListenerState(hList,state)
% Sets 'on' or 'off' the Enabled property of a listener or a set of
% listeners grouped either in a cell or array.

%  Copyright 2010 The MathWorks, Inc.


if islogical(state)
    logicalState = state(1);
else
    logicalState = strcmpi(state,'on');
end
if logicalState
    state = 'on';
else
    state = 'off';
end
    
if feature('HGUsingMATLABClasses')
    if iscell(hList)
        for i = 1:numel(hList)
            listener = hList{i};
            if isa(listener,'handle.listener')
                set(listener,'Enabled', state); 
            else
                listener.Enabled = logicalState;
            end
        end
    else
        for i = 1:numel(hList)
            listener = hList(i);
            if isa(listener,'handle.listener')
                set(listener,'Enabled', state); 
            else
                listener.Enabled = logicalState;
            end
        end
    end
else
    if iscell(hList)
        for i = 1:numel(hList)
            listener = hList{i};
            set(listener,'Enabled', state);
        end
    else
        set(hList,'Enabled',state);
    end
end
