function result = loadSingleVariableAnonymously(file_name, variable_name) 
% Load a single variable from a .mat file, but just return it, don't assign it
% to a variable.  Handy if you want to rename the variable in the .mat when
% you load it.

load('-mat', file_name, variable_name) ;
result = eval(variable_name) ;

end

