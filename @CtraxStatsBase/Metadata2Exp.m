% [ns,expdirs] = obj.Metadata2Exp(fns,conditions)
% Find all experiments where the conditions specified by fns and
% conditions are satisfied. For each field in fns, we check if
% the corresponding element in conditions holds. conditions can either
% be a function handle, cell of possible char values, array of possible
% number values, single string, single number.

function [ns,expdirs] = Metadata2Exp(obj,fns,conditions)

% make arguments cells
if ischar(fns),
  if numel(conditions) > 1,
    error('Sizes of fields and conditions must match');
  end
  fns = {fns};
  if ~iscell(conditions),
    conditions = {conditions};
  end
end
if numel(conditions) ~= numel(fns),
  error('Sizes of fields and conditions must match');
end

% start with all exps included
idx = true(1,obj.nexpdirs);

% loop through metadata fields
for i = 1:numel(fns),

  fn = fns{i};
  condition = conditions{i};
  
  % get this field for all experiments
  val = obj.getMetaDataField(fn,'n',1:obj.nexpdirs);
  
  % make sure the values are in a cell
  if ~iscell(val),
    val = {val};
  end
  
  if strcmpi(class(condition),'function_handle'),
    % if condition is a function handle, apply it to val
    idx = idx & reshape(condition(val),size(idx));
  elseif iscell(condition),
    % if condition is a cell, then see if val is within condition
    if isnumeric(val),
      idx = idx & reshape(ismember(val,cell2mat(condition)),size(idx));
    else
      idx = idx & reshape(ismember(val,condition),size(idx));
    end
  elseif ischar(condition),
    % if condition is a string, compare with strcmp
    idx = idx & reshape(strcmp(val,condition),size(idx));
  elseif isnumeric(condition),
    % if condition is a number, then use ismember
    idx = idx & reshape(ismember(cell2mat(val),condition),size(idx));
  end
  
end

ns = find(idx);
expdirs = obj.expdir_bases(ns);