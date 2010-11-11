function simple_behavior_statistics(expdir,varargin)

[matfilestr] = myparse(varargin,'movie.mat');

matname = fullfile(expdir,matfilestr);
if ~exist(matname,'file'),
  error('Mat file %s does not exist',matname);
end