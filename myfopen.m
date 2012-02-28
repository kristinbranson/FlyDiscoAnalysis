function fid = myfopen(filename,varargin)

if isempty(filename),
  fid = 1;
else
  fid = fopen(filename,varargin{:});
  if fid <= 0,
    error('Could not open file %s',filename);
  end
end