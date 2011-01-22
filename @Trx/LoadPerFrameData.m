function data = LoadPerFrameData(obj,fn,n)

filename = obj.GetPerFrameFile(fn,n);
if exist(filename,'file')
  x = load(filename,'data','units');
else
  [x.data,x.units] = obj.ComputerPerFrameData(fn,n);
end
obj.units.(fn) = x.units;
obj.SetPerFrameData(fn,x.data,n);
data = x.data;