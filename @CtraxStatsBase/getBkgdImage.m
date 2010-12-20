function [bkgdImage,xdata,ydata] = getBkgdImage(obj,varargin)

[expdirs,ns] = myparse(varargin,'expdir',{},'n',[]);

if isempty(ns),
  ns = obj.expdir2n(expdirs);
end

bkgdImage = cell(size(ns));
for i = 1:numel(ns),
  n = ns(i);
  bkgdImage0 = obj.anns{n}.background_center;
  tform = maketform('affine',obj.registrationData{n}.affine);
  [bkgdImage1,xdata,ydata] = imtransform(bkgdImage0,tform,'linear',...
    'UData',[1,obj.ncs(n)],'VData',[1,obj.nrs(n)],...
    'XData',obj.arena_center_mm(1)+obj.arena_radius_mm*[-1,1],...
    'YData',obj.arena_center_mm(2)+obj.arena_radius_mm*[-1,1],...
    'Size',[obj.ncs(n),obj.nrs(n)]);
  bkgdImage{i} = bkgdImage1;
end

  
if numel(ns) == 1,
  bkgdImage = bkgdImage{1};
end