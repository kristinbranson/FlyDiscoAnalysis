function meanim = AverageLineAnatomy(line_name,matfilename)

% [useallims,maxqi,filsizexy,filsizez,savedir,cm,hfig] = ...
%   myparse(varargin,'useallims',false,'maxqi',.1,...
%   'filsizexy',5,'filsizez',2,'savedir','',...
%   'colormap',kjetsmooth(256),'hfig',[]);

meanim = [];
load(matfilename);

nresamplexy = round(filsizexy*2);
nresamplez = round(filsizez*2);

% [xgrid,ygrid,zgrid] = ndgrid(floor(-2.5*sigmaxy):ceil(2.5*sigmaxy),...
%   floor(-2.5*sigmaxy):ceil(2.5*sigmaxy),...
%   floor(-2.5*sigmaz):ceil(2.5*sigmaz));
% fil = single(exp(-.5*(xgrid.^2/sigmaxy^2+ygrid.^2/sigmaxy^2+zgrid.^2/sigmaz^2)));
fil = ones([filsizexy,filsizez,2],'single');
fil = fil / sum(fil(:));

isallowed = [anndata.imdata.qi] <= maxqi & ...
  strcmpi({anndata.imdata.gender},'Female');

i = find(strcmp(anndata.line_names,line_name));
if isempty(i),
  fprintf('No images found for line %s\n',line_name);
  return;
end

idxcurr = find(isallowed & strcmp({anndata.imdata.line},line_name));
if isempty(idxcurr),
  idxcurr = anndata.imagechoice(i);
  if isnan(idxcurr) || isempty(idxcurr),
    fprintf('No images found for line %s\n',line_name);
    return;
  end
end

nimscurr = 0;
for jj = 1:numel(idxcurr),
  j = idxcurr(jj);
  im = single(loadRaw2Stack(anndata.imdata(j).raw_file_system_path));
    
  if useallims,
    imsmall = im(1:nresamplexy:end,1:nresamplexy:end,1:nresamplez:end,:);
    sz = size(imsmall);
    imsmall = reshape(imsmall,[prod(sz(1:3)),sz(4)]);
    n255 = sum(abs(imsmall-255)<=3,1);
    [~,maski] = max(n255);
    if maski ~= 3,
      fprintf('Image %s has a bad mask channel\n',anndata.imdata(j).name);
      continue;
    end
    
    % make sure the reference mask is channel 2
    sumin = sum(imsmall(imsmall(:,3)<123,1:2),1);
    sumout = sum(im(im(:,3)>=123,1:2),1);
    refsum = sumin-sumout;
    [~,refi] = max(refsum);
    if refi ~= 2,
      fprintf('Image %s has a bad reference channel\n',anndata.imdata(j).name);
      continue;
    end
  end

  imgreen = im(:,:,:,1)/2^12;
  imgreen(anndata.mask==0) = 0;
  imgreen = convn(imgreen,fil,'same');
  imnorm = min(1,max(0,(imgreen-anndata.imdata(j).bg_thresh_fit)/...
    (anndata.imdata(j).fg_thresh_fit-anndata.imdata(j).bg_thresh_fit)));
  
  if nimscurr == 0,
    meanim = imnorm;
  else
    meanim = meanim + imnorm;
  end

  nimscurr = nimscurr+1;
end

if nimscurr == 0,
  return;
end

meanim = meanim / nimscurr;
  
outfile = fullfile(savedir,sprintf('meanim_%s.png',line_name));
Irgb = im2uint8(colormap_image(max(meanim/max(meanim(:)),[],3)',cm));
imwrite(Irgb,outfile,'png');
    
outfile = fullfile(savedir,sprintf('meanim_%s.mat',line_name));
save(outfile,'meanim');