function meanim = AverageAnatomyImage(line_names,anndata,varargin)

[useallims,maxqi,filsizexy,filsizez,savedir,cm,hfig] = ...
  myparse(varargin,'useallims',false,'maxqi',.1,...
  'filsizexy',5,'filsizez',2,'savedir','',...
  'colormap',kjetsmooth(256),'hfig',[]);

if isempty(hfig),
  hfig = figure;
else
  if ishandle(hfig),
    set(0,'CurrentFigure',hfig);
    clf;
  else
    figure(hfig);
  end
end
hax = createsubplots(2,1,.05);
colormap(cm);

if ~isempty(savedir) && ~exist(savedir,'dir'),
  mkdir(savedir);
end

nresamplexy = round(filsizexy*2);
nresamplez = round(filsizez*2);

% [xgrid,ygrid,zgrid] = ndgrid(floor(-2.5*sigmaxy):ceil(2.5*sigmaxy),...
%   floor(-2.5*sigmaxy):ceil(2.5*sigmaxy),...
%   floor(-2.5*sigmaz):ceil(2.5*sigmaz));
% fil = single(exp(-.5*(xgrid.^2/sigmaxy^2+ygrid.^2/sigmaxy^2+zgrid.^2/sigmaz^2)));
fil = ones([filsizexy,filsizez,2],'single');
fil = fil / sum(fil(:));


[~,lineidx] = ismember(line_names,anndata.line_names);

if any(lineidx==0),
  fprintf('Lines %s are not imaged\n',sprintf('%s ',line_names{lineidx==0}));
  line_names = line_names(lineidx~=0);
  lineidx = lineidx(lineidx~=0);
end
nlinescurr = numel(lineidx);
imidx_perline = cell(1,nlinescurr);

isallowed = [anndata.imdata.qi] <= maxqi & ...
  strcmpi({anndata.imdata.gender},'Female');


if useallims,
    
  for ii = 1:nlinescurr,
    i = lineidx(ii);
    idxcurr = find(isallowed & strcmp({anndata.imdata.line},anndata.line_names{i}));
    if isempty(idxcurr),
      idxcurr = anndata.imagechoice(i);
      if isnan(idxcurr) || isempty(idxcurr),
        fprintf('No images found for line %s\n',line_names{i});
        continue;
      end
    end
    imidx_perline{ii} = idxcurr;    
  end
  
else
  
  for ii = 1:nlinescurr,
    i = lineidx(ii);
    idxcurr = anndata.imagechoice(i);
    if isnan(idxcurr) || isempty(idxcurr),
      fprintf('No images found for line %s\n',line_names{i});
      continue;
    end
    imidx_perline{ii} = idxcurr;
  end  
end

nlinescount = 0;
for ii = 1:nlinescurr,
  nimscurr = 0;
  for jj = 1:numel(imidx_perline{ii}),
    j = imidx_perline{ii}(jj);
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
      meanimcurr = imnorm;
    else
      meanimcurr = meanimcurr + imnorm;
    end

    imagesc(max(imnorm,[],3)','Parent',hax(1));
    title(hax(1),anndata.imdata(j).name,'Interpreter','none');
    imagesc(max(meanimcurr,[],3)'/(nimscurr+1),'Parent',hax(2));
    title(hax(2),line_names{ii},'Interpreter','none');
    axis(hax,'image','off');
    set(hax,'CLim',[0,1]);
    drawnow;
    
    if ~isempty(savedir) && useallims,
      [~,name] = fileparts(anndata.imdata(j).name);
      outfile = fullfile(savedir,sprintf('imnorm_%s.png',name));
      Irgb = im2uint8(colormap_image(max(imnorm,[],3)',cm));
      imwrite(Irgb,outfile,'png');
    end

    nimscurr = nimscurr+1;
  end
  meanimcurr = meanimcurr / nimscurr;
  meanimcurr = meanimcurr / max(meanimcurr(:));
  
  if ~isempty(savedir),
    outfile = fullfile(savedir,sprintf('meanim_%s.png',line_names{ii}));
    Irgb = im2uint8(colormap_image(max(meanimcurr,[],3)',cm));
    imwrite(Irgb,outfile,'png');
  end
    
  if nlinescount == 0,
    meanim = meanimcurr;
  else
    meanim = meanim + meanimcurr;
  end

  nlinescount = nlinescount + 1;

  if ~isempty(savedir),
    outfile = fullfile(savedir,sprintf('meanim%02d.png',nlinescount));
    Irgb = im2uint8(colormap_image(max(meanim/nlinescount,[],3)',cm));
    imwrite(Irgb,outfile,'png');
  end

  
end

meanim = meanim / nlinescount;

if ~isempty(savedir),
  outfile = fullfile(savedir,'meanim.png');
  Irgb = im2uint8(colormap_image(max(meanim,[],3)',cm));
  imwrite(Irgb,outfile,'png');
end
