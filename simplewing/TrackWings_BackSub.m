function [iswing,isfore_thresh,idxfore_thresh,npxfore_thresh,fore2body,fore2dbkgd,debugdata] = ...
  TrackWings_BackSub(im,bgmodel,isarena,params,debugdata)

dbkgd = bgmodel - im;
dbkgd(~isarena) = 0;

%% morphology to get foreground, body wing pixels

isbody_thresh = dbkgd >= params.mindbody;
isbody = imdilate(isbody_thresh,params.se_dilate_body);

iswing_high = ~isbody & dbkgd >= params.mindwing_high;
iswing_low = ~isbody & dbkgd >= params.mindwing_low;
iswing = imreconstruct(iswing_high,iswing_low,4);
iswing = imopen(iswing,params.se_open_wing);
iswing = imclose(iswing,params.se_open_wing);

isfore_thresh = imclose(isbody_thresh | iswing,params.se_dilate_body);

%isfore_thresh = dbkgd > params.mindwing_high;
idxfore_thresh = find(isfore_thresh);
npxfore_thresh = nnz(isfore_thresh);

fore2dbkgd = dbkgd(isfore_thresh);
fore2body = fore2dbkgd >= params.mindbody;

if debugdata.DEBUG > 1,
  DebugPlot_BackSub();
end

  function DebugPlot_BackSub()
    
    if debugdata.DEBUG > 1,
      axcurr = 1;
      if isnan(debugdata.hims(axcurr)),
        cla(debugdata.hax(axcurr));
      end
      [nr,nc,~] = size(debugdata.im);
      imtmp = repmat(debugdata.im(:),[1,3]);
      imtmp(isbody,1) = min(imtmp(isbody,1)+100,255);
      imtmp(iswing,2) = min(imtmp(iswing,2)+100,255);
      imtmp = uint8(reshape(imtmp,[nr,nc,3]));
      if isnan(debugdata.hims(1)),
        debugdata.hims(axcurr) = image(imtmp,'Parent',debugdata.hax(axcurr)); axis(debugdata.hax(axcurr),'image','off');
      else
        set(debugdata.hims(axcurr),'CData',imtmp);
      end
      drawnow;
      title(debugdata.hax(axcurr),sprintf('Segmentation of frame %d into bg, body and wing',debugdata.t));
    end
    
  end

end