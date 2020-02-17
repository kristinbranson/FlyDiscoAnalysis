function [wingtrx,debugdata,istouching] = TrackWingsOneFrame(im,bgmodel,isarena,trxcurr,wingtrxprev,params,XGRID,YGRID,debugdata)


%% morphology to get foreground, body wing pixels
[iswing,isfore_thresh,idxfore_thresh,npxfore_thresh,fore2body,fore2dbkgd,debugdata] = ...
  TrackWings_BackSub(im,bgmodel,isarena,params,debugdata);

[fore2fly,xgrid_isfore,ygrid_isfore,debugdata,istouching] = TrackWings_SegmentFlies(im,isfore_thresh,idxfore_thresh,npxfore_thresh,fore2body,...
  trxcurr,params,XGRID,YGRID,debugdata);

[wingtrx,debugdata] = TrackWings_FitWings(fore2fly,xgrid_isfore,ygrid_isfore,isfore_thresh,iswing,fore2dbkgd,idxfore_thresh,trxcurr,wingtrxprev,params,debugdata);


