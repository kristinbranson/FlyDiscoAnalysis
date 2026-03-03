function [mask_count, trimmed_mask_from_mask_index, trimmed_mask_lo_bound_in_frame_coords, trimmed_mask_hi_bound_in_frame_coords, trimmed_mask_shape] = ...
  compute_trimmed_masks_and_bounds(has_explicit_masks, mask_from_mask_index, ledIndicatorPointsXY, nr, nc, indicator_params)

% Compute some things we'll need for extracting the LED signals from each
% movie frame.
% If we're using masks, that means to compute the overlap of the mask and the
% movie frame.
% If not using masks, things are simpler.
% (Values in ledIndicatorPoints are guaranteed to be integral.)
if has_explicit_masks ,
  mask_width = size(mask_from_mask_index,2) ;
  mask_height = size(mask_from_mask_index,1) ;
  mask_count = size(mask_from_mask_index,3) ;
else
  ledwindowr = indicator_params.indicatorwindowr ;
  mask_height = 2*ledwindowr + 1 ;
  mask_width = 2*ledwindowr + 1 ;
  mask_count = 1 ;
end
mask_shape = [ mask_height mask_width ] ;
frame_shape = [ nr nc ] ;  % yx order
% There are three rectangles to deal with: the frame, the mask, and the
% trimmed mask.  The center of the mask in frame coordinates is given by
% ledIndicatorPoints, which are all integers.  The center of the mask in mask
% coordinates is the middle of the mask.  The trimmed mask is the part of the
% mask that overlaps with the frame.  In what follows, the bounds of an image
% in x/y, in a particular reference frame, are given by a "bound", a
% two-element row vector in yx order.  The bounds endpoints are inclusive.  An
% image has a lower bound and an upper bound.  We also talk about "anchors",
% which are two-element row vectors in yx order that determine the location of
% something.  The center of an image is one kind of anchor, but an anchor is
% not always a center.  They are often the center of that thing, but not
% always.  We also deal with "extents".  The rule is that anchor + extent =
% bound.  And image has a lower extent and an upper extent.  Bounds, anchors,
% and centers are all relative to a coordinate system, but extents are not.
% Lower extents are typically all-negative, upper extents are typically all-postive.
% All elements of anchors, bounds, and extents are integers.
% The lower extent of an image with odd dimensions will be the negative of the
% upper extent, when the image center is the anchor.
mask_center_in_mask_coords = [ round((mask_height+1)/2) round((mask_width+1)/2) ] ;  % yx order
mask_hi_extent = mask_shape - mask_center_in_mask_coords ;  % yx order
mask_lo_extent = [1 1] - mask_center_in_mask_coords ;  % yx order
mask_center_in_frame_coords = [ ledIndicatorPointsXY(2) ledIndicatorPointsXY(1) ] ;  % yx order
mask_lo_bound_in_frame_coords = mask_center_in_frame_coords + mask_lo_extent ;
mask_hi_bound_in_frame_coords = mask_center_in_frame_coords + mask_hi_extent ;
trimmed_mask_lo_bound_in_frame_coords = max([1 1], mask_lo_bound_in_frame_coords) ;
trimmed_mask_hi_bound_in_frame_coords = min(mask_hi_bound_in_frame_coords, frame_shape) ;
trimmed_mask_lo_extent = trimmed_mask_lo_bound_in_frame_coords - mask_center_in_frame_coords ;
trimmed_mask_hi_extent = trimmed_mask_hi_bound_in_frame_coords - mask_center_in_frame_coords ;
trimmed_mask_shape = trimmed_mask_hi_extent - trimmed_mask_lo_extent + [1 1] ;
trimmed_mask_anchor_in_frame_coords = mask_center_in_frame_coords ;
trimmed_mask_anchor_in_mask_coords = mask_center_in_mask_coords ;
trimmed_mask_lo_bound_in_mask_coords = trimmed_mask_anchor_in_mask_coords + trimmed_mask_lo_extent ;
trimmed_mask_hi_bound_in_mask_coords = trimmed_mask_anchor_in_mask_coords + trimmed_mask_hi_extent ;
trimmed_mask_lo_bound_in_frame_coords = trimmed_mask_anchor_in_frame_coords + trimmed_mask_lo_extent ;
trimmed_mask_hi_bound_in_frame_coords = trimmed_mask_anchor_in_frame_coords + trimmed_mask_hi_extent ;
if has_explicit_masks ,
  trimmed_mask_from_mask_index = ...
    mask_from_mask_index(trimmed_mask_lo_bound_in_mask_coords(1):trimmed_mask_hi_bound_in_mask_coords(1), ...
                         trimmed_mask_lo_bound_in_mask_coords(2):trimmed_mask_hi_bound_in_mask_coords(2), ...
                         :) ;
else
  trimmed_mask_from_mask_index = true(trimmed_mask_shape) ;
end

end  % function
