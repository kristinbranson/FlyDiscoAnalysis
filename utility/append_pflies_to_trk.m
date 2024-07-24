function result = append_pflies_to_trk(input_trk, trx)
% Appends the pflies in the Ctrax-style trx struct to the tracks in input_trk (a FlyTracker-style trk
% struct).

% Examples of the inputs typically look like:
%
% input_trk = 
% 
%   struct with fields:
% 
%                names: {1×35 cell}
%                 data: [9×109837×35 double]
%     flies_in_chamber: {[1]  [2]  [3]  [4]  [5]  [6]  [7]  [8]  [9]}
%                flags: [0×6 double]
%
%
% input_trk.names
% 
% ans =
% 
%   1×35 cell array
% 
%   Columns 1 through 15
% 
%     {'pos x'}    {'pos y'}    {'ori'}    {'major axis len'}    {'minor axis len'}    {'body area'}    {'fg area'}    {'img contrast'}    {'min fg dist'}    {'wing l x'}    {'wing l y'}    {'wing r x'}    {'wing r y'}    {'wing l ang'}    {'wing l len'}
% 
%   Columns 16 through 32
% 
%     {'wing r ang'}    {'wing r len'}    {'leg 1 x'}    {'leg 1 y'}    {'leg 2 x'}    {'leg 2 y'}    {'leg 3 x'}    {'leg 3 y'}    {'leg 4 x'}    {'leg 4 y'}    {'leg 5 x'}    {'leg 5 y'}    {'leg 6 x'}    {'leg 6 y'}    {'leg 1 ang'}    {'leg 2 ang'}    {'leg 3 ang'}
% 
%   Columns 33 through 35
% 
%     {'leg 4 ang'}    {'leg 5 ang'}    {'leg 6 ang'}
%
% trx = 
% 
%   1×18 struct array with fields:
% 
%     moviename
%     moviefile
%     firstframe
%     off
%     endframe
%     nframes
%     fps
%     pxpermm
%     id
%     sex
%     timestamps
%     dt
%     x
%     y
%     theta
%     a
%     b
%     xwingl
%     ywingl
%     xwingr
%     ywingr
%     x_mm
%     y_mm
%     a_mm
%     b_mm
%     theta_mm
%     is_pfly

% Break out the feature names
feature_names = input_trk.names ;
feature_count = numel(feature_names) ;

% Break out the input data array
input_data = input_trk.data ;  % fly_count x frame_count x feature_count
real_fly_count = size(input_data, 1) ;
frame_count = size(input_data, 2) ;
feature_count_check = size(input_data, 3) ;
if feature_count ~= feature_count_check ,
  error(strcat('Dimensions of input_trk.names and input_trk.data are not compatible.  ', ...
               'numel(input_trk.names) and size(input_trk.data,3) must be equal.  ', ...
               'First is %d, second is %d.'), ...
        feature_count, ...
        feature_count_check) ;
end

% Build a struct that maps from feature names to indices into third dimension
% of input_data
ind = generate_feature_index(feature_names) ;

% Count how many pflies are present
is_pfly = [ trx.is_pfly ] ;
pfly_count = sum(is_pfly) ;

% Create pfly_data, like input_data except for the pflies
pfly_data = nan(pfly_count, frame_count, feature_count) ;
for pfly_index = 1 : pfly_count ,
  trxi = trx(pfly_index) ;
  firstframe = trxi.firstframe ;
  endframe = trxi.endframe ;
  pfly_data(pfly_index, firstframe:endframe, ind.pos_x) = trxi.x ;
  pfly_data(pfly_index, firstframe:endframe, ind.pos_y) = trxi.y ;
  pfly_data(pfly_index, firstframe:endframe, ind.ori)   = -trxi.theta ;  % different conventions
  pfly_data(pfly_index, firstframe:endframe, ind.major_axis_len) = 4*trxi.a ;
  pfly_data(pfly_index, firstframe:endframe, ind.minor_axis_len) = 4*trxi.b ;
end

% Check that the number of frames has not changed
if size(pfly_data,2) ~= frame_count ,
  error('After adding pfly tracks to input_trk.data, the frame count has changed.  This means that at least one pfly had a bad firstframe or endframe value.') ;
end

% input_trk.flies_in_chamber is a cell array with arena_count elements.  Each
% element is a row vector listing the fly indices in that arena.  So a (too
% long) name might be fly_index_from_arena_fly_index_from_arena_index.
% 
% Example value: {[1]  [2]  [3]  [4]  [5]  [6]  [7]  [8]  [9]}
%
% So we want to add one pfly to each per-arena list.  
%
% The flytracker calibration lists the arena centers in a particular order.
% The pflies are guaranteed to be in the same order in trx, so it's safe to
% assume that the pfly_index equals the arena_index.

% Break out input flies_in_chamber, do a sanity check
input_flies_in_chamber = input_trk.flies_in_chamber ;
arena_count = numel(input_flies_in_chamber) ;
if arena_count ~= pfly_count ,
  error('The number of arenas (%d) does not match the number of pflies (%d)', arena_count, pfly_count) ;
end

% Add the pfly to each chamber
flies_in_chamber = cell(size(input_flies_in_chamber)) ;
for arena_index = 1 : arena_count ,
  input_flies_in_this_chamber = input_flies_in_chamber{arena_index} ;
  pfly_index = arena_index ;  % Guaranteed by how pflies are generated
  fly_index = real_fly_count + pfly_index ;  % pflies come after real flies
  flies_in_chamber{arena_index} = horzcat(input_flies_in_this_chamber, fly_index) ;
end

% Assemble the result trk struct
result = struct() ;
result.names = feature_names ;
result.data = cat(1, input_data, pfly_data) ;
result.flies_in_chamber = flies_in_chamber ;
result.flags = input_trk.flags ;  % no idea what this is...

end
