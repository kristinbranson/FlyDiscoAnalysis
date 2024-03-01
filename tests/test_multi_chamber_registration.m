rng(0);  % Want 'frozen' random numbers

% Make sure get same numbers for a single chamber
sample_count = 100 ;
x = normrnd(0,1,[sample_count 1]) ;
y = normrnd(0,1,[sample_count 1]) ;
offX = normrnd(0,1) ;
offY = normrnd(0,1) ;
scale = exprnd(1) ;
offTheta = unifrnd(-pi,+pi) ;
[x_reg_single,y_reg_single] = registerForSingleChamber(x,y,offX,offY,offTheta,scale) ;
[x_reg_multi,y_reg_multi] = registerForMultipleChambers(x,y,offX,offY,offTheta,scale) ;
A_single = affineTransformMatrixFromOffsetsAndScale(offX,offY,offTheta,scale) ;
A_multi = affineTransformMatricesFromOffsetsAndScale(offX,offY,offTheta,scale) ;
assert(all(abs(x_reg_single-x_reg_multi)<eps)) ;
assert(all(abs(y_reg_single-y_reg_multi)<eps)) ;
assert(all(all(abs(A_single-A_multi)<eps))) ;

% Make sure it gets the chamber assignments right
chamber_sample_count = 10 ;  % samples per chamber
chamber_count = 9 ;
sample_count = chamber_count * chamber_sample_count ;
r_center_from_chamber_index = 3*[ 1 1 ; 1 2 ; 1 3 ; 2 1 ; 2 2 ; 2 3 ; 3 1 ; 3 2 ; 3 3 ]' ;  % 2 x chamber_count
r_center_from_chamber_index = r_center_from_chamber_index(:,1:chamber_count) ;
radius_full = unifrnd(0,1,[1 chamber_count chamber_sample_count]) ;
theta_full = unifrnd(-pi,+pi,[1 chamber_count chamber_sample_count]) ;
dx_full = radius_full .* cos(theta_full) ;
dy_full = radius_full .* sin(theta_full) ;
dr_full = vertcat(dx_full,dy_full) ;  % 2 x chamber_count x chamber_sample_count
r_full = r_center_from_chamber_index + dr_full ;  % broadcast, 2 x chamber_count x chamber_sample_count
chamber_index_from_chamber_index_and_chamber_sample_index = repmat(1:chamber_count, [1 1 chamber_sample_count]) ;
r = reshape(r_full, [2 sample_count]) ;
chamber_index_from_sample_index = reshape(chamber_index_from_chamber_index_and_chamber_sample_index, [1 sample_count]) ;
% Shuffle the points
sample_index_from_shuffled_sample_index = randperm(sample_count) ;
r_shuffled = r(:,sample_index_from_shuffled_sample_index) ;
chamber_index_from_shuffled_sample_index = chamber_index_from_sample_index(sample_index_from_shuffled_sample_index) ;
scale = exprnd(1) ;
offTheta = unifrnd(-pi,+pi) ;
x_shuffled = r_shuffled(1,:) ;
y_shuffled = r_shuffled(2,:) ;
offX = -r_center_from_chamber_index(1,:) ;
offY = -r_center_from_chamber_index(2,:) ;
[x_reg_multi,y_reg_multi,output_chamber_index_from_shuffled_sample_index] = registerForMultipleChambers(x_shuffled,y_shuffled,offX,offY,offTheta,scale) ;
assert(isequaln(output_chamber_index_from_shuffled_sample_index,chamber_index_from_shuffled_sample_index)) ;

% Do single-chamber registration for each chamber's samples, compare with
% multi-chamber registration results.
x_reg_single = zeros(1,sample_count) ;
y_reg_single = zeros(1,sample_count) ;
for chamber_index = 1 : chamber_count ,
  r_center = r_center_from_chamber_index(:,chamber_index) ;
  offX_single = -r_center(1) ;
  offY_single = -r_center(2) ;
  is_in_chamber_from_shuffled_sample_index = (chamber_index_from_shuffled_sample_index == chamber_index) ;
  x_shuffled_chamber = x_shuffled(is_in_chamber_from_shuffled_sample_index) ;
  y_shuffled_chamber = y_shuffled(is_in_chamber_from_shuffled_sample_index) ;
  [x_shuffled_chamber_reg,y_shuffled_chamber_reg] = registerForSingleChamber(x_shuffled_chamber,y_shuffled_chamber,offX_single,offY_single,offTheta,scale) ;
  x_reg_single(is_in_chamber_from_shuffled_sample_index) = x_shuffled_chamber_reg ;
  y_reg_single(is_in_chamber_from_shuffled_sample_index) = y_shuffled_chamber_reg ;  
end
assert(all(abs(x_reg_single-x_reg_multi)<eps)) ;
assert(all(abs(y_reg_single-y_reg_multi)<eps)) ;

% Assemble all the single-chamber transforms, compare with multi
A_multi = affineTransformMatricesFromOffsetsAndScale(offX,offY,offTheta,scale) ;
A_single = zeros(3,3,chamber_count) ;
for chamber_index = 1 : chamber_count ,
  r_center = r_center_from_chamber_index(:,chamber_index) ;
  offX_single = -r_center(1) ;
  offY_single = -r_center(2) ;
  A_single(:,:,chamber_index) = affineTransformMatrixFromOffsetsAndScale(offX_single,offY_single,offTheta,scale) ;
end
assert(all(all(all(abs(A_single-A_multi)<eps)))) ;

% Declare success
fprintf('All tests passed.\n') ;

