function result = are_registrationdata_basically_same(test_data, reference_data)
% Compares the two inputs to see if they are the same in all the important
% ways.

% Don't compare the bkgdImages.  Those are computed by FlyTracker, and
% FlyTracker does not seem to be 100% deterministic.
test_data_1 = rmfield(test_data, 'bkgdImage') ;
reference_data_1 = rmfield(reference_data, 'bkgdImage') ;
result = isequaln(test_data_1, reference_data_1) ;
end
