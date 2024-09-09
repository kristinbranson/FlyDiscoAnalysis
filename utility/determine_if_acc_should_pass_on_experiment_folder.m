function result = determine_if_acc_should_pass_on_experiment_folder(expdir)
% For testing, it's nice to have a way to indicate that we *want* ACC checks
% to fail for a particular experiment.  In the test suite, we indicate this by
% having a file named ACC_SHOULD_FAIL in the experiment folder.  If
% this file is missing (the more common case), we expect the ACC checks to
% pass.

file_name = fullfile(expdir, 'ACC_SHOULD_FAIL') ;
result = ~logical(exist(file_name, 'file')) ;
