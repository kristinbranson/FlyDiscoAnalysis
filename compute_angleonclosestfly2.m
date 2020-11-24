% AR 20200825
%  RecomputeAngleonClosestFly
% SET npoints in this file and its passed along to ellipsedist_hack
% calls
% compute_closestfly_nose2ell2
% calls
% dnose2ell_pair2 (i set nsamples to 150, was 20)
% calls
% ellipsedist_hack2 (nsample = 100, if not passed in) - in jaaba code
% calls
% ellipsepoints_even (npoints = 60, if not passed in) - in jaaba code

function [data,units] = compute_angleonclosestfly2(trx,n)

npoints = 150;
% don't save in code (see below)
dosave_d = false;
[~,~,~,angle] = compute_closestfly_nose2ell2(trx,n,dosave_d,npoints);


% now done in ellipsepoints_even
% change angleonclosestfly to be in new range | 0 : pi |
% angle2 = cell(size(angle));
% for i = 1:numel(angle)
%     angle2{i} = pi - abs(angle{i});
% end


% saves angle as angleonclosestfly2 
data = angle; %#ok<NASGU>
units = parseunits('rad'); %#ok<NASGU>
filename = trx.GetPerFrameFile('angleonclosestfly',n);
[path,filename,extension] = fileparts(filename);
filename = fullfile(path,[filename,'2',extension])
if exist(filename,'file'),
    try
        delete(filename);
    catch ME
        warning('Could not delete file %s: %s',filename,getReport(ME));
    end
end
try
    save(filename,'data','units');
catch ME
    warning('Could not save to file %s: %s',filename,getReport(ME));
end






