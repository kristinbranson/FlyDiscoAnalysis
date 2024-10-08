% closest fly, based on dell2ell
% AR 3/8/2018
function [data,units,mind] = compute_closestfly_ell2ell(trx,n,dosave_d)

if nargin < 3
    dosave_d = true;
end

flies = trx.exp2flies{n};
nflies = numel(flies);
closestfly = cell(1,nflies);
mind = cell(1,nflies);

% parfor i1 = 1:nflies,
for i1 = 1:nflies,
    fly1 = flies(i1);
    fprintf('fly1 = %d\n',fly1);
    flies2 = flies(trx(fly1).chamber_index==trx.chamber_index(flies));
    d = nan(numel(flies2),trx(fly1).nframes);
    
    % use dcenter2nose and major axis length to computer upper and lower
    % bounds on dell2ell
    
    [mindupper,dlower] = dell2ell_bounds(trx,fly1,flies2);
    
    for i2 = 1:numel(flies2)
        fly2 = flies2(i2);
        if fly1 == fly2
            continue
        end
        %only try for frames where lower bound is smaller than min upper
        %bound
        idx1try = find(mindupper >= dlower(i2,:));
        d(i2,:) = dell2ell_pair(trx,fly1,fly2,idx1try);
    end
    [mind{i1},closesti] = min(d,[],1);
    closestfly{i1} = flies2(closesti);
    closestfly{i1}(isnan(mind{i1})|isinf(mind{i1})) = nan;
end

% so that we don't compute dcenter twice
if dosave_d,
    data = mind; %#ok<NASGU>
    units = parseunits('mm'); %#ok<NASGU>
    filename = trx.GetPerFrameFile('dell2ell',n);
    try
        save(filename,'data','units')
    catch ME,
        warning('Could not save file %s:\%s',filename,getReport(ME));
    end
    
end


data = closestfly;
units = parseunits('unit');
