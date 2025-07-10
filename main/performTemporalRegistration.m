function [trx, timestamps, registration_data, newid2oldid] = ...
        performTemporalRegistration(dotemporalreg, trx, timestamps, registration_data, newid2oldid, fns_notperframe, expdir, dataloc_params)

% If not doing temporal registration, exit early
if ~dotemporalreg ,
    fprintf('NOT applying temporal registration.\n');
    return
end

% As of 2024-02-16, none of the analysis-protocol files in
% FlyDiscoAnalysis/settings-internal use temporal registration, and apparently
% it hasn't been used for several years.  So just note that this code has not
% been under any selective pressure for a while.

% how long did we actually record for?
timestamps = timestamps - timestamps(1);

% how long did we want to record for? read from config file
configfile = dir(fullfile(expdir,dataloc_params.configfilepattern));
if isempty(configfile),
    error('Could not find config file for experiment %s',expdir);
end
configfile = fullfile(expdir,configfile(1).name);
config_params = ReadParams(configfile);

% read flies loaded time -- could use metadata tree loader, but this is
% pretty simple
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
if ~exist(metadatafile,'file'),
    error('Could not find metadata file %s',metadatafile);
end
fid = fopen(metadatafile,'r');
if fid < 0,
    error('Could not open metadata file %s for reading',metadatafile);
end
cleaner1 = onCleanup(@()(fclose(fid))) ;
while true,
    s = fgetl(fid);
    if ~ischar(s),
        fclose(fid);
        error('Could not find "seconds_fliesloaded" in metadata file %s',metadatafile);
    end
    m = regexp(s,'seconds_fliesloaded\s*=\s*"([^"]*)"','tokens','once');
    if ~isempty(m),
        seconds_fliesloaded = str2double(m{1});
        break;
    end
end
fclose(fid);

% how much time should we crop from the beginning?
timeCropStart = registration_params.maxFliesLoadedTime - seconds_fliesloaded;
if timeCropStart < 0,
    warning('Load time = %f seconds, greater than max allowed load time = %f seconds.',...
        seconds_fliesloaded,registration_params.maxFliesLoadedTime);
    timeCropStart = 0;
end

% find closest timestamp to timeCropStart
i0 = find(timestamps >= timeCropStart,1);
if isempty(i0),
    error('No timestamps occur after timeCropStart = %f. Cannot crop start.',timeCropStart);
end
if i0 > 1 && (timeCropStart-timestamps(i0-1)) < (timestamps(i0)-timeCropStart),
    i0 = i0 - 1;
end

% how long is the video currently?
recordLengthCurr = timestamps(end)-timestamps(i0);

% how long should the video be?
recordLengthIdeal = config_params.RecordTime - registration_params.extraBufferFliesLoadedTime - ...
    (registration_params.maxFliesLoadedTime - registration_params.minFliesLoadedTime);

% how much time should we crop from the end?
if recordLengthCurr < recordLengthIdeal,
    warning('Cropped video is %f seconds long, shorter than ideal length %f seconds.',recordLengthCurr,recordLengthIdeal);
    i1 = numel(timestamps);
else
    i1 = find(timestamps - timestamps(i0) >= recordLengthIdeal,1);
    if isempty(i1),
        warning('No timestamps occur after timeCropEnd = %f. Cannot crop end.',timestamps(i0)+recordLengthIdeal);
        i1 = numel(timestamps);
    else
        if i1 > 1 && ...
                (recordLengthIdeal - (timestamps(i1-1)-timestamps(i0))) < ...
                ((timestamps(i1)-timestamps(i0)) - recordLengthIdeal),
            i1 = i1 - 1;
        end
    end
end
registration_data.seconds_crop_start = timestamps(i0);
registration_data.start_frame = i0;
registration_data.seconds_crop_end = timestamps(end)-timestamps(i1);
registration_data.end_frame = i1;

fns = setdiff(fieldnames(trx),fns_notperframe);
isperframe = true(1,numel(fns));
nperfn = nan(1,numel(fns));
ndelete = 0;
if ~isempty(trx),
    for j = 1:numel(fns),
        fn = fns{j};
        for i = 1:numel(trx),
            if ~isnumeric(trx(i).(fn)),
                isperframe(j) = false;
                break;
            end
            if i == 1,
                ncurr = trx(i).nframes - numel(trx(i).(fn));
            else
                if trx(i).nframes - numel(trx(i).(fn)) ~= ncurr,
                    isperframe(j) = false;
                    break;
                end
            end
        end
        if all([trx.nframes]) == trx(1).nframes && ...
                trx(1).nframes > 1 && numel(trx(1).(fn)) == 1,
            isperframe(j) = false;
        end
        if isperframe(j),
            nperfn(j) = ncurr;
        end
    end

    nperfn = nperfn(isperframe);
    fns = fns(isperframe);
    ncropright = ceil(nperfn/2);
    ncropleft = nperfn - ncropright;

    trxdelete = false(1,numel(trx));
    for i = 1:numel(trx),
        if trx(i).firstframe > i1,
            trxdelete(i) = true;
            continue;
        end

        if trx(i).endframe < i0,
            trxdelete(i) = true;
            continue;
        end

        trx(i).nframes = min(i1,trx(i).endframe)-max(i0,trx(i).firstframe)+1;

        if trx(i).firstframe < i0,
            off = i0 - trx(i).firstframe;
            for j = 1:numel(fns),
                fn = fns{j};
                trx(i).(fn) = trx(i).(fn)(off+1+ncropleft(j):end);
            end
            trx(i).firstframe = i0;
        end

        if trx(i).endframe > i1,
            for j = 1:numel(fns),
                fn = fns{j};
                trx(i).(fn) = trx(i).(fn)(1:trx(i).nframes-nperfn(j));
            end
            trx(i).endframe = i1;
        end

        trx(i).off = -trx(i).firstframe + 1;
    end
    trx(trxdelete) = [];
    newid2oldid(trxdelete) = [];
    ndelete = nnz(trxdelete);
end
fprintf('Applied temporal registration, deleted %d trajectories.\n',ndelete);
