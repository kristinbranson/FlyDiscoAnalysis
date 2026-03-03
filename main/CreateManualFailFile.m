function [success] = CreateManualFailFile(expdir,failure_reason,varargin)
%creates a manual fail file in experiment directory
% if the fail file doesn't exist, write new file
% if the fail already exits, defaults to appending, but if Replace = true will
% rewrite with new value.
% expdir can be a cell array of expdirs with the same error. 
% % current failure categories:
% 'wrong eye color in cross';
% 'crud in bubble';
% 'testing';
% 'bad eye color effector';
% 'line name error';
% 'long load time';
% 'cold plate too cold';
% 'bubble placement error';
% 'tracked dead or damaged fly';
% 'starvation too long';
% 'bad video';
% 'bubble occluded';
% 'tracked fly outside bubble';

[file_name,Replace] = myparse(varargin,'file_name','manual_fail.txt','Replace',false);
success = false(1,numel(expdir));
if ~iscell(expdir)
    cnt = 1;
else
    cnt = numel(expdir);
end


for i = 1:cnt
    if ~iscell(expdir)
        failfile = fullfile(expdir,file_name);
    else
    failfile = fullfile(expdir{i},file_name);
    end

    if ~exist(failfile,'file')
        % if failfile doesn't exist - write file
        fid = fopen(failfile,'w');
        if fid < 0,
            sprintf('Cound not write to %s',failfile')
            success(i) = false;
            continue;
        end
        fprintf(fid,'%s',failure_reason);
        fclose(fid);
        if exist(failfile,'file')
            success(i) = true;
        end
    else
        if Replace
            % if file exists and replace = true, overwrite failfile
            fid = fopen(failfile,'w+');
            if fid < 0,
                sprintf('Cound not write to %s',failfile')
                success(i) = false;
                continue;
            end
            fprintf(fid,'%s',failure_reason);
            fclose(fid);
            success(i) = true;
        else
            % if file exists, appends another failure_reason to fail file
            fid = fopen(failfile,'a');
            if fid < 0,
                sprintf('Cound not write to %s',failfile')
                success(i) = false;
                continue;
            end
            fprintf(fid,'\n%s',failure_reason);
            fclose(fid);
            success(i) = true;
        end
    end

end