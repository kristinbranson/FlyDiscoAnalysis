function CleanOutExpDirs_leavetracking(folder_path_from_experiment_index)
% delete pipeline files before rerunning pipeline
%list of pipeline generated files, add more if needed

todeletefiles = {'automatic_checks_complete_info.mat',...
    'automatic_checks_complete_results.txt',...
    'automatic_checks_incoming_info.mat',...
    'automatic_checks_incoming_results.txt',...
    'indicatordata.mat',...
    'indicator_log.txt',...
    'ledregistrationimage.png',...
    'perframefeatures_info.mat',...
    'registered_trx.mat',...
    'registrationdata.mat',...
    'registrationdata.txt',...
    'registrationimage.png',...
    'sexclassifier.mat',...
    'sexclassifier_diagnostics.txt',...
    'wingtracking_results.mat'};

todeletewildcard = {'ANALYSIS*',...
    'ctrax_results_movie_*.mp4','apt*','ctrax*'};

todeletedir = {'/perframe'};

for i = 1:numel(folder_path_from_experiment_index)
    expdir = folder_path_from_experiment_index{i};
    % delete files
    for j = 1:numel(todeletefiles)
        filestr = fullfile(expdir,todeletefiles{j});
        if exist(filestr,'file')
            delete(filestr)
            if exist(filestr,'file')
                fprintf('not deleted %s\n',filestr)
                return
            end
        end
    end
    %remove dirs
    for j = 1:numel(todeletedir)
        dirstr = fullfile(expdir,todeletedir{j});

        if exist(dirstr,'dir')
            filelist = dir(dirstr);
            for k = 3:numel(filelist)
                filestr = fullfile(expdir,todeletedir{j},filelist(k).name);

                if exist(filestr,'file')
                    delete(filestr)
                    if exist(filestr,'file')
                        fprintf('not deleted %s\n',filestr)
                        return
                    end
                end
            end

            rmdir(dirstr);
            if exist(dirstr,'file')
                fprintf('not deleted %s\n',dirstr)
                return
            end
        end
    end
    % delete files with variable names
    for j = 1:numel(todeletewildcard)
        filelist = dir(fullfile(expdir,todeletewildcard{j}));
        for k = 1:numel(filelist)
            filestr = fullfile(expdir,filelist(k).name);
            if exist(filestr,'file')
                delete(filestr)
                if exist(filestr,'file')
                    fprintf('not deleted %s\n',filestr)
                    return
                end
            end
        end
    end

end