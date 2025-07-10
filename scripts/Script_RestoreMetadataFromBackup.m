%% set up path
modpath

%% find dirs with bad metadata 
% rootdatadir = '/groups/branson/bransonlab/flydisco_data';
% tmp = dir(rootdatadir);
% metadatafile = 'Metadata.xml';
% fid1 = fopen('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/nometadataexplist_bransonlab.txt','w');
% fid2= fopen('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/malformedmetadataexplist_bransonlab.txt','w');
% 
% 
% for i = 1:numel(tmp)
%     if ismember(tmp(i).name,{'.','..'})
%         continue;
%     end
%     expdir = tmp(i).name;
%     metadatafilestr = fullfile(rootdatadir,expdir,metadatafile);
%     if exist(metadatafilestr,'file')
%         try
%             tmpM = ReadMetadataFile(metadatafilestr);
%         catch
%             fprintf(fid2,'%s\n',expdir);
%         end
%     else
%         fprintf(fid1,'%\n',expdir);
%     end
% end

%% run restore metadata from back up code on these experiemnts
% 
% explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/malformedmetadataexplist_bransonlab.txt','%s');
% % explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210804_testingMetadatafix/VNC_JHS_K_85321_RigD_20210322T164809',...
% %     '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210804_testingMetadatafix/VNC_JRC_SS25466_RigC_20210322T143411'};
% 
% for i = 2:numel(explist)
%     expdir = explist{i};
%     expdir = fullfile(rootdatadir,expdir)
%     metadatafile = fullfile(expdir,'Metadata.xml');
%     tmp2 = dir(fullfile(expdir,'Metadata.xml_*'));
%     if isempty(tmp2)
%         fprintf('no backup metadata file %s\n',expdir);
%     else
%         if numel(tmp2) > 1
%             dn = [];
%             %find most recent
%             for j = 1:numel(tmp2)
%                 dn(j) = datenum(tmp2(j).date);
%             end
%             [~,irecent] = max(dn);
%             metadatabackupfile = tmp2(irecent).name;
%         else
%             metadatabackupfile = tmp2(1).name;
%         end
%         
%         metadatabackupfile = fullfile(expdir,metadatabackupfile);
%         
%         [success] = RestoreMetadataFromBackup(metadatafile,metadatabackupfile);
%         if ~success
%             fprintf('restoring metadata failed %s\n',expdir)
%         end
%         
%     end
% end
%% run restore metadata from back up
rootdatadir = '/groups/branson/bransonlab/flydisco_data';
explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/wk40_experimentermetadatafix.txt','%s');
for i = 1:numel(explist)
    expdir = explist{i};
    expdir = fullfile(rootdatadir,expdir)
    metadatafile = fullfile(expdir,'Metadata.xml');
    tmp2 = dir(fullfile(expdir,'Metadata.xml_*'));
    if isempty(tmp2)
        fprintf('no backup metadata file %s\n',expdir);
    else
        if numel(tmp2) > 1
            dn = [];
            %find most recent
            for j = 1:numel(tmp2)
                dn(j) = datenum(tmp2(j).date);
            end
            [~,irecent] = max(dn);
            metadatabackupfile = tmp2(irecent).name;
        else
            metadatabackupfile = tmp2(1).name;
        end
        
        metadatabackupfile = fullfile(expdir,metadatabackupfile);
        
        [success] = RestoreMetadataFromBackup(metadatafile,metadatabackupfile);
        if ~success
            fprintf('restoring metadata failed %s\n',expdir)
        end
        
    end
end