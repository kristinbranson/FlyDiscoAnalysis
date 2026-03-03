rootdatadir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlData';
sciservdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

old_line_names = {'GMR_MB101B','GMR_MB161B','GMR_MB044B'};
new_line_names = {'GMR_MB316B','GMR_MB077B','GMR_MB315B'};

expnames = dir(fullfile(rootdatadir,'*20*'));
expnames(~[expnames.isdir]) = [];
expnames = {expnames.name};

for i = 1:numel(old_line_names),
  
  ismatch = ~cellfun(@isempty,regexp(expnames,['^',old_line_names{i}],'once'));
  idxmatch = find(ismatch);
  for jj = 1:numel(idxmatch),
    j = idxmatch(jj);
    newexpname = regexprep(expnames{j},['^',old_line_names{i}],new_line_names{i});
    inexpdir = fullfile(rootdatadir,expnames{j});
    outexpdir = fullfile(rootdatadir,newexpname);
    insciservexpdir = fullfile(sciservdir,expnames{j});
    outsciservexpdir = fullfile(sciservdir,newexpname);
    
    % fix all broken soft-links
    files = dir(inexpdir);
    for k = 1:numel(files),
      if ismember(files(k).name,{'.','..'}),
        continue;
      end
      filename = fullfile(inexpdir,files(k).name);
      if exist(filename,'file'),
        continue;
      end
      [res,s] = unix(sprintf('readlink %s',filename));
      s = strtrim(s);
      if res ~= 0 || isempty(s),
        continue;
      end
      [~,name] = myfileparts(s);
      cmd = sprintf('rm %s',filename);
      unix(cmd);
      newlink = fullfile(outsciservexpdir,name);
      cmd = sprintf('ln -s %s %s',newlink,filename);
      unix(cmd);
      if ~exist(filename,'file') && ~isempty(regexp(name,old_line_names{i},'once')),
        name = regexprep(name,old_line_names{i},new_line_names{i});
        filename = fullfile(inexpdir,name);
        newlink = fullfile(outsciservexpdir,name);
        if exist(newlink,'file'),
          fprintf('Renaming file %s to %s\n',s,newlink);
          cmd = sprintf('ln -s %s %s',newlink,filename);
          unix(cmd);
        end
      end
     
      if ~exist(filename,'file')        
        error('File %s still does not exist',filename);
      end
    end
    
    files = dir(fullfile(inexpdir,'perframe'));
    for k = 1:numel(files),
      if ismember(files(k).name,{'.','..'}),
        continue;
      end
      filename = fullfile(inexpdir,'perframe',files(k).name);
      if exist(filename,'file'),
        continue;
      end
      [res,s] = unix(sprintf('readlink %s',filename));
      s = strtrim(s);
      if res ~= 0 || isempty(s),
        continue;
      end
      [~,name] = myfileparts(s);
      cmd = sprintf('rm %s',filename);
      unix(cmd);
      newlink = fullfile(outsciservexpdir,'perframe',name);
      cmd = sprintf('ln -s %s %s',newlink,filename);
      unix(cmd);
      if ~exist(filename,'file') && ~isempty(regexp(name,old_line_names{i},'once')),
        name = regexprep(name,old_line_names{i},new_line_names{i});
        filename = fullfile(inexpdir,'perframe',name);
        newlink = fullfile(outsciservexpdir,'perframe',name);
        if exist(newlink,'file'),
          fprintf('Renaming file %s to %s\n',s,newlink);
          cmd = sprintf('ln -s %s %s',newlink,filename);
          unix(cmd);
        end
      end
    end
      
    
    % rename directory
    cmd = sprintf('mv %s %s',inexpdir,outexpdir);
    unix(cmd);
  
  end
  
end
  
  
  
