function MakePerFrameStatsWebpage(filename,basename,groups,statfiles,statperflyfiles,...
  stimfeatures,stimfiles,stimperflyfiles,...
  ions,fliesplot,stimtrajfile,vidfiles,...
  histfiles,histgroups)

fid = fopen(filename,'w');
assert(fid >= 1);
    
fprintf(fid,'<html>\n<title>%s</title>\n<body>\n',basename);
% fprintf(fid,'<head>\n');
% fprintf(fid,'<style>\n');
% fprintf(fid,'table\n');
% fprintf(fid,'{\n');
% fprintf(fid,'border-collapse:collapse;\n');
% fprintf(fid,'}\n');
% fprintf(fid,'table, td, th\n');
% fprintf(fid,'{\n');
% fprintf(fid,'border:1px solid black;\n');
% fprintf(fid,'}\n');
% fprintf(fid,'</style>\n');
% fprintf(fid,'</head>\n');
    

fprintf(fid,'<h1>%s</h1>\n',basename);

% move hist to begining AR 
ishist = ~isempty(histfiles) && ~all(cellfun(@isempty,histfiles));
if ishist,
  fprintf(fid,'<a id="hist"><h2>Histograms</h2></a>\n');
  fprintf(fid,'<center><table width="50%%">\n');
  for i = 1:numel(histfiles),
    name = histfiles{i};
    if isempty(name),
      continue;
    end
    fprintf(fid,'<tr><td><a id="hist_%s" href="%s"><img src="%s" width="100%%"></a></td></tr>\n',histgroups{i},name,name);
  end
  fprintf(fid,'</table>\n');
end


isstatfile = ~isempty(statfiles) && ~all(cellfun(@isempty,statfiles(:)));

if isstatfile,
  fprintf(fid,'<a id="stats"><h2>Mean statistics</h2></a>\n');
  fprintf(fid,'<table width="100%%">\n');
  for i = 1:numel(statfiles),
    fprintf(fid,'<tr>\n');
    name = statfiles{i};
    if ~isempty(name),
      fprintf(fid,'  <td><a id="stat_%s" href="%s"><img src="%s" width="100%%"></a></td>\n',groups{i},name,name);
    end
    if ~isempty(statperflyfiles),
      name = statperflyfiles{i};
      fprintf(fid,'  <td><a href="%s"><img src="%s" width="100%%"></a></td>\n',name,name);
    end
    fprintf(fid,'</tr>\n');
  end
  fprintf(fid,'</table>\n');
end

isstimfile = ~isempty(stimfiles) && ~all(cellfun(@isempty,stimfiles(:)));
if isstimfile,
  fprintf(fid,'<a id="stim"><h2>Mean time series during stimulation</h2></a>\n');
  
  fprintf(fid,'<table width="100%%">\n');
  fprintf(fid,'<tr><th></th>');
  for i = 1:size(stimfiles,1),
    fprintf(fid,'<th>Period %d</th> ',i);
  end
  fprintf(fid,'</tr>\n');
  for fi = 1:size(stimfiles,2),
    fprintf(fid,'<tr>\n');
    fprintf(fid,'<th><span style="writing-mode: vertical-lr; -ms-writing-mode: tb-rl; transform: rotate(180deg);">%s</span></th>\n',stimfeatures{fi});
    %fprintf(fid,'<a id="stim_%s"><h3>%s</h3></a>\n',stimfeatures{fi},stimfeatures{fi});
    for ion = 1:size(stimfiles,1),
      name = stimfiles{ion,fi};
      if isempty(name),
        continue;
      end
      fprintf(fid,'  <td><a id="stim_%d_%s" href="%s"><img src="%s" width="100%%"></a></td>\n',ion,stimfeatures{fi},name,name);
    end
    fprintf(fid,'</tr>\n');
  end
  fprintf(fid,'</table>\n');
end

plotflies = ~isempty(stimperflyfiles) && ~all(cellfun(@isempty,stimperflyfiles(:)));

if plotflies,
  fprintf(fid,'<a id="stimperfly"><h2>Per-fly time series during stimulation</h2></a>\n');
  
  fprintf(fid,'<table width="100%%">\n');
  fprintf(fid,'<tr><th></th>');
  for i = 1:size(stimfiles,1),
    fprintf(fid,'<th>Period %d</th> ',i);
  end
  fprintf(fid,'</tr>\n');
  for fi = 1:size(stimfiles,2),
    fprintf(fid,'<tr>\n');
    fprintf(fid,'<th><span style="writing-mode: vertical-lr; -ms-writing-mode: tb-rl; transform: rotate(180deg);">%s</span></th>\n',stimfeatures{fi});
    %fprintf(fid,'<a id="stim_%s"><h3>%s</h3></a>\n',stimfeatures{fi},stimfeatures{fi});
    for ion = 1:size(stimfiles,1),
      name = stimperflyfiles{ion,fi};
      if isempty(name),
        continue;
      end
      fprintf(fid,'  <td><a id="stimperfly_%d_%s" href="%s"><img src="%s" width="100%%"></a></td>\n',ion,stimfeatures{fi},name,name);
    end
    fprintf(fid,'</tr>\n');
  end
  fprintf(fid,'</table>\n');
end

if ~isempty(stimtrajfile),
  fprintf(fid,'<p>Trajectories at the onset of activation</p>\n');
  fprintf(fid,'<p><a id="stimtraj" href="%s"><img src="%s" width="100%%"></a></p>\n',stimtrajfile,stimtrajfile);
end


% ishist = ~isempty(histfiles) && ~all(cellfun(@isempty,histfiles));
% if ishist,
%   fprintf(fid,'<a id="hist"><h2>Histograms</h2></a>\n');
%   fprintf(fid,'<center><table width="50%%">\n');
%   for i = 1:numel(histfiles),
%     name = histfiles{i};
%     if isempty(name),
%       continue;
%     end
%     fprintf(fid,'<tr><td><a id="hist_%s" href="%s"><img src="%s" width="100%%"></a></td></tr>\n',histgroups{i},name,name);
%   end
%   fprintf(fid,'</table>\n');
% end

if ~isempty(vidfiles),
  fprintf(fid,'<a id="stim"><h2>Stimulation onset videos</h2></a>\n');
  fprintf(fid,'<table width="100%%">\n');
  %fprintf(fid,'<tr><th></th>');
  colw = 100/size(vidfiles,1);
  for i = 1:size(vidfiles,1),
    fprintf(fid,'<th width="%f%%">Period %d</th> ',colw,i);
  end
  fprintf(fid,'</tr>\n');
  for fi = 1:size(vidfiles,2),
    fprintf(fid,'<tr>\n');
    %fprintf(fid,'<th><span style="writing-mode: vertical-lr; -ms-writing-mode: tb-rl; transform: rotate(180deg);">Fly %d</span></th>\n',fliesplot(fi));
    %fprintf(fid,'<a id="stim_%s"><h3>%s</h3></a>\n',stimfeatures{fi},stimfeatures{fi});
    for ion = 1:size(vidfiles,1),
      name = vidfiles{ion,fi};
      if isempty(name),
        continue;
      end
      fprintf(fid,'  <td><a id="stimvideo_period%d_fly%d" href="%s"><img src="%s" width="100%%"></a></td>\n',ions(ion),fliesplot(fi),name,name);
    end
    fprintf(fid,'</tr>\n');
  end
  fprintf(fid,'</table>\n');
end

fprintf(fid,'</body>\n</html>\n');
fclose(fid);