function [isproblem,notes] = ExamineVideoDiagnostics(data,expidx,varargin)

[isproblem,notes] = myparse(varargin,...
  'isproblem',nan(1,numel(data)),...
  'notes',cell(1,numel(data)));

notefns = {'notes_technical','notes_behavioral','notes_curation','notes_keyword'};
[~,~,setidx] = unique({data.set});

for ii = 1:numel(expidx),
  i = expidx(ii);
  
  % show the video diagnostics fimage
  imfilename = fullfile(data(i).file_system_path,'video_diagnostics.png');
  clf;
  imshow(imfilename);
  
  % put some diagnostic text on it
  ax = axis;
  mn = {};
  for j = 1:numel(notefns),
    if ~isempty(data(i).(notefns{j})) && ~ismember(data(i).(notefns{j}),{'None','NULL'}),
      mn{end+1} = sprintf('%s: %s',notefns{j},data(i).(notefns{j})); %#ok<AGROW>
    end
  end
  idxcurr = find(setidx == setidx(i) & isproblem == 1);
  for j = idxcurr,
    mn{end+1} = sprintf('%s: %s',data(j).experiment_name(9:end),notes{j}); %#ok<AGROW>
  end
  if ~isempty(mn),
    text(ax(2),ax(4),mn,'HorizontalAlignment','right','VerticalAlignment','bottom',...
      'BackgroundColor','k','Color','w','Interpreter','none');
  end
  
  % title
  [~,name] = fileparts(data(i).file_system_path);
  hti = title(name);
  set(hti,'Interpreter','none');
  
  res = questdlg('Is there a problem?','Problem?','Yes','See UFMF','No','See UFMF');
  
  if strcmpi(res,'See UFMF'),
    uiwait(showufmf('UFMFName',fullfile(data(i).file_system_path,'movie.ufmf')));
    res = questdlg('Is there a problem?','Problem?');
  end
  
  if strcmpi(res,'Yes'),
    isproblem(i) = 1;
    
    res = inputdlg('Notes: ','Add note',1);
    if isempty(res),
      return;
    end
    if iscell(res),
      res = res{1};
    end
    if ~isempty(res),
      notes{i} = res;
    end

    
  elseif strcmpi(res,'No'),
    isproblem(i) = 0;
  else
    fprintf('Annotated %d experiments\n',ii-1);
    return;
  end
  
  res = questdlg('Next experiment?','Next experiment','Yes','Quit','Yes');
  if strcmpi(res,'Quit'),
    return;
  end
  
end

fprintf('Annotated ALL %d experiments\n',numel(expidx));