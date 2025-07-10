function outfilenames = SaveFigLotsOfWays(hfig,basename,formats)

if nargin < 3,
  formats = {'fig','png','pdf','svg'};
end
[path,basename] = myfileparts(basename);

if exist('savefig1')
  savefigfun = @savefig1;
elseif exist('savefig_pa')
  savefigfun = @savefig_pa;
else
  savefigfun = @savefig;
end

outfilenames = {};
if ismember('svg',formats),
  filename = [basename,'.svg'];
  plot2svg(fullfile(path,filename),hfig,'png');
  outfilenames{end+1} = fullfile(path,filename);
end
formats = setdiff(formats,{'svg'});

if ismember('ai',formats),
  filename = [basename,'.ai'];
  saveas(hfig,fullfile(path,filename),'ai');
  outfilenames{end+1} = fullfile(path,filename);
end
formats = setdiff(formats,{'ai'});

if ismember('fig',formats),
  filename = [basename,'.fig'];
  saveas(hfig,fullfile(path,filename),'fig');
  outfilenames{end+1} = fullfile(path,filename);
end
formats = setdiff(formats,{'fig'});

if isequal(class(hfig),'matlab.ui.Figure'),
  hfig = hfig.Number;
end

tmpbasename = regexprep(basename,'[^a-zA-Z_0-9]','_');
for i = 1:numel(formats),
  filename = [basename,'.',formats{i}];
  tmpfilename = [tmpbasename,'.',formats{i}];
  try
    savefigfun(tmpfilename,hfig,formats{i});
    if ~isempty(path),
      movefile(tmpfilename,fullfile(path,filename));
    end
  catch ME,
    warning('Error while creating %s: %s',filename,getReport(ME));
    continue;
  end
  outfilenames{end+1} = fullfile(path,filename);
end