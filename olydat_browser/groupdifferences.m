function [h p grpDat grpleg0] = groupdifferences(hFigOrAx,dat,grp,datname,grpname,titleLbl)
% [h p] = groupdifferences(hFigOrAx,dat,grp,datname,grpname,titleLbl)
% Examine differences of data across groups.
%
% Group data, plot distributions for each group, compute KS statistic
% matrix between groups.
%
% dat: column double vec, data of interest.
% grp: int or cellstr grouping vec, same size as dat (see groupdata)
% datname: label for data
% grpname: label for group
% titleLbl: label for title.
%


assert(numel(dat)==numel(grp));
KSALPHA = 0.01; % KS significant threshold

[grpDat grpleg] = groupdata(dat,grp);

% get rid of groups with < 5 els
tfFewerThanFive = cellfun(@(x)numel(x)<5,grpDat);
if any(tfFewerThanFive)
  if isnumeric(grpleg),
    s = sprintf(' %d',grpleg(tfFewerThanFive));
  else
    s = sprintf(' %s',grpleg{tfFewerThanFive});
  end
  warning('groupdifferences:fewerThanFiveDataPoints',...
    'Ignoring groups with fewer than five data points for:%s.',s);
end
grpDat = grpDat(~tfFewerThanFive);
grpleg = grpleg(~tfFewerThanFive);

grpleg0 = grpleg;
if isnumeric(grpleg)
    grpleg = cellstr(num2str(grpleg(:)));
end

if strcmp(get(hFigOrAx,'Type'),'figure')
    clf(hFigOrAx);
    ax = gca;
else
    ax = hFigOrAx;
end
distplot(ax,grpDat,[],[],grpleg);
ylabel(datname,'interpreter','none');
xlabel(grpname,'interpreter','none');

if ~isempty(grpname)
    title(sprintf('%s: %s by %s',titleLbl,datname,grpname),'interpreter','none');
else
    title(sprintf('%s: %s',titleLbl,datname),'interpreter','none');
end

if isempty(grpDat)
    h = [];
    p = [];
else
    [h p] = ksmatrix(grpDat,KSALPHA);
end

end
