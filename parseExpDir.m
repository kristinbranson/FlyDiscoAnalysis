function [res,success] = parseExpDir(expdir)

effectors = {'TrpA','CTRL_CantonS_1101243'};
expr = ['^(?<pathstr>(.*[/\\])*)(?<line>[^/^\\]+)_(?<effector>(([^_]*)',sprintf('|(%s)',effectors{:}),'))_Rig(?<rig>[0-9]+)Plate(?<plate>[0-9]+)Bowl(?<bowl>[A-Z]+)_(?<notstarted>(notstarted_)?)(?<date>[^_/]+)/?$'];
res = regexp(expdir,expr,'names','once');
success = ~isempty(res);
if success,
  res.notstarted = strcmp(res.notstarted,'notstarted_');
end
