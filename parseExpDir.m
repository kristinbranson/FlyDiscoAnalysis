function [res,success] = parseExpDir(expdir)

res = regexp(expdir,'^(?<pathstr>(.*[/\\])*)(?<line>[^/^\\]+)_(?<effector>[^_]*)_Rig(?<rig>[0-9]+)Plate(?<plate>[0-9]+)Bowl(?<bowl>[A-Z]+)_(?<notstarted>(notstarted_)?)(?<date>[^_/]+)/?$','names','once');
success = ~isempty(res);
if success,
  res.notstarted = strcmp(res.notstarted,'notstarted_');
end
