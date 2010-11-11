function res = parseExpDir(expdir)

res = regexp(expdir,'(?<line>[^/]+)_TrpA_Rig(?<rig>[0-9]+)Plate(?<plate>[0-9]+)Bowl(?<bowl>[A-Z]+)_(?<notstarted>(notstarted_)?)(?<date>[^_/]+)/?$','names','once');
if ~isempty(res),
  res.notstarted = strcmp(res.notstarted,'notstarted_');
end