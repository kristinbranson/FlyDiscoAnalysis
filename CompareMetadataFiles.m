function CompareMetadataFiles(file1,file2)

metadata1 = ReadMetadataFile(file1);
metadata2 = ReadMetadataFile(file2);
fns1 = fieldnames(metadata1);
fns2 = fieldnames(metadata2);

for i = 1:numel(fns1),
  fn = fns1{i};
  if ~ismember(fn,fns2),
    fprintf('%s missing in %s\n',fn,file2);
    continue;
  end
  if ischar(metadata1.(fn)) ~= ischar(metadata2.(fn)),
    fprintf('%s: data types do not match',fn);
    continue;
  end
  if ischar(metadata1.(fn)),
    if ~strcmp(metadata1.(fn),metadata2.(fn)),
      fprintf('%s: %s ~= %s\n',fn,metadata1.(fn),metadata2.(fn));
    end
  else
    if ~all(metadata1.(fn) == metadata2.(fn)),
      fprintf('%s: %s ~= %s\n',fn,mat2str(metadata1.(fn)),mat2str(metadata2.(fn)));
    end
  end
end

missingfns = setdiff(fns2,fns1);
for i = 1:numel(missingfns),
  fprintf('%s missing in %s\n',missingfns{i},file1);
end