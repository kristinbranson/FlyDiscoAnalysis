function escaped_string = escape_string(raw_string)

cell_array = cellstr(raw_string);
escaped_string = '';
for i = 1:size(cell_array, 1)
  escaped_string = [escaped_string cell_array{i, 1} char(10)]; %#ok
end
escaped_string = escaped_string(1:end-1);
escaped_string = strrep(escaped_string, '\', '\\');
escaped_string = strrep(escaped_string, '''', '\''');