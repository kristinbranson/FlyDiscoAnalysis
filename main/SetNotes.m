function SetNotes(data,db)

%if isempty(db), return; end

experiment_id = data.experiment_id;
property_names = {'notes_behavioral','notes_technical'};
property_values = {data.notes_behavioral,data.notes_technical};

%% do the update

for i = 1:numel(property_names),
  property_name = property_names{i};
  property_value = property_values{i};

  query = sprintf('select id from experiment_property where experiment_id = %d and type_id = %s',experiment_id,Sage_cv_term(property_name));
  if isempty(db),
    fprintf([query,'\n']);
    property_id = nan;
  else
    curs = exec(db, query);
    curs = fetch(curs);
    close(curs);
    property_id = curs.Data{1};
  end
  
  if strcmp(property_id,'No Data'),
    query = ['insert into experiment_property (experiment_id, type_id, value) values (' num2str(experiment_id) ', ' Sage_cv_term(property_name) ', ''' escape_string(property_value) ''')'];
  else
    query = sprintf('update experiment_property set value = ''%s'' where id = %d',escape_string(property_value),property_id);
  end
  
  if isempty(db),
    fprintf([query,';\n']);
  else
    curs = exec(db, query);
    if ~isempty(curs.Message)
      error(['Could not update experiment #' num2str(experiment_id) ' ''' property_name ''' property: ' curs.Message]);
    end
    close(curs);
  end
  
end