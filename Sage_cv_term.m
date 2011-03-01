function term_lookup = Sage_cv_term(term_name,table)

if nargin < 2,
  table = 'fly_olympiad';
end

% Create the SQL lookup string for the CV term.
term_lookup = ['getCvTermId(''' table ''', ''' escape_string(term_name) ''', NULL)'];