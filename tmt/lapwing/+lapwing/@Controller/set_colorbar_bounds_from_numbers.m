function set_colorbar_bounds_from_numbers(self,cb_min,cb_max)

% Set the view colorbar bounds given max and min values as numbers.  This
% calls self.set_colorbar_bounds_from_strings().  It is assumed that cb_min
% and cb_max are doubles that happen to be integral.

% convert to strings
cb_min_string=sprintf('%d',cb_min);
cb_max_string=sprintf('%d',cb_max);

% call the method that does the real work
self.set_colorbar_bounds_from_strings(cb_min_string,cb_max_string);

end
