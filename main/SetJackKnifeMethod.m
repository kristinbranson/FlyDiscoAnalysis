function jackknife = SetJackKnifeMethod(jackknife,nflies,nexpdirs)

if isempty(jackknife),
  if nexpdirs > 1,
    jackknife = 'perexp';
  elseif nflies == 1,
    jackknife = 'none';
  else
    jackknife = 'perfly';
  end
end

% to split over experiments, we need more than one experiment
% if there isn't set to perfly
if nexpdirs == 1 && strcmpi(jackknife,'perexp'),
  warning('Only one experiment selected, but jackknife = ''perexp''. Splitting over flies instead.');
  jackknife = 'perfly';
end
% to split over flies, we need more than one fly. if there isn't set to
% none
if nflies == 1 && strcmpi(jackknife,'perfly'),
  warning('Only one fly selected, but jackknife = ''perfly''. Splitting over frames not implemented. Not jackknifing.');
  jackknife = 'none';
end

