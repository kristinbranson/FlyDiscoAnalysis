function LoadIndicatorDataFromFile(obj,n)

try
  indicatordatafile = fullfile(obj.expdirs{n},obj.dataloc_params.indicatordatafilestr);
  if ~exist(indicatordatafile,'file'),
    return;
  end
  id = load(indicatordatafile);
catch ME
  warning(getReport(ME));
  return;
end
indicatorLED = convertIndicatorData(id.indicatorLED) ;
obj.indicatorLED{n} = indicatorLED;

end  % function



function result = convertIndicatorData(led)
% Converts the indicator LED onset/offset data in the form it is held in the
% data file to the form we want it in.  This is a pure function.

result = struct;
result.starton = [];
result.endon = [];
result.startoff = [];
result.endoff = [];

result.starton = led.startframe;
result.endon = led.endframe;
% if led.StartEndStatus(1),
%   result.starton = [1,result.starton];
% end
% if led.StartEndStatus(2),
%   result.endon(end+1) = obj.movie_nframes(n);
% end
result.startoff = led.endframe+1;
result.endoff = led.startframe-1;
% if led.StartEndStatus(1)==0,
%   result.startoff = [1,result.startoff];
% end
% if led.StartEndStatus(2)==0,
%   result.endoff(end+1) = obj.movie_nframes(n);
% end
assert(numel(result.starton) == numel(result.endon));
assert(numel(result.startoff) == numel(result.endoff));

end  % function
