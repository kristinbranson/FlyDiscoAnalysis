function LoadIndicatorDataFromFile(obj,n)

id = load(fullfile(obj.expdirs{n},obj.dataloc_params.indicatordatafilestr));
indicatorLED = struct;
indicatorLED.starton = id.indicatorLED.startframe;
indicatorLED.endon = id.indicatorLED.endframe;
if id.indicatorLED.StartEndStatus(1),
  indicatorLED.starton = [1,indicatorLED.starton];
end
if id.indicatorLED.StartEndStatus(2),
  indicatorLED.endon(end+1) = obj.movie_nframes(n);
end
indicatorLED.startoff = id.indicatorLED.endframe+1;
indicatorLED.endoff = id.indicatorLED.startframe-1;
if id.indicatorLED.StartEndStatus(1)==0,
  indicatorLED.startoff = [1,indicatorLED.startoff];
end
if id.indicatorLED.StartEndStatus(2)==0,
  indicatorLED.endoff(end+1) = obj.movie_nframes(n);
end
assert(numel(indicatorLED.starton) == numel(indicatorLED.endon));
assert(numel(indicatorLED.startoff) == numel(indicatorLED.endoff));
obj.indicatorLED{n} = indicatorLED;
