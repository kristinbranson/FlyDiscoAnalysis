load detect-circles-2.mat

input = [ zeros([512 1024]) zeros([512 1024]) ; ...
          zeros(size(circleim)) circleim              ] ;

imglance(input)
title('input') ;

[circleRadius,circleCenterX,circleCenterY,featureStrengths,circleDetectParams] = ...
  detectcircles(input,...
                'cannythresh',circleCannyThresh,'cannysigma',circleCannySigma,...
                'binedgesa',1024+binedgesa,'bincentersb',512+bincentersb,'bincentersr',bincentersr,...
                'maxncircles',1,'doedgedetect',strcmpi(circleImageType,'canny'));

circleRadius
circleCenterX
circleCenterY
