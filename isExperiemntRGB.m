function result = isExperiemntRGB(metadata)

result = ( strcmp(metadata.assay,'FlyBubbleRGB') || strcmp(metadata.assay,'FlyBowlRGB') || ...
    strcmp(metadata.assay,'FlyBowlFRGB') || strcmp(metadata.assay,'FlyBubbleFRGB')  ) ;
end