function result = isExperimentRGB(metadata)
    result = ( strcmp(metadata.assay,'FlyBubbleRGB') || strcmp(metadata.assay,'FlyBowlRGB') || ...
               strcmp(metadata.assay,'FlyBowlFRGB') || strcmp(metadata.assay,'FlyBubblefRGB') ) ;
end
