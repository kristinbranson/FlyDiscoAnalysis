function aug_names = compute_flytracker_augmented_jaaba_input_feature_names()
    % Computes the names of FlyTracker augmented JAABA input features, using
    % functions from FlyTracker.  We assume that relative features are not being
    % used, which should be guaranteed by the way we call FlyTracker from FlyDisco
    % code.  Output is a row vector cell array of old-style strings.
    [personal_feat, enviro_feat] = feat_names_and_units() ;  % feat_names_and_units() is part of FlyTracker
    names = [personal_feat enviro_feat];   
    aug_names = feat_augment_names(names) ;  % feat_augment_names() is part of FlyTracker
end
