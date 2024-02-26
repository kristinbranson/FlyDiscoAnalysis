function extractRegistrationFunctionAndTransform(registration)
  registration.registerfn = @(x,y) registerForSingleAffine(x,y,registration.offX,...
    registration.offY,registration.offTheta,registration.scale);
  registration.affine = affineTransformMatrixFromOffsetsAndScale(registration.offX,...
    registration.offY,registration.offTheta,registration.scale);
end

