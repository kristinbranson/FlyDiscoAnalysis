function registration = extractRegistrationFunctionAndTransform(registration)
  registration.registerfn = @(x,y) registerForMultipleChambers(x,y,registration.offX,...
    registration.offY,registration.offTheta,registration.scale);
  registration.affine = affineTransformMatricesFromOffsetsAndScale(registration.offX,...
    registration.offY,registration.offTheta,registration.scale);
end

