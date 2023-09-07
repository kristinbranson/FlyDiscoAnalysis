function extractRegistrationFunctionAndTransform(registration)
  registration.registerfn = @(x,y) register(x,y,registration.offX,...
    registration.offY,registration.offTheta,registration.scale);
  registration.affine = affineTransform(registration.offX,...
    registration.offY,registration.offTheta,registration.scale);
end

