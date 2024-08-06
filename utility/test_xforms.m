x = normrnd(0,1) 
y = normrnd(0,1) 
offX = normrnd(0,1) 
offY = normrnd(0,1) 
offTheta = normrnd(0,1) 
scale = normrnd(0,1)

Aref = affineTransformMatrixFromOffsetsAndScale(offX,offY,offTheta,scale)

Atest = affineTransformMatricesFromOffsetsAndScale(offX,offY,offTheta,scale)

[xref,yref] = registerForSingleChamber(x,y,offX,offY,offTheta,scale)

[xtest,ytest] = registerForMultipleChambers(x,y,offX,offY,offTheta,scale)
