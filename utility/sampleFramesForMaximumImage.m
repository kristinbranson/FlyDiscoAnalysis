function result = sampleFramesForMaximumImage(readfcn, firstFrameIndex, stride, lastFrameIndex)
result = readfcn(1) ;
for j = firstFrameIndex:stride:lastFrameIndex
  frame = readfcn(j) ;
  result = max(result,frame) ;
end
end
