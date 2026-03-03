function [bIsImageJBigStack, bIsImageJHyperStack, vnStackDims, vnInterleavedFrameDims] = IsImageJBigStack(sInfo, nAparrentSize)

   % - Set up default return arguments
   bIsImageJBigStack = false;
   bIsImageJHyperStack = false;
   vnStackDims = [];
   vnInterleavedFrameDims = [];
   
   % - Check for ImageDescription field
   if (~isfield(sInfo, 'ImageDescription'))
      return;
   end
   
   % - Get image description
   strImageDesc = sInfo(1).ImageDescription;
   
   % - Look for ImageJ version information
   strImageJVer = sscanf(strImageDesc(strfind(strImageDesc, 'ImageJ='):end), 'ImageJ=%s');
   
   % - Look for stack size information
   if (~isempty(strImageJVer))
      nNumImages = sscanf(strImageDesc(strfind(strImageDesc, 'images='):end), 'images=%d');
      
      % - Does ImageJ report a greater number of images than sInfo?
      if (~isempty(nNumImages) && (nAparrentSize ~= nNumImages))
         bIsImageJBigStack = true;
      end
      
      % - Is this a hyperstack?
      strHyperStack = sscanf(strImageDesc(strfind(strImageDesc, 'hyperstack='):end), 'hyperstack=%s');
      
      if (strcmpi(strHyperStack, 'true'))
         bIsImageJHyperStack = true;
         
         % - Extract information about the stack size for a hyperstack
         nNumChannels = sscanf(strImageDesc(strfind(strImageDesc, 'channels='):end), 'channels=%d');
         nNumSlices = sscanf(strImageDesc(strfind(strImageDesc, 'slices='):end), 'slices=%d');
         nNumFrames = sscanf(strImageDesc(strfind(strImageDesc, 'frames='):end), 'frames=%d');
         
         if (isempty(nNumChannels))
            nNumChannels = 1;
         end
         
         if (isempty(nNumSlices))
            nNumSlices = 1;
         end
         
         if (isempty(nNumFrames))
            nNumFrames = 1;
         end
         
         % - Deinterleave stack
         vnStackDims = [sInfo(1).Width sInfo(1).Height nNumFrames*nNumSlices*nNumChannels 1];
         vnInterleavedFrameDims = [nNumChannels nNumSlices nNumFrames];
         
      else
         % - Extract information about the stack size for a fake big stack
         vnStackDims = [sInfo(1).Width sInfo(1).Height nNumImages sInfo(1).SamplesPerPixel];
         vnInterleavedFrameDims = [];
      end
   end
end
