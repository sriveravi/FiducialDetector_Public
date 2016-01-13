% Samuel Rivera
% Feb 3, 2011
% Notes: This function will input an image, and output all the possible
% crops for doing a further classification
% 
% centerPos = [ height; width ]
% searchSize = [ height; width ] of maximal area over which to look
% idealCropSize = [h; w] of the expected size of the object
% 
% targetSize  = [h;w] to scale to for classification

function [ cropImages cropParams ] = sampleRegionsInImage( inputImage, ...
                        centerPos, searchSize, idealCropSize, cropSizePerturbations, ...
                        posnStepSize, targetSize )

showCropping = 0;

cropImages = [];
cropParams = [];

targetSize = round(targetSize(:));
idealCropSize = round(idealCropSize);                   
numScales = length( cropSizePerturbations );


for i1 = 1:numScales
    
    cropSize =  round(idealCropSize.*cropSizePerturbations(i1));    
    
    % generate all possible offsets from center for cropping the fiducial
    [ yOffset, xOffset] = ndgrid( 1:posnStepSize:searchSize(1), 1:posnStepSize:searchSize(2));
    yOffset = yOffset(:) - searchSize(1)/2;
    xOffset = xOffset(:) - searchSize(2)/2;
    numCrops = length(yOffset); 
    offset = (xOffset+ 1i.*yOffset )';
    
    [ X2 Y2 oriCropPosn] = cropFiducial( repmat(inputImage, [1,1,numCrops]), ...
                                repmat( centerPos(2)+1i*centerPos(1), [1, numCrops] ), ...
                                1, cropSize, 0, offset, showCropping, centerPos );   
    clear Y2
    
    % find
    % in [ y; x] format
    centerCropPosn = [imag(oriCropPosn); real(oriCropPosn)] ...
                      + repmat(.5*cropSize(:), [1, size(oriCropPosn,2)] );
    
    
    % resize images if necessary
    if cropSize(:) ~= targetSize(:)
        X2rightSize = zeros( targetSize(1), targetSize(2), numCrops);
        for k1 = 1:numCrops
            X2rightSize(:,:,k1) = imresize( X2(:,:,k1), [targetSize(1) targetSize(2)] );
        end
    else
        X2rightSize = X2;
    end
      
    % need to add a check, since imresize does not work properly sometimes
    % :(
    
    

    
    % stack images/centers for output
    if i1 == 1  % for first one
        cropImages = X2rightSize;
    else
        cropImages = cat( 3, cropImages, X2rightSize);
    end
    cropParams = [ cropParams centerCropPosn ];
    
end

