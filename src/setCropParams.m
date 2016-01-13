function [ imgFeatureVector cropSizeMatrix targetSize modelType] = setCropParams(numFids)

    originalCropSize = [ 20; 20];  % height; width
    downScaleWindowSize =  [10; 10];   %height; width
    
     
    %----------------------------------------------
    %----------------------------------------------
    originalCropSize = originalCropSize(:);
    downScaleWindowSize = downScaleWindowSize(:);
    imgFeatureVector = 9*ones( 1, numFids);       % 9 for centered unit pixels        
    cropSizeMatrix = [ repmat( originalCropSize, [1,numFids] )]; 
    cropSizeMatrix = round(cropSizeMatrix);  
    targetSize = repmat( downScaleWindowSize, [1, numFids]);



    %modelType is 2 for Adaboost/ 1 for kLDA classifier ( will us kSDA for the
    modelType = 1*ones(1,numFids) ;%repmat( 2, [1, numFids]);    
    

% set to smaller value if necessary
for i1 = 1:size(targetSize,2) 
    if round(cropSizeMatrix(1,i1)) < targetSize(1,i1) 
        targetSize(:,i1) = round(cropSizeMatrix(:,i1));
    end
end
