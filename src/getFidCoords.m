%  Samuel Rivera
%  March 7, 2011
% 
% fiducialCoord are the fiducials that are detected using the low level features.
%      It is a cell containing the indices for the original landmark system.
% 
% fiducialCoordAugmented has the extra fiducials which we are going to add
% to that original set of detected fiducials
% 
% associatedFid: indexes from the list 'fiducialCoord' so that we can steal
%   the raw detections from those fiducial detections easily to fill in the
%   augmentedCoordinates
%   
function  [ fiducialCoord  fiducialCoordAugmented  associatedFid ] = getFidCoords( trainMarkings )

    % just 1:N
    fiducialCoord = cell(1,size(trainMarkings,1));
    for i1 = 1:size(trainMarkings,1)
       fiducialCoord{i1} = i1;
    end
      fiducialCoordAugmented = [];
      associatedFid = [];
    
