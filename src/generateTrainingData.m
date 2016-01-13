% Samuel Rivera
% Feb 2, 2011
%
% notes: this function generates training data for particular fiducial
% coordinates.  This just finds perturbed positions, and does not do any
% perturbations in scale or rotation.
% 
% syntax: generateTrainingData( fidCoord, imageMatrix, coordMatrix)
%

function [posCropped negCropped negOffset] = generateTrainingData( fiducialCoord, imageMatrix,...
                                                coordinateMatrix, minCropSize)

% fiducialCoord = 1:13; % 1:13 for left eye
% load( 'AR_Origi.mat', 'coordinateMatrix', 'imageMatrix');  %original true data


% align the images by specific fiducial coordinates 
showCropping = 0;
offset = [];
percentSize = [4; 4] ; % [1.5;1.5]; % [4; 4] for heart and face;
cropSize = findCropSizeUltimate( fiducialCoord, coordinateMatrix, percentSize );
[ X2 Y2 ] = cropFiducialMissing( imageMatrix, coordinateMatrix, fiducialCoord, cropSize, 0, ...
                          offset, showCropping, fiducialCoord );
clear coordinateMatrix imageMatrix                   

% crop the positive based on the centered position     
showCropping = 0;
offset = [];
% minCropSize = findCropSizeUltimate( 1:size(Y2,1), Y2, [1;1] );  % smallest window with whole fiducial                  
posCropped = cropFiducialMissing( X2, Y2, 1:size(Y2,1), minCropSize, 0, ...
                          offset, showCropping, 1:size(Y2,1) );

% ------  crop some negative samples very near to the fiducial
showCropping = 0;

baseOffsetParams = 3.*([ -1 1 -1 1 ] + 1i*[ -1 -1 1 1 ]);  % +- 3 from center
numRepetitions = size( baseOffsetParams,2);
offset = repmat( baseOffsetParams, [ size(Y2,2), 1 ]);
offsetNear = offset(:)';  % offset = offset(1,:) +1i*offset(2,:);
X2near = repmat( X2, [1,1,numRepetitions] );
Y2near = repmat( Y2, [1, numRepetitions] );

negCroppedNear = cropFiducialMissing( X2near, Y2near, 1:size(Y2near,1), minCropSize, 0, ...
                          offsetNear, showCropping, 1:size(Y2near,1) );
clear X2near Y2near

% ---------   crop the negative based on random perturbed position -------

showCropping = 0;

% perturb samples in scale and rotation (for negative samples)
% this will generate many more samples
[X2, Y2]=fbSrPerturbScaleAndRotation( X2, Y2, showCropping);

offsetScalar = 30;
offset = offsetScalar*(rand( 2, size(Y2,2)) -.5);

% make sure not cropping too near the pos fiducial in both x and y coordinate
while sum(  sum(abs(offset),1) < 3 )
    offset( abs(offset) < 2 ) = 3.*offset( abs(offset) < 2 );
end
offset = offset(1,:) +1i*offset(2,:);
negCroppedPerturbed = cropFiducialMissing( X2, Y2, 1:size(Y2,1), minCropSize, 0, ...
                          offset, showCropping, 1:size(Y2,1) );
clear X2 Y2                      
                      
                      
negCropped = cat(3, negCroppedPerturbed, negCroppedNear);


negOffset = [ offsetNear offset ];
