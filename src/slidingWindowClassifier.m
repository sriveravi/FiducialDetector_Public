% Samuel Rivera
% Feb 3, 2011
% notes: This function is main wrapper for sliding window classifier


function [classEstimates cropParams ] = slidingWindowClassifier( testImage, centerPos, ...
                                            slidingParameters, kSDAModel )

% load( 'AR_Origi.mat', 'coordinateMatrix', 'testImagesMatrix');  %original true data

% posnStepSize =2;  integet specifying the window position stepsize
% cropSizePerturbations = [ .8 1 1.2]; % scale the idealCropSize by this
                                         % to get new crop sizes
% idealCropSize = [ 10; 20];  in [ heigh; width ]
% targetSize = %size for classifier  in [ heigh; width ]
% searchWindowSize = [10; 15];  in [ heigh; width ]
% centerPos   in [ heigh; width ]

% center position given as mean of true coords for preliminary work
% centerPos = mean( coordinateMatrix(desiredFid,i1),1 );
% centerPos = [ imag(centerPos); real(centerPos) ];
%
% kSDAModel is a struct containing the KSDA model data

idealCropSize = slidingParameters.idealCropSize(:);
targetSize = slidingParameters.targetSize(:); 
searchWindowSize= slidingParameters.searchWindowSize(:);
posnStepSize  = slidingParameters.posnStepSize(:);
cropSizePerturbations= slidingParameters.cropSizePerturbations(:);

% crop the possible positions
[ cropImages cropParams ] = sampleRegionsInImage( testImage, ...
                    centerPos, searchWindowSize, idealCropSize, ...
                    cropSizePerturbations, posnStepSize, targetSize ); 


% [ kSDAModel ] = srTrainKSDA( trainingImages, nc, imgFeatureMode )


% run classifier
test_label = ones( 1, size(cropImages,3));

if kSDAModel.type == 1 % 1 for kLDA
    [ classEstimates ] = srTestKSDA( kSDAModel, cropImages, test_label );
elseif kSDAModel.type == 2 % 2 for adaboost
    [ classEstimates ] = srTestKSDA_ADA( kSDAModel, cropImages, test_label );
end
                   
