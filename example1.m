% Samuel Rivera

load( 'trainTestImageMarkingAR.mat', 'trainImages','trainMarkings', 'testImages' );

outputFolder = 'Output/';  
saveFile = [ outputFolder 'detections.mat' ];  % to store the detections

% trainImages:  height x width x N training images
% trainMarkings:  numLandmarks x N complex landmarks  ( x + yi)
% testImages height x width x Ntest images

detectFiducialsFull( outputFolder, trainImages, trainMarkings(2, :), testImages(:,:,1:10), saveFile);