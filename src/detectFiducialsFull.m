% Samuel Rivera
% feb 2, 2011
%
% This function does the KSDA shape detection on a set of testing images.  
% You can give it one at a time if you want
% 
% modelFolder: end with a '/' for KSDA model directory
% trainImages: ( H x W x numSamples )
% trainMarkings, complex ( numLandmarks x numSamples)
% testImages: ( H x W x numSamples )
% saveFile:  a string with output .mat file directory

function [detectedPositionsCell detectTimeCell ]= detectFiducialsFull( ...
                modelFolder, trainImages, trainMarkings, testImages, saveFile)

if ~exist( modelFolder, 'dir')
    mkdir(modelFolder);
end

fast = 0;

tic  %start the clock!
fiducialCoord = getFidCoords( trainMarkings ); % defineds permutation of landmarks
% style = getStyle();   

% initialize the feature types/window crop size/ etc. ---------
[ imgFeatureVector cropSizeMatrix ...
                targetSize modelType] = setCropParams(length( fiducialCoord ) );            
            
            
% pre-allocate and load detections if some exist, expand if cell is not big enough
detectedPositionsCell = cell(length( fiducialCoord ) , size(testImages,3));
detectTimeCell  = cell(length( fiducialCoord ) , size(testImages,3));
if exist( saveFile, 'file')
    load( saveFile )
end
if (size(detectedPositionsCell,1) < length( fiducialCoord )) || ...
        (size(detectedPositionsCell,2) < size(testImages,3))
    detectedPositionsCell{ length(fiducialCoord), size(testImages,3)} = [];  
    detectTimeCell{ length(fiducialCoord), size(testImages,3)} = [];
end
%-----------------------------------------
%-----------------------------------------

% load up all pre-trained models
kSDAModelStruct = cell(length( fiducialCoord )  ,1);
for fidIndex = 1:length( fiducialCoord )  % [1:13 15:length( fiducialCoord ) ]    
    %initialize some things
    fidCoords = fiducialCoord{ fidIndex };
    posNegName = [ modelFolder 'posNegDataFid' int2str(fidIndex) '.mat'];
    kSDAModelFile  = [modelFolder 'kSDAModelFid' int2str(fidIndex)  '.mat'];
    
    %----- generate the training data if needed, or load--------
    if ~exist( kSDAModelFile, 'file')
        if exist( posNegName, 'file')
            load( posNegName )            
        else
            % otherwise generate training data
            cropSize = cropSizeMatrix(:,fidIndex); % [h;w]
            display( ['Generating fiducial ' int2str(fidIndex) ' training data.']);
            [posCropped negCropped ] = generateTrainingData( fidCoords, ...
                trainImages, trainMarkings, cropSize);                                    
            %save( posNegName, 'posCropped', 'negCropped' );   % HERE
        end
    end
    % ------------ train model --------------------        
    % load model if already exists
    if exist( kSDAModelFile, 'file')
        load( kSDAModelFile, 'kSDAModel' )         
    else
        display( ['Training fiducial ' int2str(fidIndex) ' classifier.']);
        imgFeatureMode = imgFeatureVector( fidIndex );  % 1 for unit pixel, 10 for Haar, 9 for center shifted unit pixel
        nc = [ size( posCropped,3) size( negCropped,3) ]; % 1 for class positive, 2 for class negative
        trainStartTime = toc;  
        if modelType(fidIndex) == 1 %1 for kLDA
            [ kSDAModel ] = srTrainKSDA( cat( 3,posCropped, negCropped), ...
                nc, imgFeatureMode, targetSize(:,fidIndex)  );
        elseif modelType(fidIndex) == 2 %1 for adaboost
             [ kSDAModel ] = srTrainKSDA_ADA( cat( 3,posCropped, negCropped), ...
                nc, imgFeatureMode, targetSize(:,fidIndex)  );
        else
                error( 'SR: Unknown model type');
        end
        save( kSDAModelFile, 'kSDAModel' );
        trainEndTime = toc;
        display( [ 'Training time is ' num2str( (trainEndTime-trainStartTime)/60) ' minutes.' ]);
    end    
    % store up the models beforehehand
    kSDAModelStruct{fidIndex} = kSDAModel;
end
%-------------------------------

% loop through each image
for i1 = 1:size( testImages,3)  %1:min( 500,size( testImages,3))

    clf(gcf);
    imagesc( testImages(:,:,i1) ), colormap gray; axis equal; hold on
    display( ['Detecting image ' int2str(i1) '.']);
    % do the detection for each fiducial    
    for fidIndex =  1:length( fiducialCoord ) 
        
        %important initlialize this
        fidCoords = fiducialCoord{ fidIndex };
        
        %------ specify the window search paramters        
        centerPos = mean( trainMarkings(fidCoords,:),1);
        centerPos(isnan(centerPos)) = []; %for missing data
        centerPos = mean( centerPos,2);
        
        centerPos = [imag(centerPos); real(centerPos) ];
        trainFidPositions = mean( trainMarkings(fidCoords,:),1);  % search window params
        xRange = abs( max(real(trainFidPositions)) - min(real(trainFidPositions)));
        yRange = abs( max(imag(trainFidPositions)) - min(imag(trainFidPositions)));
        searchWindowSize = round( 1.2*[ yRange; xRange ] );  % add 20% to the size
        slidingParameters.targetSize =targetSize(:,fidIndex);  
        slidingParameters.idealCropSize = kSDAModelStruct{fidIndex}.oriImageSize;
        slidingParameters.searchWindowSize = searchWindowSize;
        slidingParameters.posnStepSize =  1; 
        slidingParameters.cropSizePerturbations = 1; %[ .9 1 1.1 ];       

        % if not already detected, detect the thing
        if isempty(detectedPositionsCell{fidIndex,i1} )        
            startTime = toc;
            [classEstimates cropParams ] = slidingWindowClassifier( testImages(:,:,i1), ...
                        centerPos, slidingParameters, kSDAModelStruct{fidIndex} );
            % detectedPositions = cropParams(2, classEstimates == 1 ) + 1i*(cropParams(1, classEstimates == 1 ));
            detectedPositions = cropParams(2, classEstimates == 1 )-1 + 1i*(cropParams(1, classEstimates == 1 )-1);
            endTime = toc;      
            detectTime = endTime-startTime;       
            display( ['Fid ' int2str(fidIndex)  ', time: ' num2str( detectTime) ' seconds, ' num2str(length(classEstimates)) ' search positions.' ] );   
            if isempty( detectedPositions );  detectedPositions = -1; end  % parity bid so we don't repeat

            % put the detections in output cell and save
            detectedPositionsCell{fidIndex,i1} = detectedPositions;
            detectTimeCell{fidIndex,i1} = detectTime;                        
        else
            detectedPositions = detectedPositionsCell{fidIndex,i1} ;
        end
        
        % plot the result           
        rectangle('Position', [ centerPos(2)-searchWindowSize(2)/2 ...
            centerPos(1)-searchWindowSize(1)/2 searchWindowSize(2) searchWindowSize(1)] );             
        if detectedPositions ~= -1 % if none detected, parity bit here                
            plot( real(detectedPositions), imag(detectedPositions), 'b+');     %     style{fidIndex}       
        end           
        
        if ~fast % pause to display
            pause(.01); 
        end 
    
    end  % for all fiducials     
    
    
 
        
    save( saveFile, 'detectedPositionsCell', 'detectTimeCell' );
end     % for all images


