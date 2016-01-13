

function [ nX ] = formatImages( parameterFile, X, doItFast, imgFeatureMode, trainThis, testThis) 

% parameterFile: stores the PCAcomponents/parameters for the normalization
% X: matrix with images to read
% trainThis: vector specifying which are the training images (will be used to calculate PCA)
% doItFast: will check if parameters already stores, and use them if they exist
%
% imgFeatureMode = 1 for unit length pixels
%                  2 for C1 features
%                  3 for LBP ( local binary pattern)
%                  4 for Berkeley Wavelet Transform
%                  5 for unit PCA
%                  6 for Laplacian Eigenmaps
%                  7 for Pixel
%                  8 for PCA
%                  9 for center shifted unit pixel
%                 10 for Haar features
%                 11 for Gabor features (not implemented)
%                 12 for Self Similarity descriptors      
% 
% imageDimensions = hxw of image for reshape and getting the C1 features if
% the images are not already in vectorized format, I'll check of course

opts.disp=0;
xMu = []; 
Vimg = []; 
Dimg = [];

imagePercent = .99;

%-----------------------------------------------------------------

% vectorize images if necessary
trainThis = trainThis(:);
testThis = testThis(:);
s = size( X );
if  size( size(X),2) == 3 || trainThis==1, % vectorize images if necessary
    if trainThis==1,
        s(3)=1;
    end
    X = reshape( X, size( X,1)*size(X,2), s(3) );
    imageDimensions = [s(1), s(2)];
end

N = size(X,2);

% feature selection
if imgFeatureMode == 1 %unit normalized pixel intensity
    nX =  normc(X ); 
    
elseif imgFeatureMode == 2 %C1 features
    c1file = [ 'C1File.mat' ];
    [r c lib] = calcC1( c1file, reshape( X(:,1) ,imageDimensions ));
    X2 = zeros( size(r,1), N);
    if ~exist( c1file, 'file' )
        save( c1file, 'c', 'lib' );
    end
    for i2 = 1:N
        [r c lib] = calcC1( c1file, reshape( X(:,i2) ,imageDimensions));
        X2(:,i2) = r;
    end
    nX = X2; 
    
    clear X2%cleanup
%     delete(c1file)

    X =  normc(nX ); 
    needLearnParams = 1;
    if doItFast && exist( parameterFile, 'file') 

        load( parameterFile, 'xMu', 'Vimg', 'Dimg' );
        needLearnParams = 0;
        if size( X,1) ~= size( xMu,1)
            needLearnParams = 1;
            clear  xMu Vimg Dimg
        else
            nX = X - repmat( xMu, 1, N );
            nX =  diag(diag(Dimg).^(-1/2))*Vimg'*nX;
        end
    end

    if needLearnParams  % whiten and reduce
            xMu = mean(X( :, trainThis),2);
            [d N] = size(X); 
            Xc = X - repmat( xMu, 1, N ); 
            [V D ] = pcaSR(Xc(:,trainThis)', imagePercent);
            Dimg = D;
            Vimg = V;
            nX = diag(diag(Dimg).^(-1/2))*Vimg'*Xc;     
            % display( ['num PCA image modes is ' int2str(size(V,2)) ' out of ' ...
            %     int2str( min(size(Xc))) ]);
    end
    
    nX = normc(nX);
    
    % save PCA parameters
    if ~isempty( parameterFile )  
        save(  parameterFile,  'xMu', 'Vimg', 'Dimg', 'imageDimensions', '-v7.3' );
    end
    
    
elseif imgFeatureMode == 3 %LBP
    %LBP code image using sampling points in SP, no mapping
    SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    temp = lbp( reshape( X(:,1),imageDimensions ),SP, 0,'i'); 
    X2 = zeros( size(temp(:),1), N);
    for i2 = 1:N
        temp = lbp( reshape( X(:,i2), imageDimensions ),SP,0,'i');
        X2(:,i2) = temp(:);     
    end
    nX = X2;  %normc(X2);
    clear X2 temp
    
elseif  imgFeatureMode == 4 %BWT, berkeley wavelet transform
                            %make sure image square
                            % and length a power of 3
    n = max( 3, round( log( min( imageDimensions  ) )/log( 3) ) );                    
    tempImg = imresize( reshape( X(:,1),imageDimensions ), ...
                                     [3^n 3^n ] );
    temp = bwt( tempImg ); 
    X2 = zeros( size(temp(:),1), N);
    for i2 = 1:N
        tempImg = imresize( reshape( X(:,i2),imageDimensions ), ...
                                     [3^n 3^n ] );
        temp = bwt( tempImg );
        X2(:,i2) = temp(:);     
    end
    nX = X2;  
    
elseif  ( imgFeatureMode ==5 ) ||  ( imgFeatureMode ==8 )   % 5 whiten unit length,
                                                            % 8 no unit length 
    X =  normc(X ); 
    needLearnParams = 1;
    if doItFast && exist( parameterFile, 'file') 

        load( parameterFile, 'xMu', 'Vimg', 'Dimg' );
        needLearnParams = 0;
        if size( X,1) ~= size( xMu,1)
            needLearnParams = 1;
            clear  xMu Vimg Dimg
        else
            nX = X - repmat( xMu, 1, N );
            nX =  diag(diag(Dimg).^(-1/2))*Vimg'*nX;
        end
    end

    if needLearnParams  % whiten and reduce
            xMu = mean(X( :, trainThis),2);
            [d N] = size(X); 
            Xc = X - repmat( xMu, 1, N ); 
            [V D ] = pcaSR(Xc(:,trainThis)', imagePercent);
            Dimg = D;
            Vimg = V;
            nX = diag(diag(Dimg).^(-1/2))*Vimg'*Xc;     
            % display( ['num PCA image modes is ' int2str(size(V,2)) ' out of ' ...
            %     int2str( min(size(Xc))) ]);
    end
    
    % save PCA parameters
    if ~isempty( parameterFile )  
        save(  parameterFile,  'xMu', 'Vimg', 'Dimg', 'imageDimensions', '-v7.3' );
    end
    
    % If unit length PCA
    if imgFeatureMode ==5
        nX = normc(nX);
    end
    
elseif imgFeatureMode == 6 %Laplacian Eigenmaps
    kNN= 5;
    nDim = 10;
    sigForEmb = 2.5;
    % (rows correspond to observations, columns to dimensions)
    [nX, mapping] = compute_mapping(X', 'Laplacian', nDim, kNN, sigForEmb);  
     nX = nX';
    % visualizeEmbeddedImages( nX(1:3,:),  [] ); %'embeddingMovie.avi'
    % pause;

elseif imgFeatureMode == 7 % Pixel intensity
    nX =  X;  
    clear X
    
% note that imgFeatureMode == 8 is PCA    
 
elseif imgFeatureMode == 9 % shifted to center, unit length pixel
    nX =  X -  repmat( mean( X,1 ), size(X,1), 1 );
    nX =  normc( nX );
        
elseif imgFeatureMode == 10 % Haar features
    haarScale = 2;
    desiredNumHaar = 10000;
    nX = calcHaarFeatures( reshape(X,s), haarScale, desiredNumHaar, []);
    %nX = normc( nX);
    
elseif imgFeatureMode == 11 % Gabor filter bank
    gaborScale = 2;
    desiredNumGabor = 5000;
    nX = calcGaborFeatures( I, gaborScale, desiredNumGabor);
    
elseif imgFeatureMode == 12 || imgFeatureMode == 13,
    parms.size=5;
    parms.coRelWindowRadius=10;
    parms.numRadiiIntervals=4;
    parms.numThetaIntervals=8;
    parms.varNoise=25*3*36;
    parms.autoVarRadius=1;
    parms.saliencyThresh=0; % I usually disable saliency checking
    parms.nChannels=1;
    radius=(parms.size-1)/2; % the radius of the patch
%     marg=radius+parms.coRelWindowRadius;
    nX=zeros(parms.numRadiiIntervals*parms.numThetaIntervals,s(3));
    for k1=1:s(3),
        Xtemp=double(reshape(X(:,k1),s(1),s(2)));
        x=round(s(1)/2);
        y=round(s(2)/2);
        [resp,drawCoords,salientCoords,uniformCoords]=ssimDescriptor(Xtemp ,parms ,x,y);
        nX(:,k1)=resp;
    end   
   
    % combines SSIM, Haar, LBP
    if imgFeatureMode == 13  
        
        % Haar features
        haarScale = 2;
        desiredNumHaar = 750;  % number Haar, number LBP (same)
        nXHaar = calcHaarFeatures( reshape(X,s), haarScale, desiredNumHaar, []);
        
        %LBP code image using sampling points in SP, no mapping
        SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
        temp = lbp( reshape( X(:,1),imageDimensions ),SP, 0,'i'); 
        X2 = zeros( size(temp(:),1), N);
        for i2 = 1:N
            temp = lbp( reshape( X(:,i2), imageDimensions ),SP,0,'i');
            X2(:,i2) = temp(:);     
        end    
        
        nXLBP = X2; clear X2 temp
        if size( nXLBP,1) > desiredNumHaar
            downSampleAmount = round( size(nXLBP,1)/desiredNumHaar );         
            nXLBP = nXLBP(1:downSampleAmount:end, :); %X2;  %normc(X2);
        end
               
        % SSIM, Haar, LBP
        nX = [nX; nXHaar; nXLBP ];        
    end
    
    
else    % default case    
    %nX = normc( nX);
    error( 'SR: inappropriate image Feature mode '); 
end