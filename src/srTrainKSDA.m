% Samuel Rivera
% Feb 4, 2011
% Notes: this function takes code from 'KSDA_MaxHomo.m' written by Di You
% in order to train the KSDA.  I will store those parameters, and use them
%       later

function [ kSDAModel ] = srTrainKSDA( trainingImages, nc, imgFeatureMode, targetSize )


% This function is an implementation of Kernel Subclass Discriminant
% Analysis (KSDA) with parameters optimized by the Homoscedastic criterion
% proposed in the paper "Kernel Optimization in Discriminant Analysis".

% Input: 
%
% trainingImages, (h x w x N) black and white training image matrix
% Note that this matrix should be formatted as follows (suppose n_i is the number of 
% samples in ith class): The first n_1 columns are the n_1 samples from 1st
% class, the next n_2 columns are the n_2 samples from 2nd class, etc.
% 
% nc: 1-by-C vector, indicationg the number of samples in each class
% for the training data.  It is assumed to be ordered correctly
% 
% imgFeatureMode: integer specifying feature extraction (1, 3, 7,or 10 recommended):
%     1: unit pixel
%     2: C1 features
%     3: for LBP ( local binary pattern)
%      4 for Berkeley Wavelet Transform  (DO NOT USE)
%      5 for unit PCA
%      6 for Laplacian Eigenmaps ( DO NOT USE)
%      7 for Pixel
%      8 for PCA
%      9 for center shifted unit pixel
%     10 for Haar features
%     11 for Gabor features (not implemented)
%     12 for self similarity features
% 
% legacy code needed these parameters:
% trainingdata: The p-by-N training data matrix, where N is the number of 
% training samples and p is the number of dimensions of the data. 
% C: number of classes (now set as length(nc))


% Output:

% kSDAModel: a struct with
%     kSDAModel.v = v;
%     kSDAModel.trainingdata = trainingdata;
%     kSDAModel.sigma = op_sigma;
%     kSDAModel.nc = nc;
%     kSDAModel.imgFeatureMode = imgFeatureMode;
%     kSDAModel.oriImageSize = oriImageSize;
%     kSDAModel.targetSize= targetSize;
% 
% Copyrighted code
% (c) Di You, Onur Hamsici and Aleix M Martinez
% with additions by Samuel Rivera in Feb 2011
%
% For additional information contact the authors


% Samuel Rivera added this for  resizing, then feature extraction
targetSize = round( targetSize);
numSamples = size(trainingImages,3);

% resize if not the right size
if sum(size(trainingImages(:,:,1)) == targetSize(:)' ) ~= 2
    trainingImagesSmall = zeros(  targetSize(1),  targetSize(2), numSamples);
    % save smaller version of image
    for i1 = 1:numSamples
        trainingImagesSmall(:,:,i1) = imresize( trainingImages(:,:,i1), targetSize(:)' );
    end
else
    trainingImagesSmall = trainingImages;
end
oriImageSize = size( trainingImages(:,:,1));
clear trainingImages
        
% feature exractions
imageParamFile = 'tempKSDAImagePCAParameters.mat';
doItFast = 1;
trainThis = 1:size(trainingImagesSmall,3);
[ trainingdata ] = formatImages( imageParamFile, trainingImagesSmall, ...
                       doItFast, imgFeatureMode, trainThis, []) ;
C =  length(nc);


% reduce dimensionality if necessary
maxNumFeatures = 750;
if size(trainingdata,1) > maxNumFeatures
    downSampleAmount = round( size(trainingdata,1)/maxNumFeatures );
    trainingdata = trainingdata(1:downSampleAmount:end, :); 
end

%-----------------------------------
%make local texture for kde, resize if not the right size
numTexPdf = nc(1); %first one is the positive class
targetSizeTexPdf = [10;10];
if sum(size(trainingImagesSmall(:,:,1)) == targetSizeTexPdf(:)' ) ~= 2  
    % save smaller version of image
    trainingImagesTexPdf = zeros(  targetSizeTexPdf(1),  targetSizeTexPdf(2), numTexPdf);    
    for i1 = 1:numTexPdf
        trainingImagesTexPdf(:,:,i1) = imresize( trainingImagesSmall(:,:,i1), targetSizeTexPdf(:)' );
    end
else
    trainingImagesTexPdf = trainingImagesSmall(:,:,1:numTexPdf);
end
clear trainingImagesSmall

imageParamFile = 'tempKSDAImagePCAParameters.mat';
doItFast = 1;
imgFeatureModeTexPdf = 9;
trainThis = 1:size(trainingImagesTexPdf,3);
[ trainingdataTexPdf ] = formatImages( imageParamFile, trainingImagesTexPdf, ...
                       doItFast, imgFeatureModeTexPdf, trainThis, []) ;    
%---------------------
                   
%%% Nearest Neighbor clustering of the data

Ytrain = NNclassclustering(trainingdata',C,nc);
trainingdata=Ytrain';
l=size(trainingdata,2);

%%% get pairwise distance matrix

A = trainingdata'*trainingdata;
dA = diag(A);
DD = repmat(dA,1,l) + repmat(dA',l,1) - 2*A;

s1=sum(sum(DD,1));
num=l*(l-1)/2;
mean_DD=s1/2/num;
options = optimset('LargeScale','off', 'GradObj','off',...
    'HessUpdate','bfgs', 'TolX',1e-10, 'MaxFunEvals',5000, 'MaxIter',10000);


%%% optimize the kernel parameter and number of subclasses in KSDA

for ii=1:1  % maximum number subclasses
    H = ii*ones(1,C);
    NH = get_NH(C,H,nc);
    X0=sqrt(mean_DD/2);
    [Sigma(ii),fval(ii)] = fminunc(@(sigma)Maxhomo(H, C, NH, l, sigma,DD),X0,options);
end
[F,ind]=min(fval);
op_H=ind;
op_sigma=Sigma(ind);

%%% learn KSDA subspace after selecting optimal parameters
H = op_H*ones(1,C);
NH = get_NH(C,H,nc);
K1=exp(-DD/(2*op_sigma^2));
v=KSDA(C,trainingdata,H,NH,K1);
 
 
% store model for future use 
kSDAModel.v = v;
kSDAModel.trainingdata = trainingdata;
kSDAModel.sigma = op_sigma;
kSDAModel.nc = nc;
kSDAModel.imgFeatureMode = imgFeatureMode;
kSDAModel.oriImageSize = oriImageSize;
kSDAModel.targetSize= targetSize;
kSDAModel.maxNumFeatures = maxNumFeatures;
kSDAModel.type = 1; %1 for kLDA

kSDAModel.trainingdataTexPdf = trainingdataTexPdf;
kSDAModel.imgFeatureModeTexPdf = imgFeatureModeTexPdf;
kSDAModel.targetSizeTexturePdf = targetSizeTexPdf;