% Samuel Rivera
% Feb 4, 2011
% Notes: this function just takes the KSDA model (a struct
%  containing the trained model), and the testImages.  It returns the
% classification.  Much of this code was taken from the function
% 'KSDA_MaxHomo.m' which was written by Di You
%
% testImages in ( HxWxN), for N samples
% test_label is 1xN for integer labels
% kSDAModel returned by srTrainKSDA

function [ classEstimates rate ] = srTestKSDA( kSDAModel, testImages, test_label )


v = kSDAModel.v;
trainingdata = kSDAModel.trainingdata ;
op_sigma = kSDAModel.sigma;
nc = kSDAModel.nc;
imgFeatureMode = kSDAModel.imgFeatureMode;
oriImageSize= kSDAModel.oriImageSize;
targetSize = kSDAModel.targetSize;
maxNumFeatures= kSDAModel.maxNumFeatures;

clear kSDAModel

% Samuel Rivera added this for feature extraction

% make sure images are the size expected by the model
testImageSize = size(  testImages(:,:,1));    
if testImageSize(:) ~= targetSize(:)
    error( 'SR: the test images are not the correct size.  Perhaps you are using incorrect model.');
end

imageParamFile = 'tempKSDAImagePCAParameters.mat';
doItFast = 1;
trainThis = 1:size(testImages,3);  
[ testingdata ] = formatImages( imageParamFile, testImages, ...
                       doItFast, imgFeatureMode, trainThis, []) ;
clear testImages

% reduce dimensionality if necessary
if size(testingdata,1) > maxNumFeatures
    downSampleAmount = round( size(testingdata,1)/maxNumFeatures );
    testingdata = testingdata(1:downSampleAmount:end, :); 
end

%----------------- prep stage --------------------------

C = length(nc);

% % %%% get pairwise distance matrix
% % l=size(trainingdata,2);
% % A = trainingdata'*trainingdata;
% % dA = diag(A);
% % DD = repmat(dA,1,l) + repmat(dA',l,1) - 2*A;
% % 
% % % Gram matrix
% % K1=exp(-DD/(2*op_sigma^2));

params.sigma = op_sigma;
[K1 K2 ] = calcGramAndKappa( trainingdata, testingdata, 1, params );

  % ---------------------- testing stage -------------------------------
 train=v'*K1;
 
% %  nXtest=size(testingdata,2);
% % for i=1:nXtest
% %     B=trainingdata-repmat(testingdata(:,i),1,l);
% %     B=B.^2;
% %     dd(i,:)=sum(B,1);
% % end
% % dd=dd';
% %  K2=exp(-dd/(2*op_sigma^2));

test=v'*K2;


%%% nearest mean classifier

%         rate=Nmean(train', test', H, C, NH, test_label) 

%%% nearest neighbor classifier

[rate classEstimates ]=NearestNeighbor(train',test',test_label,C,nc);
