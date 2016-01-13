function [ rate testEtimate ] = NearestNeighbor(train, test, correct_class, category, nc)
% function: classify the testing data according to nearest neighbor rule
%
%input:
%train: training data
%test:  testing data
%correct_class: the true class label of the testing data
%category:  the number of classes 
%nc: category-by-1 matrix, indicationg the number of samples in each class
%    for training data
%
%output:
% rate: recognition rate
% testEtimate: class labels

n_test = size(test,1);
n_train = size(train,1);

testEtimate = zeros( 1, n_test);

rec = 0;
for i = 1:n_test 
    
    temp = repmat(test(i,:),n_train,1)-train;
    temp = temp.^2;
    dist = sum(temp',1);
                
    %[val,I]=min(dist);
    [val,I] = sort(dist);
    val = val(1);
    I = I(1);
    for class = 1:category
        if (I<=sum(nc(1:class)))
            found = class;
            break;
        end
    end
                
    if (found == correct_class(i)) 
        rec=rec+1;    
    end
    
    testEtimate(i) = found; % found is predicted class label
end
rate = rec/n_test;