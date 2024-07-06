function W = trainLDA(features_dataset)
% The features_dataset has the following dimensions 3 x nTrials(1168). 
% Column 1: component 1, column2: component 2, column3: label 1/0
    
    % Gathers the amount of classes (2)
    nClasses = numel(unique(features_dataset(:,3)));

    if nClasses == 2
        disp("Good");

        class1 = features_dataset(features_dataset(:,3) == 1, 1:2);
        class0 = features_dataset(features_dataset(:,3) == 0, 1:2);

        prior1 = length(class1)/(length(class1)+length(class0));
        prior0 = length(class0)/(length(class1)+length(class0));

        meanClass1 = mean(class1); % gets mean y and x value
        meanClass0 = mean(class0); % gets mean y and x value 

        class1_c = class1-meanClass1;
        class0_c = class0-meanClass0;

        cov1 = cov(class1_c);
        cov0 = cov(class0_c);

        W = (meanClass0 - meanClass1)
        
    elseif (nClasses <= 1)
        disp("ERROR: There are too little classes in column 3.");
    else
        disp("ERROR: There are too many classes in column 3")
    end
end