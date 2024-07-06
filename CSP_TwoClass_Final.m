eeglab;


close all;
clear;

participantNum = '18';

cd 'C:/Users/Rover/Documents/CCS/FYRP/all_sub_data/EEG_LAB'

%subject_all_data = importdata("S4/all_data_EEG.mat");

subject_data = load('-mat', strcat('C:/Users/Rover/Documents/CCS/FYRP/all_sub_data/EEG_LAB/Subj_data/S', participantNum, '/all_EEG.mat'));

subject_EEG = subject_data.EEG_dataset;

channelLocs = importdata("chanlocs_EEG.mat");


% BELOW THIS IS THE EEGLAB SCRIPT:
EEG.etc.eeglabvers = '2023.0'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_importdata('dataformat','array','nbchan',0,'data', 'subject_EEG','setname','EEG_Test_ForScript','srate',100,'pnts',0,'xmin',0,'chanlocs','channelLocs');
EEG = pop_editset(EEG, 'srate', [100], 'subject', participantNum);
EEG = pop_epoch( EEG, {  'TLE'  }, [0  12], 'newname', 'EEG_Test_ForScript_epochs', 'epochinfo', 'yes');
EEG = pop_reref( EEG, []);
EEG = pop_importevent( EEG, 'append','no','event', strcat('C:\\Users\\Rover\\Documents\\CCS\\FYRP\\all_sub_data\\EEG_LAB\\sub_trial_info\\TrialInfo_sub', participantNum ,'.txt'),'fields',{'epoch','type','latency'},'skipline',1,'timeunit',0.001,'align',0);

% ERDS plot (Activate if needed):
%figure; pop_newtimef( EEG, 1, 46, [0  11990], [3         0.8] , 'topovec', 53, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'CPz', 'baseline',[0], 'plotphase', 'off', 'padratio', 1);

% Keep 'EEG' unchanged. From this point, use MRCP_EEG
MRCP_EEG = EEG;

nFullTrials = MRCP_EEG.trials;
lo_cutoffFreq = 12.5;
hi_cutoffFreq = 30;
sampleRate = 100;
filterOrder = 4; % Sharpness of cutoff.

% Uses 'butterBandpass' function to apply a zero-phase order shift
% butterworth bandpass filter. between 12.5-30Hz
for i = 1:nFullTrials
    seg = MRCP_EEG.data(:,:,i);
    seg = butterBandpass(seg, lo_cutoffFreq, hi_cutoffFreq, sampleRate, filterOrder, 64);
    MRCP_EEG.data(:,:,i) = seg;
end

% Use this to check if the bandpass is working:
% for i = 1:nFullTrials
%     [pxx,f] = pwelch(MRCP_EEG.data(1,:,i),1200,200, 200, sampleRate);
%     fq_power = 10*log10(pxx);
%     plot(f,fq_power)
%     pause(1)
% end

% Create a test and train dataset using CVPartition. (10-fold cross validation)

rng('default'); % Create the same random subset each time.
cvp = cvpartition(length(MRCP_EEG.data(1,1,:)), 'kFold',10);

% To access the fold: 

%trainingData = MRCP_EEG.data(:,:,cvp.training(1));
%testData = MRCP_EEG.data(:,:,cvp.test(1));

nChannels = 64;

% Initialize the figures
f1 = figure();
sgtitle('Feature Space')
f1.Position = [1, 1, 1600, 900];
f2 = figure();
sgtitle('Confusion Matrix')
f2.Position = [1, 1, 1600, 900];

% Go through each fold of the 10-folds.
for fold = 1:10
    % Determine the training and test data.
    training_Data = MRCP_EEG.data(:,:,cvp.training(fold));
    test_Data = MRCP_EEG.data(:,:,cvp.test(fold));

    % The amount of trials depends on the fold, some differ.
    nTrainingTrials = sum(cvp.training(fold));
    nTestTrials = sum(cvp.test(fold));

    % From the training data, select the target and reference windows. 
    % Target = 5000-8000ms and Ref = 0-2000ms
    EEG_target = training_Data(:,500:800,:); % Take movement initiation as target (500-800)
    EEG_reference = training_Data(:,1:201,:); % Take fixation period as reference 1 - 201 (1-201)

    % Initiate the matrices
    covmat_reference = [zeros(nChannels)];
    covmat_target = [zeros(nChannels)];

    % Calculate the covariance matrix for the two conditions
    for i = 1:nTrainingTrials
        seg = EEG_reference(:,:,i);
        seg = seg-mean(seg,2);
        covmat_reference = covmat_reference + seg*seg'/(size(seg,2)-1);

        seg = EEG_target(:,:,i);
        seg = seg-mean(seg,2);
        covmat_target = covmat_target + seg*seg'/(size(seg,2)-1);
    end

    % Average the covariance matrix for the amount of trials in the fold.
    covmat_target = covmat_target/nTrainingTrials;
    covmat_reference = covmat_reference/nTrainingTrials;
    
    % Compute the eigen vectors to apply to the data.
    [V, D] = eig(covmat_target, covmat_reference);
    [eigenValues, sidx] = sort(diag(D), 'descend');
    eigVec = V(:,sidx);

    % Apply the CSP spatial filter to the Test and Training set:

    % TRAINING
    for i = 1:nTrainingTrials
        % Gets one trial
        seg = training_Data(:,:,i);
        % Apply CSP spatial filter
        seg = eigVec(:,:)' * seg;
        % Put it back into training_Data
        training_Data(:,:,i) = seg;
    end

    % TEST
    for i = 1:nTestTrials
        seg = test_Data(:,:,i);
        seg = eigVec(:,:)' * seg;
        test_Data(:,:,i) = seg;
    end

    % Calculate logVariance for Training & Test set

    % TRAINING

    % Initiate feature df
    interval = 200; % 200 samples = 2000ms = 2sec window that the logVariance is calculated for.
    TRAININGfeatures_logVariance = [zeros(64, nTrainingTrials, (1200/interval))];
    timeWindowInterval = interval - 1;


    % TW_int stands for timewindow_int 
    % Calculate the logVariance for the training Dataset
    for trials = 1:nTrainingTrials
        TW_int = 1;
        for timeWindow = 1:interval:1200
            TRAININGfeatures_logVariance(:,trials,TW_int) = logVarianceFeature(training_Data(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
            TW_int = TW_int + 1;
        end
    end

    % TEST

    % NExt part calculates the logVariance for Test Dataset
    % Initiate feature df,
    interval = 200;
    TESTfeatures_logVariance = [zeros(64, nTestTrials, (1200/interval))];
    timeWindowInterval = interval - 1;


    for trials = 1:nTestTrials
        TW_int = 1;
        for timeWindow = 1:interval:1200
            TESTfeatures_logVariance(:,trials,TW_int) = logVarianceFeature(test_Data(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
            TW_int = TW_int + 1;
        end
    end
    
    % features_logVariance: (channel,trials,time_window)
    %(:,:,1) = 0-2000ms REFERENCE, REST
    %(:,:,2) = 2000-4000ms
    %(:,:,3) = 4000-6000ms
    %(:,:,4) = 6000-8000ms MOVEMENT
    %(:,:,5) = 8000-10000ms
    %(:,:,6) = 10000-12000ms

    % Next part sets up the dataset for training (inc labels), so that an optimal
    % hyperplane can be established.
    trainingTargetFeatures = [];
    trainingReferenceFeatures = [];
    trainingFeatures = [];
    
    %Target
    trainingTargetFeatures(:,1) = TRAININGfeatures_logVariance(1,:,4); % First component, time window 4
    trainingTargetFeatures(:,2) = TRAININGfeatures_logVariance(64,:,4); % Last component
    
    %Ref
    trainingReferenceFeatures(:,1) = TRAININGfeatures_logVariance(1,:,1); % First Component, time window 1
    trainingReferenceFeatures(:,2) = TRAININGfeatures_logVariance(64,:,1); % Last component

    % Some manipulation of datasets to create a training dataset, with
    % labels. (:,1) = component1, (:,2) = component2 (:,3) = label 1/0
    % Ends up with trainingFeatures, which contains both condition's
    % features and labels.
    trainingFeatures(:,:) = trainingTargetFeatures(:,:);
    %set labels for target condition
    trainingFeatures(1:length(trainingFeatures),3) = 1;
    trainingFeatures = trainingFeatures';
    storeRefFeatures = trainingReferenceFeatures';
    % Set labels for ref condition
    storeRefFeatures(3,:) = 0;
    % Add it to 'trainingFeatures'
    trainingFeatures = [trainingFeatures, storeRefFeatures];
    trainingFeatures = trainingFeatures';

    % Next part calculates the optimal hyperplane for the training data.
    md1 = fitcdiscr(trainingFeatures(:,1:2), trainingFeatures(:,3), 'discrimtype', 'linear');
    K = md1.Coeffs(1,2).Const;
    L = md1.Coeffs(1,2).Linear;

    f = @(x1,x2) K + L(1)*x1 + L(2)*x2; % Determine the decision boundary

    % PLotting:
    figure(f1);
    subplot(5,2,fold) % Creates a 5x2 figure for 10 plots. One for each fold.
    title(strcat('Fold', num2str(fold)));
    hold on;
    scatter(trainingTargetFeatures(:,1), trainingTargetFeatures(:,2), "red", "x");
    scatter(trainingReferenceFeatures(:,1), trainingReferenceFeatures(:,2), "blue", "x");
    %scatter(TRAININGfeatures_logVariance(1,:,4), TRAININGfeatures_logVariance(64,:,4), "red", 'x');
    %scatter(TRAININGfeatures_logVariance(1,:,1),TRAININGfeatures_logVariance(64,:,1), "blue", 'x');
    xlabel('First Component');
    ylabel('Last Component');
    scatter(TESTfeatures_logVariance(1,:,1),TESTfeatures_logVariance(64,:,1), "black", 'o', 'filled');
    scatter(TESTfeatures_logVariance(1,:,4),TESTfeatures_logVariance(64,:,4), "black",'o');
    fimplicit(f) % Plot decision boundary
    hold off;

    testReferenceFeatures = [];

    testReferenceFeatures(:,1) = TESTfeatures_logVariance(1,:,1); % First component
    testReferenceFeatures(:,2) = TESTfeatures_logVariance(64,:,1); % Last component
    testReferenceFeatures(:,3) = 0;

    testTargetFeatures = [];

    testTargetFeatures(:,1) = TESTfeatures_logVariance(1,:,4); % First component
    testTargetFeatures(:,2) = TESTfeatures_logVariance(64,:,4); % Last component
    testTargetFeatures(:,3) = 1;

    testFeatures = [testReferenceFeatures; testTargetFeatures];

    % Determine the accuracy of each fold's test data.
    % label_predicition predicts what class the feature data coresponds to.
    % in the first column. 
    label_prediction = predict(md1, testFeatures(:,1:2));
    % add the acutal labels in the 2nd column.
    label_prediction(:,2) = testFeatures(:,3);
   
    % Create a confusion matrix
    figure(f2);
    subplot(5,2,fold)
    confusionchart(label_prediction(:,2), label_prediction(:,1));
    title(strcat('Fold', num2str(fold)));
end

% Create legend after producing all the graphs for each fold
figure(f1);
legend("Movement Condition", "Rest Condition", "Test Data - Rest", "Test Data - Movement", "Optimal Hyperplane between two classes", 'Position', [0.43 0.05 0.08 0.1])

% Save graphs into Results_TwoClass folder
saveas(f1, strcat('Results_BinaryClassifier\featureGraphs_', participantNum , '.png'))
saveas(f2, strcat('Results_BinaryClassifier\ConfMats_', participantNum , '.png'))






% BELOW HERE IS WHERE I WROTE THE SCRIPT, More of a workspace, DO NOT USE. 


% 
% % Take 5-7sec window as the tagert window (movement initiated class)
% EEG_target = pop_epoch(MRCP_EEG, {}, [5 7]);
% 
% % As reference signal, the 2-4 second window is taken (reference class).
% % Maybe 3-5 sec should be used.
% EEG_reference = pop_epoch(MRCP_EEG, {}, [2 4]);
% 
% % Then a covariance matrix is established for the 2 classes.
% covmat_reference = [zeros(nChannels)];
% 
% covmat_target = [zeros(nChannels)];
% 
% 
% 
% % Calculate the covariance matrix for the reference signal and the target
% % signal. Reference is between 2-4seconds and target 5-7 seconds.
% for i = 1:nFullTrials
%     seg = EEG_reference.data(:,:,i);
%     seg = seg-mean(seg,2);
%     covmat_reference = covmat_reference + seg*seg'/(size(seg,2)-1);
% 
%     seg = EEG_target.data(:,:,i);
%     seg = seg-mean(seg,2);
%     covmat_target = covmat_target + seg*seg'/(size(seg,2)-1);
% end
% 
% % Average the covariance matrix
% covmat_target = covmat_target/nFullTrials;
% covmat_reference = covmat_reference/nFullTrials;
% 
% % Determine the eigen vectors of the classes
% [V, D] = eig(covmat_target, covmat_reference);
% [eigenValues, sidx] = sort(diag(D), 'descend');
% eigVec = V(:,sidx);
% 
% 
% CSP_EEG = MRCP_EEG;
% 
% % Apply the eigenvectors to the data (CSP spatial filter)
% for i = 1:nFullTrials
%     for section = 1:200:1001
%         % Gets one trial
%         seg = CSP_EEG.data(:,section:section+199,i);
%         % Apply CSP spatial filter
%         seg = eigVec(:,:)' * seg;
%         CSP_EEG.data(:,section:section+199,i) = seg;
%         %disp(seg)
%     end
% end
% 
% % Calculates the logVariance of the full trial.
% interval = 200; % 200 samples = 2000ms = 2sec
% tempVar = 1;
% 
% features_logVariance = [zeros(64, nFullTrials, (1200/interval))];
% 
% % To correct the values according to the sizes of the matrices
% timeWindowInterval = interval - 1;
% 
% % Uses 'LogVarianceFeature.m' to calculate the logVariance of the data
% for trials = 1:nFullTrials
%     tempVar = 1;
%     for timeWindow = 1:interval:1200
%         features_logVariance(:,trials,tempVar) = logVarianceFeature(CSP_EEG.data(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
%         %disp(tempVar)
%         tempVar = tempVar + 1;
%     end
% end
% 
% % features_logVariance: (channel,trials,time_window)
% %(:,:,1) = 0-200ms
% %(:,:,2) = 200-400ms
% %(:,:,3) = 400-600ms
% %(:,:,4) = 600-800ms
% %(:,:,5) = 800-1000ms
% %(:,:,6) = 1000-1200ms
% 
% % THe CSP is trained to increase the variance between 200-400ms (ref) and
% % 500-700ms (target)
% 
% % Each window is 200ms in length. So, one axis is showing the logVariance
% % during the movement and the other axis the logVariance activity in the
% % 200ms window before. 
% figure(100)
% scatter(features_logVariance(1,:,3), features_logVariance(1,:,4));
% hold on;
% scatter(features_logVariance(64,:,3),features_logVariance(64,:,4));
% xlabel('LogVariance of the CSP filtered window T-1');
% ylabel("LogVariance of the CSP filtered EEG signal window T");
% 
% % However, I think this is the usual way of doing it (THIS WAY MAKES SENSE):
% figure(101)
% scatter(features_logVariance(1,:,3), features_logVariance(64,:,3));
% hold on;
% scatter(features_logVariance(1,:,4),features_logVariance(64,:,4));
% xlabel('First Component');
% ylabel("Last Component");

% The following part is now part of preprocess_raw_data 
%
% % Combining all the files containing the list of trials for each condition.
% % In combine_cond; the row indicates the condition. Below is the guide:
% % row 1 =  LC power                 row 7 = SC_power
% % row 2 = LC precision 2            row 8 = SC_precision2
% % row 3 = LC precision 5            row 9 = SC_precision5
% % row 4 = LS_Power                  row 10 = SS_power
% % row 5 = LS_precision2             row 11 = SS_precision2
% % row 6 = LS_precision5             row 12 = SS_precision5
% combine_cond = [subject_data.indexes_LC_power; subject_data.indexes_LC_precision2; 
%     subject_data.indexes_LC_precision5; subject_data.indexes_LS_power;
%     subject_data.indexes_LS_precision2; subject_data.indexes_LS_precision5;
%     subject_data.indexes_SC_power; subject_data.indexes_SC_precision2;
%     subject_data.indexes_SC_precision5; subject_data.indexes_SS_power;
%     subject_data.indexes_SS_precision2; subject_data.indexes_SS_precision5];
% 
% 
% cond_list(1,3) = "latency"; 
% cond_list(1,2) = "type";
% cond_list(1,1) = "epoch";
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(1,i);
%     cond_list(tempCondPos + 1, 2) = "LC_Power";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(2,i);
%     cond_list(tempCondPos + 1, 2) = "LC_Precision_2";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(3,i);
%     cond_list(tempCondPos + 1, 2) = "LC_Precision_5";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(4,i);
%     cond_list(tempCondPos + 1, 2) = "LS_Power";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(5,i);
%     cond_list(tempCondPos + 1, 2) = "LS_Precision_2";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(6,i);
%     cond_list(tempCondPos + 1, 2) = "LS_Precision_5";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(7,i);
%     cond_list(tempCondPos + 1, 2) = "SC_Power";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(8,i);
%     cond_list(tempCondPos + 1, 2) = "SC_Precision_2";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(9,i);
%     cond_list(tempCondPos + 1, 2) = "SC_Precision_5";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(10,i);
%     cond_list(tempCondPos + 1, 2) = "SS_Power";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(11,i);
%     cond_list(tempCondPos + 1, 2) = "SS_Precision_2";
% end
% 
% for i = 1:length(combine_cond)
%     tempCondPos = combine_cond(12,i);
%     cond_list(tempCondPos + 1, 2) = "SS_Precision_5";
% end
% 
% % Create a new column for the placement of the events, start from row 2 to keep
% % the column names.
% p = 1;
% for i = 1:length(cond_list) - 1
%      cond_list(i + 1,1) = p;
%      p = p + 1;
% end
% 
% % include latencies of each trial for event import:
% p = 0;
% 
% for i = 1:length(cond_list) - 1
%     cond_list(i + 1, 3) = p;
%     p = p + 12000;
% end
% 
% 
% % Save to txt file, so that it can be accessed by EEGlab
% %writematrix(cond_list, 'TrialInfo.txt');

