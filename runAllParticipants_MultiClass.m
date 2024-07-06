eeglab;


close all;
clear;

% Participant 5 is missing data
allParticipants = [3,4,6,7,8,9,10,11,12,13,14,15,16,17,18];

%allParticipants = [3,4,6];

save_all_results_knn = [];
save_all_results_lda = [];

for participants = 1:length(allParticipants)
    participantNum = num2str(allParticipants(participants)); % Selects the correct participant from the list on line 14 and changes it into a string.

    cd 'C:/Users/Rover/Documents/CCS/FYRP/all_sub_data/EEG_LAB'

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

    % Create a test and train dataset using CVPartition. (10-fold cross validation)

    rng('default'); % Create the same random subset each time.
    cvp = cvpartition(length(MRCP_EEG.data(1,1,:)), 'kFold',10);

    % To access the fold:

    %trainingData = MRCP_EEG.data(:,:,cvp.training(1));
    %testData = MRCP_EEG.data(:,:,cvp.test(1));

    nChannels = 64;


    % Creates a new dataset that splits conditions into three: Power, Prec5,
    % Prec2. Had to use this instead of pop_eeg from matlab as managing the
    % struct in the fold did not work.

    % 1 = power
    % 2 = prec5
    % 3 = prec2

    condition = zeros(1, length(MRCP_EEG.urevent));

    for index = 1:length(MRCP_EEG.data(1,1,:))
        if (strcmp(MRCP_EEG.urevent(index).type, 'LC_Power') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'SC_Power') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'LS_Power') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'SS_Power') == 1)
            condition(1, index) = 1;
        elseif (strcmp(MRCP_EEG.urevent(index).type, 'LC_Precision_5') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'SC_Precision_5') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'LS_Precision_5') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'SS_Precision_5') == 1)
            condition(1, index) = 2;
        elseif (strcmp(MRCP_EEG.urevent(index).type, 'LC_Precision_2') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'SC_Precision_2') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'LS_Precision_2') == 1) || (strcmp(MRCP_EEG.urevent(index).type, 'SS_Precision_2') == 1)
            condition(1, index) = 3;
        else
            condition(1, index) = 4;
        end
    end

    % I have to setup a dataset that contains information of the condition of
    % each trial. So that I do not have to use pop_EEG function, which requires
    % a EEGLab dataset, which 'training_Data' are not. With that dataset,
    % training_Data can be split up so that the training data can be setup.

    % % Making it so that channel 65 in the MRCP_EEG contains the condition
    % information. This is to make sure that with the scrabbling of trials in the
    % k-folds, this information is not lost.

    % To access the trial condition in data_with_cond, see channel 65
    % (65,1,trial). Fill in trial for which trial.

    data_with_cond = MRCP_EEG.data;

    for trial = 1:length(data_with_cond(1,1,:))
        data_with_cond(65,1,trial) = condition(1, trial);
    end

    % Initialize the figures
    f1 = figure();
    sgtitle('Feature Space')
    f1.Position = [1, 1, 1600, 900];

    % Go through each fold of the k-fold.
    for fold = 1:10

        EEG_power = zeros(65,1200,length(data_with_cond(1,1,:)));
        EEG_prec5 = zeros(65,1200,length(data_with_cond(1,1,:)));
        EEG_prec2 = zeros(65,1200,length(data_with_cond(1,1,:)));
        EEG_wo_power = zeros(65,1200,length(data_with_cond(1,1,:)));
        EEG_wo_prec5 = zeros(65,1200,length(data_with_cond(1,1,:)));
        EEG_wo_prec2 = zeros(65,1200,length(data_with_cond(1,1,:)));


        % Determine the training and test data.
        training_Data = data_with_cond(1:64,:,cvp.training(fold));
        test_Data = data_with_cond(1:64,:,cvp.test(fold));

        % The amount of trials depends on the fold, some differ in num of trials.
        nTrainingTrials = sum(cvp.training(fold));
        nTestTrials = sum(cvp.test(fold));


        % Setting up the datasets for training; EEG_power, EEG_prec5 &
        % EEG_prec2 (target). The 'wo' (e.g. EEG_wo_prec5) is the reference data.
        for trial = 1:length(data_with_cond(1,1,:))
            if data_with_cond(65,1,trial) == 1
                EEG_power(:,:,trial) = data_with_cond(:,:,trial);

                EEG_wo_prec5(:,:,trial) = data_with_cond(:,:,trial);
                EEG_wo_prec2(:,:,trial) = data_with_cond(:,:,trial);
                disp("1 - power")

            elseif data_with_cond(65,1,trial) == 2
                EEG_prec5(:,:,trial) = data_with_cond(:,:,trial);

                EEG_wo_power(:,:,trial) = data_with_cond(:,:,trial);
                EEG_wo_prec2(:,:,trial) = data_with_cond(:,:,trial);
                disp("2 - prec5")

            elseif data_with_cond(65,1,trial) == 3
                EEG_prec2(:,:,trial) = data_with_cond(:,:,trial);

                EEG_wo_prec5(:,:,trial) = data_with_cond(:,:,trial);
                EEG_wo_power(:,:,trial) = data_with_cond(:,:,trial);
                disp("3 - prec2")
            else
                disp("Organising datasets for training did not work.")
            end
        end

        % Clear out the empty trials in the condition datasets and save it to a
        % new dataset, without zeros.

        counter_power = 1;
        counter_prec5 = 1;
        counter_prec2 = 1;
        counter_power_ref = 1;
        counter_prec5_ref = 1;
        counter_prec2_ref = 1;

        for trial = 1:length(data_with_cond(1,1,:))
            if EEG_power(65,1,trial) == 1
                EEG_power_target(:,:,counter_power) = EEG_power(:,:,trial);
                counter_power = counter_power + 1;

            elseif EEG_prec5(65,1,trial) == 2
                EEG_prec5_target(:,:,counter_prec5) = EEG_prec5(:,:,trial);
                counter_prec5 = counter_prec5 + 1;

            elseif EEG_prec2(65,1,trial) == 3
                EEG_prec2_target(:,:,counter_prec2) = EEG_prec2(:,:,trial);
                counter_prec2 = counter_prec2 + 1;

            else
                disp("Something went wrong #1")

            end

            if EEG_wo_power(65,1,trial) ~= 0
                EEG_power_reference(:,:,counter_power_ref) = EEG_wo_power(:,:,trial);
                counter_power_ref = counter_power_ref + 1;
            end

            if EEG_wo_prec5(65,1,trial) ~= 0
                EEG_prec5_reference(:,:,counter_prec5_ref) = EEG_wo_prec5(:,:,trial);
                counter_prec5_ref = counter_prec5_ref + 1;
            end

            if EEG_wo_prec2(65,1,trial) ~= 0
                EEG_prec2_reference(:,:,counter_prec2_ref) = EEG_wo_prec2(:,:,trial);
                counter_prec2_ref = counter_prec2_ref + 1;
            end
        end

        % From the training data, select the movement initiation window.
        % Target = with condition (e.g. only power, only prec5 or only prec2)
        % and Ref = without condition (e.g. Power reference: prec5+prec2,
        % Prec5 reference: power+prec2 & Prec2 reference: power+prec5)
        % Movement was prompted at 5000-8000ms, which is 500-800 in the
        % dataset.
        EEG_target_power = EEG_power_target(1:64,500:800,:); % Take movement initiation as target of power condition
        EEG_reference_power = EEG_power_reference(1:64,500:800,:); % Take movement initiation non-power condition as reference

        EEG_target_prec5 = EEG_prec5_target(1:64,500:800,:);
        EEG_reference_prec5 = EEG_prec5_reference(1:64,500:800,:);

        EEG_target_prec2 = EEG_prec2_target(1:64,500:800,:);
        EEG_reference_prec2 = EEG_prec2_reference(1:64,500:800,:);


        % Initiate the matrices
        covmat_reference_power = [zeros(nChannels)];
        covmat_target_power = [zeros(nChannels)];

        covmat_reference_prec5 = [zeros(nChannels)];
        covmat_target_prec5 = [zeros(nChannels)];

        covmat_reference_prec2 = [zeros(nChannels)];
        covmat_target_prec2 = [zeros(nChannels)];

        % Calculate the covariance matrix for power
        for i = 1:length(EEG_reference_power(1,1,:))
            seg = EEG_reference_power(:,:,i);
            seg = seg-mean(seg,2);
            covmat_reference_power = covmat_reference_power + seg*seg'/(size(seg,2)-1);
        end

        for i = 1:length(EEG_target_power(1,1,:))
            seg = EEG_target_power(:,:,i);
            seg = seg-mean(seg,2);
            covmat_target_power = covmat_target_power + seg*seg'/(size(seg,2)-1);
        end

        % Calculate the covariance matrix for precision 5 (prec5)
        for i = 1:length(EEG_reference_prec5(1,1,:))
            seg = EEG_reference_prec5(:,:,i);
            seg = seg-mean(seg,2);
            covmat_reference_prec5 = covmat_reference_prec5 + seg*seg'/(size(seg,2)-1);
        end

        for i = 1:length(EEG_target_prec5(1,1,:))
            seg = EEG_target_prec5(:,:,i);
            seg = seg-mean(seg,2);
            covmat_target_prec5 = covmat_target_prec5 + seg*seg'/(size(seg,2)-1);
        end

        % Calculate the covariance matrix for precision 2 (prec2)
        for i = 1:length(EEG_reference_prec2(1,1,:))
            seg = EEG_reference_prec2(:,:,i);
            seg = seg-mean(seg,2);
            covmat_reference_prec2 = covmat_reference_prec2 + seg*seg'/(size(seg,2)-1);
        end

        for i = 1:length(EEG_target_prec2(1,1,:))
            seg = EEG_target_prec2(:,:,i);
            seg = seg-mean(seg,2);
            covmat_target_prec2 = covmat_target_prec2 + seg*seg'/(size(seg,2)-1);
        end


        % Average the covariance matrix for the amount of trials in the fold.
        % For each condition, power, prec5 & prec2
        covmat_target_power = covmat_target_power/length(EEG_target_power(1,1,:));
        covmat_reference_power = covmat_reference_power/length(EEG_reference_power(1,1,:));

        covmat_target_prec5 = covmat_target_prec5/length(EEG_target_prec5(1,1,:));
        covmat_reference_prec5 = covmat_reference_prec5/length(EEG_reference_prec5(1,1,:));

        covmat_target_prec2 = covmat_target_prec2/length(EEG_target_prec2(1,1,:));
        covmat_reference_prec2 = covmat_reference_prec2/length(EEG_reference_prec2(1,1,:));

        % Compute the eigen vectors to apply to the data.

        %Power
        [V_power, D_power] = eig(covmat_target_power, covmat_reference_power);
        [eigenValues_power, sidx_power] = sort(diag(D_power), 'descend');
        eigVec_power = V_power(:,sidx_power);

        % Precision 5
        [V_prec5, D_prec5] = eig(covmat_target_prec5, covmat_reference_prec5);
        [eigenValues_prec5, sidx_prec5] = sort(diag(D_prec5), 'descend');
        eigVec_prec5 = V_prec5(:,sidx_prec5);

        % Precision 2
        [V_prec2, D_prec2] = eig(covmat_target_prec2, covmat_reference_prec2);
        [eigenValues_prec2, sidx_prec2] = sort(diag(D_prec2), 'descend');
        eigVec_prec2 = V_prec2(:,sidx_prec2);

        % Apply the power, precision5 and precision2 CSP spatial filter to the
        % Test and Training set. This creates 6 new datasets, 2 for each
        % condition, namely test and training dataset. To be clear: the CSP is
        % being applied to the training and test data; including all condition
        % types.


        % TRAINING - power csp
        for i = 1:nTrainingTrials
            % Gets one trial
            seg = training_Data(:,:,i);
            % Apply CSP spatial filter
            seg = eigVec_power(:,:)' * seg;
            % Put it back into training_Data
            training_Data_power(:,:,i) = seg;
        end

        % TEST - power csp
        for i = 1:nTestTrials
            seg = test_Data(:,:,i);
            seg = eigVec_power(:,:)' * seg;
            test_Data_power(:,:,i) = seg;
        end


        % TRAINING - precision5 csp
        for i = 1:nTrainingTrials
            % Gets one trial
            seg = training_Data(:,:,i);
            % Apply CSP spatial filter
            seg = eigVec_prec5(:,:)' * seg;
            % Put it back into training_Data
            training_Data_precision5(:,:,i) = seg;
        end

        % TEST - precision5 csp
        for i = 1:nTestTrials
            seg = test_Data(:,:,i);
            seg = eigVec_prec5(:,:)' * seg;
            test_Data_precision5(:,:,i) = seg;
        end

        % TRAINING - precision2 csp
        for i = 1:nTrainingTrials
            % Gets one trial
            seg = training_Data(:,:,i);
            % Apply CSP spatial filter
            seg = eigVec_prec2(:,:)' * seg;
            % Put it back into training_Data
            training_Data_precision2(:,:,i) = seg;
        end

        % TEST - precision2 csp
        for i = 1:nTestTrials
            seg = test_Data(:,:,i);
            seg = eigVec_prec2(:,:)' * seg;
            test_Data_precision2(:,:,i) = seg;
        end

        % Calculate logVariance for Training, Test & Reference set

        % TRAINING

        % Initiate feature df
        interval = 200; % 200 samples = 2000ms = 2sec window that the logVariance is calculated for.
        train_feat_logVar_power = [zeros(64, nTrainingTrials, (1200/interval))];
        train_feat_logVar_prec5 = [zeros(64, nTrainingTrials, (1200/interval))];
        train_feat_logVar_prec2 = [zeros(64, nTrainingTrials, (1200/interval))];

        timeWindowInterval = interval - 1;

        %% Training:
        % TW_int stands for timewindow_int
        % Calculate the logVariance for the training Datasets. Uses
        % 'logVarianceFeature.m' function.

        % Power

        for trials = 1:nTrainingTrials
            TW_int = 1;
            for timeWindow = 1:interval:1200
                train_feat_logVar_power(:,trials,TW_int) = logVarianceFeature(training_Data_power(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
                TW_int = TW_int + 1;
            end
        end

        % Precision 5

        for trials = 1:nTrainingTrials
            TW_int = 1;
            for timeWindow = 1:interval:1200
                train_feat_logVar_prec5(:,trials,TW_int) = logVarianceFeature(training_Data_precision5(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
                TW_int = TW_int + 1;
            end
        end


        % Precision 2

        for trials = 1:nTrainingTrials
            TW_int = 1;
            for timeWindow = 1:interval:1200
                train_feat_logVar_prec2(:,trials,TW_int) = logVarianceFeature(training_Data_precision2(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
                TW_int = TW_int + 1;
            end
        end

        %% TEST data

        % NExt part calculates the logVariance for Test Dataset (all trials are
        % put through the filter for each condition.
        % Initiate feature df,
        interval = 200;
        test_feat_logVar_power = [zeros(64, nTestTrials, (1200/interval))];
        test_feat_logVar_prec5 = [zeros(64, nTestTrials, (1200/interval))];
        test_feat_logVar_prec2 = [zeros(64, nTestTrials, (1200/interval))];
        timeWindowInterval = interval - 1;

        % Power

        for trials = 1:nTestTrials
            TW_int = 1;
            for timeWindow = 1:interval:1200
                test_feat_logVar_power(:,trials,TW_int) = logVarianceFeature(test_Data_power(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
                TW_int = TW_int + 1;
            end
        end

        % Precision 5

        for trials = 1:nTestTrials
            TW_int = 1;
            for timeWindow = 1:interval:1200
                test_feat_logVar_prec5(:,trials,TW_int) = logVarianceFeature(test_Data_precision5(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
                TW_int = TW_int + 1;
            end
        end

        % Precision 2

        for trials = 1:nTestTrials
            TW_int = 1;
            for timeWindow = 1:interval:1200
                test_feat_logVar_prec2(:,trials,TW_int) = logVarianceFeature(test_Data_precision2(:,timeWindow:timeWindow+timeWindowInterval,trials),64);
                TW_int = TW_int + 1;
            end
        end

        % access the features datasets: (channel,trials,time_window)
        %(:,:,1) = 0-2000ms
        %(:,:,2) = 2000-4000ms
        %(:,:,3) = 4000-6000ms
        %(:,:,4) = 6000-8000ms MOVEMENT
        %(:,:,5) = 8000-10000ms
        %(:,:,6) = 10000-12000ms

        %% Next part sets up the dataset for training (inc labels), so that an optimal
        % hyperplane can be established.
        % Could be simplified TODO
        trainingTargetFeatures_power = [];
        trainingTargetFeatures_prec5 = [];
        trainingTargetFeatures_prec2 = [];

        % Power

        trainingTargetFeatures_power(:,1) = train_feat_logVar_power(1,:,4); % First component
        trainingTargetFeatures_power(:,2) = train_feat_logVar_power(64,:,4); % Last component

        % Precision 5

        trainingTargetFeatures_prec5(:,1) = train_feat_logVar_prec5(1,:,4); % First component
        trainingTargetFeatures_prec5(:,2) = train_feat_logVar_prec5(64,:,4); % Last component

        % Precision 2

        trainingTargetFeatures_prec2(:,1) = train_feat_logVar_prec2(1,:,4); % First component
        trainingTargetFeatures_prec2(:,2) = train_feat_logVar_prec2(64,:,4); % Last component

        %% Test target features

        testTargetFeatures_power = [];
        testTargetFeatures_prec5 = [];
        testTargetFeatures_prec2 = [];

        % Power

        testTargetFeatures_power(:,1) = test_feat_logVar_power(1,:,4);
        testTargetFeatures_power(:,2) = test_feat_logVar_power(64,:,4);

        % Precision 5

        testTargetFeatures_prec5(:,1) = test_feat_logVar_prec5(1,:,4);
        testTargetFeatures_prec5(:,2) = test_feat_logVar_prec5(64,:,4);

        % Precision 2

        testTargetFeatures_prec2(:,1) = test_feat_logVar_prec2(1,:,4);
        testTargetFeatures_prec2(:,2) = test_feat_logVar_prec2(64,:,4);

        %% In the Features of the Training Dataset, include the condition of the trial.
        trainingTargetFeatures_power(:,3) = data_with_cond(65,1,cvp.training(fold));

        trainingTargetFeatures_prec5(:,3) = data_with_cond(65,1,cvp.training(fold));

        trainingTargetFeatures_prec2(:,3) = data_with_cond(65,1,cvp.training(fold));

        testTargetFeatures_power(:,3) = data_with_cond(65,1,cvp.test(fold));

        testTargetFeatures_prec5(:,3) = data_with_cond(65,1,cvp.test(fold));

        testTargetFeatures_prec2(:,3) = data_with_cond(65,1,cvp.test(fold));

        % The following adds the correct label in the 3rd column. This is used
        % for the confusion matrix to determine if the predicted class is
        % correct/incorrect. And for the graphing of the different conditions.
        for trial = 1:nTrainingTrials
            % power
            if trainingTargetFeatures_power(trial,3) == 1
                trainingTargetFeatures_power(trial,3) = 1;
            else
                trainingTargetFeatures_power(trial,3) = 0;
            end
            %prec5
            if trainingTargetFeatures_prec5(trial,3) == 2
                trainingTargetFeatures_prec5(trial,3) = 1;
            else
                trainingTargetFeatures_prec5(trial,3) = 0;
            end
            %prec2
            if trainingTargetFeatures_prec2(trial,3) == 3
                trainingTargetFeatures_prec2(trial,3) = 1;
            else
                trainingTargetFeatures_prec2(trial,3) = 0;
            end

        end

        % test trials
        for trial = 1:nTestTrials
            %power
            if testTargetFeatures_power(trial,3) == 1
                testTargetFeatures_power(trial,3) = 1;
            else
                testTargetFeatures_power(trial,3) = 0;
            end
            %prec5
            if testTargetFeatures_prec5(trial,3) == 2
                testTargetFeatures_prec5(trial,3) = 1;
            else
                testTargetFeatures_prec5(trial,3) = 0;
            end
            %prec2
            if testTargetFeatures_prec2(trial,3) == 3
                testTargetFeatures_prec2(trial,3) = 1;
            else
                testTargetFeatures_prec2(trial,3) = 0;
            end
        end

        % 3 models for each CSP-filter, that creates a decision boundary
        % between the conditions.
        md_power = fitcdiscr(trainingTargetFeatures_power(:,1:2), trainingTargetFeatures_power(:,3));
        md_prec5 = fitcdiscr(trainingTargetFeatures_prec5(:,1:2), trainingTargetFeatures_prec5(:,3));
        md_prec2 = fitcdiscr(trainingTargetFeatures_prec2(:,1:2), trainingTargetFeatures_prec2(:,3));

        K_power = md_power.Coeffs(1,2).Const;
        L_power = md_power.Coeffs(1,2).Linear;

        K_prec5 = md_prec5.Coeffs(1,2).Const;
        L_prec5 = md_prec5.Coeffs(1,2).Linear;

        K_prec2 = md_prec2.Coeffs(1,2).Const;
        L_prec2 = md_prec2.Coeffs(1,2).Linear;

        f_power = @(x1,x2) K_power + L_power(1)*x1 + L_power(2)*x2;

        f_prec5 = @(x1,x2) K_prec5 + L_prec5(1)*x1 + L_prec5(2)*x2;

        f_prec2 = @(x1,x2) K_prec2 + L_prec2(1)*x1 + L_prec2(2)*x2;

        %% K-NN
        md_knn_power_csp = fitcknn(trainingTargetFeatures_power(:,1:2), trainingTargetFeatures_power(:,3), 'NumNeighbors', 5, 'Standardize',1);
        md_knn_prec5_csp = fitcknn(trainingTargetFeatures_prec5(:,1:2), trainingTargetFeatures_prec5(:,3), 'NumNeighbors', 5, 'Standardize',1);
        md_knn_prec2_csp = fitcknn(trainingTargetFeatures_prec2(:,1:2), trainingTargetFeatures_prec2(:,3), 'NumNeighbors', 5, 'Standardize',1);


        %% Graphing results
        % The tables of EEG trial data from the 3 different CSP filters are split into the
        % condition group or into the 'rest' group. The rest group are the
        % other two conditions. This is done so that for plotting, these can be
        % coloured differently in the graph.

        power_csp_power_cond = trainingTargetFeatures_power(trainingTargetFeatures_power(:,3) == 1,:);
        power_csp_rest_cond = trainingTargetFeatures_power(trainingTargetFeatures_power(:,3) == 0,:);

        prec5_csp_prec5_cond = trainingTargetFeatures_prec5(trainingTargetFeatures_prec5(:,3) == 1,:);
        prec5_csp_rest_cond = trainingTargetFeatures_prec5(trainingTargetFeatures_prec5(:,3) == 0,:);

        prec2_csp_prec2_cond = trainingTargetFeatures_prec2(trainingTargetFeatures_prec5(:,3) == 1,:);
        prec2_csp_rest_cond = trainingTargetFeatures_prec2(trainingTargetFeatures_prec5(:,3) == 0,:);

        % Plot:
        figure(f1);
        subplot(5,2,fold) % Creates a 5x2 figure for 10 plots. One for each fold.
        title(strcat('Fold', num2str(fold)));
        hold on;
        %scatter(power_csp_power_cond(:,1), power_csp_power_cond(:,2), "red", "o");
        %scatter(power_csp_rest_cond(:,1), power_csp_rest_cond(:,2), "red", 'o', "filled");
        hold on;
        xlabel('First Component');
        ylabel("Last Component");
        % plot power condition as filled circles and non-power as empty circles
        scatter(testTargetFeatures_power(testTargetFeatures_power(:,3) == 1,1),testTargetFeatures_power(testTargetFeatures_power(:,3) == 1,2), "red", 'o', 'filled');
        scatter(testTargetFeatures_power(testTargetFeatures_power(:,3) == 0,1),testTargetFeatures_power(testTargetFeatures_power(:,3) == 0,2), "red", 'o');

        %scatter(prec5_csp_prec5_cond(:,1), prec5_csp_prec5_cond(:,2), "blue", "diamond");
        %scatter(prec5_csp_rest_cond(:,1), prec5_csp_rest_cond(:,2), "blue", 'diamond', "filled");
        hold on;
        xlabel('First Component');
        ylabel("Last Component");
        % plot prec5 condition as filled circles and non-power as empty circles
        scatter(testTargetFeatures_prec5(testTargetFeatures_prec5(:,3) == 1,1),testTargetFeatures_prec5(testTargetFeatures_prec5(:,3) == 1,2), "blue", 'o', 'filled');
        scatter(testTargetFeatures_prec5(testTargetFeatures_prec5(:,3) == 0,1),testTargetFeatures_prec5(testTargetFeatures_prec5(:,3) == 0,2), "blue", 'o');

        %scatter(prec2_csp_prec2_cond(:,1), prec2_csp_prec2_cond(:,2), "green", "square");
        %scatter(prec2_csp_rest_cond(:,1), prec2_csp_rest_cond(:,2), "green", 'square', "filled");
        hold on;
        xlabel('First Component');
        ylabel("Last Component");
        % plot prec2 condition as filled circles and non-power as empty circles
        scatter(testTargetFeatures_prec2(testTargetFeatures_prec2(:,3) == 1,1),testTargetFeatures_prec2(testTargetFeatures_prec2(:,3) == 1,2), "green", 'o', 'filled');
        scatter(testTargetFeatures_prec2(testTargetFeatures_prec2(:,3) == 0,1),testTargetFeatures_prec2(testTargetFeatures_prec2(:,3) == 0,2), "green", 'o');
        %hold off;

        % Plot decision boundaries
        fimplicit(f_power, 'Color', 'red');
        fimplicit(f_prec5, 'Color', 'blue');
        fimplicit(f_prec2, 'Color', 'green');

        md_power_csp_prediction = predict(md_power, testTargetFeatures_power(:,1:2));
        md_prec5_csp_prediction = predict(md_prec5, testTargetFeatures_prec5(:,1:2));
        md_prec2_csp_prediction = predict(md_prec2, testTargetFeatures_prec2(:,1:2));

        % Run the k-nn classifier over the test data
        md_knn_power_prediction = predict(md_knn_power_csp, testTargetFeatures_power(:,1:2));
        md_knn_prec5_prediction = predict(md_knn_prec5_csp, testTargetFeatures_prec5(:,1:2));
        md_knn_prec2_prediction = predict(md_knn_prec2_csp, testTargetFeatures_prec2(:,1:2));

        % Add true labels to matrix
        md_knn_power_prediction(:,2) = testTargetFeatures_power(:,3);
        md_knn_prec5_prediction(:,2) = testTargetFeatures_prec5(:,3);
        md_knn_prec2_prediction(:,2) = testTargetFeatures_prec2(:,3);

        multiclass_prediction = [];
        lda_multiclass_prediction = [];

        for trial = 1:nTestTrials

            % K-NN
            predict_power = md_knn_power_prediction(trial, 1);
            predict_prec5 = md_knn_prec5_prediction(trial, 1);
            predict_prec2 = md_knn_prec2_prediction(trial, 1);

            if predict_power == 1
                multiclass_prediction(trial,1) = 1;
            elseif predict_prec5 == 1
                multiclass_prediction(trial,1) = 2;
            elseif predict_prec2 == 1
                multiclass_prediction(trial,1) = 3;
            else
                multiclass_prediction(trial,1) = 0;
            end


            % Following part makes sure that the true labels are in the same
            % coding as the multiclass_prediction ((0, 1, 2, 3) instead of (0,1))
            % 1 = power
            % 2 = prec5
            % 3 = prec2
            if md_knn_power_prediction(trial,2) == 1
                multiclass_prediction(trial,2) = 1;
                lda_multiclass_prediction(trial, 2) = 1;
            elseif md_knn_prec5_prediction(trial,2) == 1
                multiclass_prediction(trial,2) = 2;
                lda_multiclass_prediction(trial, 2) = 2;
            elseif md_knn_prec2_prediction(trial,2) == 1
                multiclass_prediction(trial,2) = 3;
                lda_multiclass_prediction(trial, 2) = 3;
            else
                multiclass_prediction(trial,2) = 0;
                lda_multiclass_prediction(trial, 2) = 0;
            end


            % LDA Classification rates setup:
            predict_power_lda = md_power_csp_prediction(trial,1);
            predict_prec5_lda = md_prec5_csp_prediction(trial,1);
            predict_prec2_lda = md_prec2_csp_prediction(trial,1);

            % Collect the information of predicted trial type from the 3
            % classifiers.
            if predict_power_lda == 1
                lda_multiclass_prediction(trial,1) = 1;
            elseif predict_prec5_lda == 1
                lda_multiclass_prediction(trial,1) = 2;
            elseif predict_prec2_lda == 1
                lda_multiclass_prediction(trial,1) = 3;
            end


        end

        save_all_results_knn = cat(1, save_all_results_knn, multiclass_prediction);
        save_all_results_lda = cat(1, save_all_results_lda, lda_multiclass_prediction);
    end
end
figure(18);
confusionchart(save_all_results_knn(:,2), save_all_results_knn(:,1))
saveas(18, 'Results_MultiClassClassifier\ConfusionMat_Multi_kNN.png');

figure(19);
confusionchart(save_all_results_lda(:,2), save_all_results_lda(:,1))
saveas(19, 'Results_MultiClassClassifier\ConfusionMat_Multi_lda.png');


