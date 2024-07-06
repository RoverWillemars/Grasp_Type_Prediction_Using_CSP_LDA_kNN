subject_data = load('-mat', strcat('Subj_data\S3\all_EEG.mat'));

% This file gets the raw data of the EEG and produces a txt file for each
% subject that will be used by the main script to determine which condition
% the trial is. Restructuring for the importing of events in EEGLab.
for subject = 3:18
    if(subject == 4 || subject == 8)
        
        subject_data = load('-mat', strcat('Subj_data\S', string(subject), '\all_EEG.mat'));

        combine_cond = [subject_data.indexes_LC_power; subject_data.indexes_LC_precision2; 
            subject_data.indexes_LC_precision5; subject_data.indexes_LS_power;
            subject_data.indexes_LS_precision2; subject_data.indexes_LS_precision5;
            subject_data.indexes_SC_power; subject_data.indexes_SC_precision2;
            subject_data.indexes_SC_precision5; subject_data.indexes_SS_power;
            subject_data.indexes_SS_precision2; subject_data.indexes_SS_precision5];

        cond_list(1,3) = "latency"; 
        cond_list(1,2) = "type";
        cond_list(1,1) = "epoch";

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(1,i);
            cond_list(tempCondPos + 1, 2) = "LC_Power";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(2,i);
            cond_list(tempCondPos + 1, 2) = "LC_Precision_2";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(3,i);
            cond_list(tempCondPos + 1, 2) = "LC_Precision_5";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(4,i);
            cond_list(tempCondPos + 1, 2) = "LS_Power";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(5,i);
            cond_list(tempCondPos + 1, 2) = "LS_Precision_2";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(6,i);
            cond_list(tempCondPos + 1, 2) = "LS_Precision_5";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(7,i);
            cond_list(tempCondPos + 1, 2) = "SC_Power";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(8,i);
            cond_list(tempCondPos + 1, 2) = "SC_Precision_2";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(9,i);
            cond_list(tempCondPos + 1, 2) = "SC_Precision_5";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(10,i);
            cond_list(tempCondPos + 1, 2) = "SS_Power";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(11,i);
            cond_list(tempCondPos + 1, 2) = "SS_Precision_2";
        end

        for i = 1:length(combine_cond)
            tempCondPos = combine_cond(12,i);
            cond_list(tempCondPos + 1, 2) = "SS_Precision_5";
        end

        % Create a new column for the placement of the events, start from row 2 to keep
        % the column names.
        p = 1;
        for i = 1:length(cond_list) - 1
            cond_list(i + 1,1) = p;
            p = p + 1;
        end

        % include latencies of each trial for event import:
        p = 0;

        for i = 1:length(cond_list) - 1
            cond_list(i + 1, 3) = p;
            p = p + 12000;
        end

        % Save to txt file, so that it can be accessed by EEGlab
        writematrix(cond_list, strcat('sub_trial_info\TrialInfo_sub', string(subject), '.txt'));

    elseif(subject == 5)
        % Empty file, no data on the subject.
    else
        subject_data = load('-mat', strcat('Subj_data\S', string(subject), '\all_EEG.mat'));

        subject_trial_info = load('-mat', strcat('Subj_data\S', string(subject), '\indexes_trials.mat'));

        trial_time = 0;

        cond_list(1,3) = "latency"; 
        cond_list(1,2) = "type";
        cond_list(1,1) = "epoch";

        % Collect the index of the trial from
        % subject_trial_info.all_indexes and place it into cond_list
        for trial = 1:length(subject_trial_info.all_indexes)
            cond_list(trial+1, 2) = string(subject_trial_info.all_indexes(trial));
            cond_list(trial+1, 3) = string(trial_time);
            trial_time = trial_time + 12000;
        end
        
        %Set first column to trial number
        p = 1;
        for i = 1:length(cond_list) - 1
            cond_list(i + 1,1) = p;
            p = p + 1;
        end

        % Change the numerical value representating the condition into the
        % string form. Easier to read and it is how the rest of the script
        % expects the condition to be represented.
        for trial = 2:length(cond_list)
            if cond_list(trial, 2) == '1'
                cond_list(trial, 2) = "LC_Power";

            elseif cond_list(trial,2) == '2'
                cond_list(trial, 2) = 'LC_Precision_5';

            elseif cond_list(trial,2) == '3'
                cond_list(trial, 2) = 'LC_Precision_2';

            elseif cond_list(trial,2) == '4'
                cond_list(trial, 2) = 'SC_Power';

            elseif cond_list(trial,2) == '5'
                cond_list(trial, 2) = 'SC_Precision_5';

            elseif cond_list(trial,2) == '6'
                cond_list(trial, 2) = 'SC_Precision_2';

            elseif cond_list(trial,2) == '7'
                cond_list(trial, 2) = 'LS_Power';

            elseif cond_list(trial,2) == '8'
                cond_list(trial, 2) = 'LS_Precision_5';

            elseif cond_list(trial,2) == '9'
                cond_list(trial, 2) = 'LS_Precision_2';

            elseif cond_list(trial,2) == '10'
                cond_list(trial, 2) = 'SS_Power';

            elseif cond_list(trial,2) == '11'
                cond_list(trial, 2) = 'SS_Precision_2';

            elseif cond_list(trial,2) == '12'
                cond_list(trial, 2) = 'SS_Precision_5';

            else
                cond_list(trial, 2) = 'NA';
            
            end
        end

        writematrix(cond_list, strcat('sub_trial_info\TrialInfo_sub', string(subject), '.txt'));
    end
end 


% indexes_LC_power = 1
% indexes_LC_precision2 = 3
% indexes_LC_precision5 = 2
% indexes_LS_power = 7
% indexes_LS_precision2 = 9
% indexes_LS_precision5 = 8
% indexes_SC_power = 4
% indexes_SC_precision2 = 6
% indexes_SC_precision5 = 5
% indexes_SS_power = 10
% indexes_SS_precision2 = 12
% indexes_SS_precision5 = 11
