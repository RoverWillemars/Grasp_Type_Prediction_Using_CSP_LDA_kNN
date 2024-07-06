function logVar = logVarianceFeature(trials, Nchannels)

% Returns the log of the variance of the EEG signal.
% Input a 2d-array (channels x samples). This means that in the main
% script, you have to run through each trial 1:nTrials.
% Outputs a 2d-array (Channels x trials)
    logVar = [zeros(Nchannels, 1)];

    for channel = 1:Nchannels
        logVar(channel) = log(var(trials(channel,:)));
    end

end