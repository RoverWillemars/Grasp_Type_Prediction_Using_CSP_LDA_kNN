function filteredTrial = butterBandpass(trials, low, high, sample_rate, order, Nchannels)

% When two elements are specified in the [], the filter type becomes
% bandpass. Calculates the coeffients for the bandpass filter.

    [b, a] = butter(order, [low, high]/(sample_rate/2), 'bandpass');

    filteredTrial = [zeros(64, 1200)];

    for channel = 1:Nchannels
        seg = trials(channel,:,:);
        seg = filtfilt(b,a,double(seg));

        filteredTrial(channel, :) = seg;
    end


end