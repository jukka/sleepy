function power = supersleep(channel)
    bands = [
         0.5   4.0
         5.0   8.0
        11.0  15.0
        30.0  45.0
         1.0  35.0
         0.0   1.0
        ];
    power = divideToEpochs(channel, 5, -1, -1);
end

function epochs = divideToEpochs(channel, secondsPerEpoch, startTime, endTime)
    if startTime < channel.start
        startTime = channel.start;
    end
    if endTime < 0 || endTime > channel.start + channel.lenght * channel.interval
        endTime = channel.start + channel.length * channel.interval;
    end
    m = ceil(secondsPerEpoch / channel.interval);
    n = ceil((endTime - startTime) / secondsPerEpoch);
    epochs = zeros(m, n);
    for i = 1 : n
        p = (i - 1) * m + 1;
        q = min(i * m + 1, channel.length);
        epochs(1 : (q - p), i) = channel.values(p : q - 1);
    end
end

function powerPerBand = calculatePowerPerBand(amplitude, bands)
    n = size(bands, 1);
    power = (real(amplitude).^2 + imag(amplitude).^2) / 2;
    powerPerBand = zeros(n, 1);
    for i = 1 : n
        powerPerBand(i) = sum(power(bands(:,i)));
    end
end

function band = createFrequencyBand(frequencies, min, max)
    % Creates a logical vector for selecting the specified frequency range
    band = frequencies >= min & frequencies <= max;
end
