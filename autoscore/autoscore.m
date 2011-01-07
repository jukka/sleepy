function [ S ] = autoscore( eeg, emg, score )
%AUTOSCORE Automatic sleep scorer based on raw EEG and EMG data
%   This function attempts to automatically classify different
%   sleep stages (wake, NREM, REM) from the given EEG and EMG
%   data and a set of initial manual scoring results.
%
%   The input variables should be EEG, EMG and score channel data
%   structures exported from the Spike2 application using the export
%   to Matlab feature. A sample EEG or EMG channel looks like this:
%
%         title: 'EEG'
%       comment: ''
%      interval: 0.0050
%         scale: 1.5259e-004
%        offset: 0
%         units: ''
%         start: 0
%        length: 3079419
%        values: [3079419x1 double]
%
%   The score channel should be already be divided into appropriately
%   sized epochs (typically a few seconds). This function will then
%   return a vector containing the automatic scoring results for each
%   epoch. A sample score channel looks like this:
%
%         title: 'Epochs'
%       comment: 'Wake|NREM|REM|ARTE|'
%    resolution: 0.0050
%        length: 19826
%         items: 44
%         times: [19826x1 double]
%         codes: [19826x4 uint8]
%          text: [19826x44 char]
%
%   The automatic scoring algorithm is based on the signal strengths
%   of each EEG frequency channel (sigma, delta, theta, beta) and the
%   total EMG signal power per each epoch. This preprocessed
%   five-dimensional data is then reduced to two dimensions using
%   principal component analysis and then classified to score classes
%   using the closes neighbour algorithm based on initial teaching data
%   from the given score channel.
%

% First define/calculate some constants
epoch_duration = score.times(2) - score.times(1);
epoch_length = floor(epoch_duration / eeg.interval);
fft_length = 2 ^ nextpow2(epoch_length);
fft_slots = linspace(0, 1 / eeg.interval, fft_length);
hamming_window = hamming(epoch_length);
epoch_count = min(score.length, floor(eeg.length / epoch_length));

% Then initialize our main variables
X = zeros(epoch_count, 5);     % Normalized EEG/EMG data per epoch
Y = zeros(epoch_count, 2);     % Data reduced to 2 dimensions by PCA
S = zeros(epoch_count, 1);     % Result of automatic scoring

% Normalize the raw EEG/EMG data to characteristic values
for epoch = (1 : epoch_count)
    % Find the measurements in the time window of this epoch
    start_time = score.times(epoch);
    start_index = floor(start_time / eeg.interval) + 1;
    end_index = start_index + epoch_length - 1;
    if end_index > eeg.length
        break
    end
    eeg_window = eeg.values(start_index : end_index);
    emg_window = emg.values(start_index : end_index);

    % Convert from amplitude per time to power per frequency and
    % calculate signal strengths per interesting EEG frequency bands
    power = abs(fft(eeg_window .* hamming_window, fft_length))' .^ 2;
    sigma = (fft_slots >= 0.4) .* (fft_slots < 4.0) .* power;
    delta = (fft_slots >= 4.0) .* (fft_slots < 10.0) .* power;
    theta = (fft_slots >= 10.0) .* (fft_slots < 15.0) .* power;
    beta = (fft_slots >= 15.0) .* (fft_slots < 20.0) .* power;
    total = sigma + delta + theta + beta;

    % Add relative powers per EEG frequency band and the total EMG power
    X(epoch, 1) = sigma / total;
    X(epoch, 2) = delta / total;
    X(epoch, 3) = theta / total;
    X(epoch, 4) = beta / total;
    X(epoch, 5) = sum(emg_window .^ 2);
end
emg_sums = X(:, 5);
X(:, 5) = emg_sums / (max(emg_sums) - min(emg_sums));

% Reduce the dimension of our data from 5 to 2 using
% principal component analysis
X = X - ones(epoch_count, 1) * (sum(X) / epoch_count);
[ E, t ] = eigs(cov(X), 2);
display(E);
display(t);
Y = X * E;

% Use a random sample of pre-scored epochs as the basis of autoscoring
% TODO: Find a better way for the user to enter a small set of base scores
manual_score = score.codes(1 : epoch_count, 1);
score_count = max(manual_score) + 1;
sample_count = 30; % Use 30 random samples per each score category
sample_epochs = [];
sample_Y = [];
sample_scores = [];
n = 1;
for s = 0 : score_count - 1
    matching_epochs = find(manual_score == s);
    matching_epoch_count = length(matching_epochs);
    if matching_epoch_count > 0
        for sample = 1 : sample_count
            sample_epochs(n) = matching_epochs(randi(matching_epoch_count));
            sample_Y(n, :) = Y(sample_epochs(n), :);
            sample_scores(n) = s;
            n = n + 1;
        end
    end
end

% Use the closest neighbour method to autoscore each epoch
distances = zeros(length(sample_epochs), 1);
for epoch = 1 : epoch_count
    epoch_vector = Y(epoch, :);
    for i = 1 : length(sample_epochs)
        sample_epoch = sample_epochs(i);
        sample_vector = Y(sample_epoch, :);
        distances(i) = sum((epoch_vector - sample_vector) .^ 2);
    end
    scores = zeros(score_count, 1);
    [ ~, index ] = sort(distances);
    for i = 1 : 15
        s = sample_scores(index(i)) + 1;
        scores(s) = scores(s) + 1;
    end
    [ ~, index ] = sort(scores);
    S(epoch) = index(score_count) - 1;
end

% colormap([ 0.5 0 0; 0 0.5 0; 0 0 0.5; 1 1 0 ]);
subplot(3, 1, 1);
scatter(Y(:, 1), Y(:, 2), [], manual_score);
subplot(3, 1, 2);
scatter(sample_Y(:, 1), sample_Y(:, 2), [], sample_scores);
subplot(3, 1, 3);
scatter(Y(:, 1), Y(:, 2), [], S);

return;
subplot(3, 1, 1);
scatter3(Y(:, 1), Y(:, 2), Y(:, 3), [], manual_score);
subplot(3, 1, 2);
scatter3(sample_Y(:, 1), sample_Y(:, 2), sample_Y(:, 3), [], sample_scores);
subplot(3, 1, 3);
scatter3(Y(:, 1), Y(:, 2), Y(:, 3), [], S);
end
