%Gets the average beta band power of the given trial and returns it in
%betaa. The width dimension of trial is channel, and the length is time
%Dependencies: none
function beta = getBeta(trial)

start = 13; %Start of the beta band
%Beta band is counted as anything 13 Hz or greater
[length width] = size(trial);

beta = zeros(1, width);

for i = 1:width
    trans = abs(fft(trial(:, i)));
    beta(1, i) = mean(trans(start:length));
    if isnan(beta(1, i))
        beta(1, i) = 0;
    end
end