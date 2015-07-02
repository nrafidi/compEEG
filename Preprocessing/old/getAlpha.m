%Gets the average alpha band power of the given trial and returns it in
%alpha. The width dimension of trial is channel, and the length is time.
%Dependencies: none
function alpha = getAlpha(trial)

start = 8; %Start of alpha band
bw = 4; %Bandwidth to be examined

[length width] = size(trial);
alpha = zeros(1, width);

for i = 1:width
    trans = abs(fft(trial(:, i)));
    alpha(1, i) = mean(trans(start:start+bw));
    if isnan(alpha(1, i))
        alpha(1, i) = 0;
    end
end