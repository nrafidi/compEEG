%Gets the average theta band power of the given trial and returns it in
%theta. The width dimension of trial is channel, and the length is time
%Dependencies: none
function theta = getTheta(trial)

start = 4; %Starting frequency
bw = 3; %Bandwidth to be examined
[length width] = size(trial);

theta = zeros(1, width);

for i = 1:width
    trans = abs(fft(trial(:, i)));
    theta(1, i) = mean(trans(start:start+bw));
    if isnan(theta(1, i))
        theta(1, i) = 0;
    end
end