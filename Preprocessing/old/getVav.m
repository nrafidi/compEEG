%Gets the average voltage of the given trial and returns it in vAv
%The width dimension of trial is channel, and the length is time
%Dependencies: none
function vAv = getVav(trial)

[length width] = size(trial);

vAv = zeros(1, width);

for i = 1:width
    vAv(1, i) = mean(trial(:, i));
end
        
end