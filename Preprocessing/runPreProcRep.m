subjects = {'YY', 'WW'};

% {'BBB', 'GGG', 'HHH', 'AAA', 'CCC', 'DDD', 'EEE', 'FFF', 'MM', ...
%      'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'WW',...
%      'YY','III', 'JJJ', 'KKK'};
%YY and WW Need to be rerun
 

numSub = length(subjects);

for s = 1:numSub
    sub = subjects{s};
    disp(sub);
    options = struct;
    options.isVis = true;
    if strcmp(sub, 'NN') 
        options.LP = 127;
    elseif strcmp(sub, 'HHH')
        options.isVis = false;
    end
%     preprocPipeline(sub, 'CompEEG', options);
%     disp('CompEEG Complete');
    preprocPipeline(sub, 'CompEEG__KR', options);
    disp('CompEEG__KR Complete');
end