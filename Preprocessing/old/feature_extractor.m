%Extracts the features of the EEG timecourse stored in data, assumed to be
%of width 15, with the first column as time in ms (0 value at stimulus
%onset), the next 14 containing data for 14 channels. Stimlabel contains
%class labels of each data point (1 for comp, 2 for non-comp). labels is a
%cell array of strings labelling each feature. If write is true,
%feature_extractor will write the extracted features to filename.xls, with
%labels as the first row, and with the class labels as letters, for direct
%input into weka. feature is an all numbers, label-less representation.
%Each column is a feature, and the rows are trials from the data set.
%Dependencies: getVav, getTheta, getAlpha, getBeta
function feature = feature_extractor(filename, stimlabel, data, labels, write)


margin_first = 200; %End of T1 for frequency band, beginning of T2
margin_second = 400; %End of T2 for frequency band, beginning of T3
margin_third = 600; %End of T3 for frequency band
margin_erp1 = 100; %Start of ERP-T1
margin_erp2 = 150; %End of ERP-T1, Start of ERP-T2
margin_erp3 = 200; %End of ERP-T2, Start of ERP-T3
margin_erp4 = 250; %End of ERP-T3, Start of ERP-T4
margin_erp5 = 300; %End of ERP-T4, Start of ERP-T5
margin_erp6 = 350; %End of ERP-T5, Start of ERP-T6
margin_erp7 = 400; %End of ERP-T6

trial_freq_first = [];
trial_freq_second = [];
trial_freq_third = [];
trial_erp_first = [];
trial_erp_second = [];
trial_erp_third = [];
trial_erp_fourth = [];
trial_erp_fifth = [];
trial_erp_sixth = [];
feature = [];
feature_inst = zeros(1, 211);
feature_inst(1, 211) = 1;
onset = false;
prev = -1;
for i = 1:length(data)
    if (data(i, 1) == 0 && prev ~= 0)
        onset = true;
    end
    prev = data(i, 1);
    if onset && (data(i,1) < margin_first)
        trial_freq_first = [trial_freq_first; data(i, :)];
    end
    if onset && (data(i,1) > margin_first) && (data(i,1) < margin_second)
        trial_freq_second = [trial_freq_second; data(i, :)];
    end
    if onset && (data(i,1) > margin_second) && (data(i,1) < margin_third)
        trial_freq_third = [trial_freq_third; data(i, :)];
    end
    if onset && (data(i,1) > margin_erp1) && (data(i,1) < margin_erp2)
        trial_erp_first = [trial_erp_first; data(i, :)];
    end
    if onset && (data(i,1) > margin_erp2) && (data(i,1) < margin_erp3)
        trial_erp_second = [trial_erp_second; data(i, :)];
    end
    if onset && (data(i,1) > margin_erp3) && (data(i,1) < margin_erp4)
        trial_erp_third = [trial_erp_third; data(i, :)];
    end
    if onset && (data(i,1) > margin_erp4) && (data(i,1) < margin_erp5)
        trial_erp_fourth = [trial_erp_fourth; data(i, :)];
    end
    if onset && (data(i,1) > margin_erp5) && (data(i,1) < margin_erp6)
        trial_erp_fifth = [trial_erp_fifth; data(i, :)];
    end
    if onset && (data(i,1) > margin_erp6) && (data(i,1) < margin_erp7)
        trial_erp_sixth = [trial_erp_sixth; data(i, :)];
    end
    if onset && (data(i,1) > margin_third)
        feature_inst(1, 1:14) = getVav(trial_erp_first(:, 2:15));
        feature_inst(1, 15:28) = getVav(trial_erp_second(:, 2:15));
        feature_inst(1, 29:42) = getVav(trial_erp_third(:, 2:15));
        feature_inst(1, 43:56) = getVav(trial_erp_fourth(:, 2:15));
        feature_inst(1, 57:70) = getVav(trial_erp_fifth(:, 2:15));
        feature_inst(1, 71:84) = getVav(trial_erp_sixth(:, 2:15));
        feature_inst(1, 85:98) = getTheta(trial_freq_first(:, 2:15));
        feature_inst(1, 99:112) = getAlpha(trial_freq_first(:, 2:15));
        feature_inst(1, 113:126) = getBeta(trial_freq_first(:, 2:15));
        feature_inst(1, 127:140) = getTheta(trial_freq_second(:, 2:15));
        feature_inst(1, 141:154) = getAlpha(trial_freq_second(:, 2:15));
        feature_inst(1, 155:168) = getBeta(trial_freq_second(:, 2:15));
        feature_inst(1, 169:182) = getTheta(trial_freq_third(:, 2:15));
        feature_inst(1, 183:196) = getAlpha(trial_freq_third(:, 2:15));
        feature_inst(1, 197:210) = getBeta(trial_freq_third(:, 2:15));
        
        feature_inst(1, 211:225) = [mean(feature_inst(1, 1:4)) mean(feature_inst(1, 15:18)) mean(feature_inst(1, 29:32)) mean(feature_inst(1, 43:46)) mean(feature_inst(1, 57:60)) mean(feature_inst(1, 71:74)) mean(feature_inst(1, 85:88)) mean(feature_inst(1, 99:102)) mean(feature_inst(1, 113:116)) mean(feature_inst(1, 127:130)) mean(feature_inst(1, 141:144)) mean(feature_inst(1, 155:158)) mean(feature_inst(1, 169:173)) mean(feature_inst(1, 183:186)) mean(feature_inst(1, 197:200))];
        feature_inst(1, 226:240) = [mean(feature_inst(1, 11:14)) mean(feature_inst(1, 25:28)) mean(feature_inst(1, 39:42)) mean(feature_inst(1, 53:56)) mean(feature_inst(1, 67:70)) mean(feature_inst(1, 81:84)) mean(feature_inst(1, 95:98)) mean(feature_inst(1, 109:112)) mean(feature_inst(1, 123:126)) mean(feature_inst(1, 137:140)) mean(feature_inst(1, 151:154)) mean(feature_inst(1, 165:168)) mean(feature_inst(1, 179:183)) mean(feature_inst(1, 193:196)) mean(feature_inst(1, 207:210))];
        feature_inst(1, 241) = stimlabel(i);
        
            
        
        trial_freq_first = [];
        trial_freq_second = [];
        trial_freq_third = [];
        trial_erp_first = [];
        trial_erp_second = [];
        trial_erp_third = [];
        trial_erp_fourth = [];
        trial_erp_fifth = [];
        trial_erp_sixth = [];
        onset = false;
        feature = [feature; feature_inst];
    end
    
end
if (write) 
    xlswrite(filename, labels, 1, 'A1');
    xlswrite(filename, feature, 1, 'A2');
    [le wi] = size(feature);
    for f = 1:le
        if (feature(f, wi) == 1)
            xlswrite(filename, 'C', 1, strcat('IG', num2str(f + 1)));
        else
            xlswrite(filename, 'N', 1, strcat('IG', num2str(f + 1)));
        end
    end
end
end