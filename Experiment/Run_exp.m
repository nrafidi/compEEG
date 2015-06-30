% original wiki: http://www.ccn.fhucla.edu/wiki/index.php/Psychtoolbox
% See README for instructions on use
% 
% try %EEG stim computer at PNI
%     ls('C:\Users\EEGLab\Desktop\CompetitionEEG\BCI\');
%     WORKING_DIR = 'C:\Users\EEGLab\Desktop\CompetitionEEG\BCI\';
% catch
%     ls('~/Desktop/CompetitionEEG/BCI/');
%     WORKING_DIR = '~/Desktop/CompetitionEEG/BCI/';
%     %error('Can''t find working directory');
% end
% cd(WORKING_DIR)

clear;
addpath expr/
KbName('UnifyKeyNames');

prompt={'Experiment Name','Subject Name', 'At Scanner? Yes=1, No=0', 'File Number',...
    'Block Number'};
def={'KR','A','0','2','1'};
title = 'SETUP EXP';
lineNo=1;

answer=inputdlg(prompt,title,lineNo,def);
exp_name = char(answer(1,:));
par.subject = char(answer(2,:));
par.atscanner = str2num(char(answer(3,:))); %#ok<ST2NM>
par.fileNum = str2num(char(answer(4,:))); %#ok<ST2NM>
par.blockNum = str2num(char(answer(5,:))); %#ok<ST2NM>
% TR = str2num(char(answer(5,:))); %#ok<ST2NM>

par.pportTime       = 0.1;
par.FixTime         = 5;
par.expname        = sprintf('%s_%s_%s',exp_name, num2str(par.blockNum),num2str(par.fileNum));
par.bgcolor        = 255;
par.txcolor        = 0;
par.shift           = 10;

%% This sets up the file naming convention and checks if you're going to
%% overwrite a subject.
filename = sprintf('%s_%s_events_redo.txt',par.subject,par.expname);
fid = fopen(filename);
while fid ~= -1
    fclose(fid);
    disp 'this file already exists ...';
    par.fileNum = input('enter new file num');
    par.expname = sprintf('%s_%s',exp_name, num2str(par.fileNum));
    filename = sprintf('%s_%s_events.txt',par.subject,par.expname);
    fid = fopen(filename);
end

%% Outputs the file (this is actually done within PresText if there is an
%% error)

%Loading stimulus files for that experiment and that subject
if strcmp(exp_name, 'Comp');
    if par.fileNum == 1
        load(['./' par.subject '/presentation.mat']);
    else
        load(['./' par.subject '/competition.mat']);
    end
else
    if par.fileNum == 1
        load(['./' par.subject '/KRpres.mat']);
    else
        load(['./' par.subject '/KRtest.mat']);
    end
end
par.story = experiment(par.blockNum).story;
if isfield(experiment(par.blockNum), 'answer')
    par.answer = experiment(par.blockNum).answer;
end
par.storyLength = experiment(par.blockNum).storyLength;
par.storyTime = experiment(par.blockNum).storyTime;
par.storyTime = par.storyTime+0.5;
par.parPort = experiment(par.blockNum).parPort;

if strcmp(exp_name, 'Comp') || par.fileNum == 1;
    RTs=PresentStim(par);
else
    [RTs, corr] = PresentStimQ(par);
    save(['./' par.subject '/KR_results.mat'], 'RTs', 'corr');
    
% Code for shrinking stimulus set: unused
%     numNew = sum(~corr)*2;
%     nextBlock = par.blockNum+1;
%     experiment(nextBlock).story = cell(1, numNew);
%     experiment(nextBlock).answer = cell(1, numNew);
%     experiment(nextBlock).storyLength = zeros(1,numNew);
%     experiment(nextBlock).parPort = zeros(1,numNew);
%     
%     ind = 3:2:(numNew-1);
%     jInd = [1 ind(randperm(numNew/2-1))];
%     k = 1;
%     for i = 1:length(corr)
%         if ~corr(i)
%             j = jInd(k);
%             experiment(nextBlock).story{j} = par.story{i}; %#ok<*SAGROW>
%             experiment(nextBlock).story{j+1} = par.story{i+1};
%             experiment(nextBlock).answer{j} = par.answer{i};
%             experiment(nextBlock).answer{j+1} = par.answer{i+1};
%             experiment(nextBlock).storyLength(j) = par.storyLength(i);
%             experiment(nextBlock).storyLength(j+1) = par.storyLength(i+1);
%             experiment(nextBlock).parPort(j) = par.parPort(i);
%             experiment(nextBlock).parPort(j+1) = par.parPort(i+1);
%             k = k+1;
%         end
%     end
%     experiment(nextBlock).storyTime = cumsum([0, experiment(nextBlock).storyLength]);
%     save(['./' par.subject '/KRtest.mat'], 'experiment');
%     
%     %Update study stimuli
%     load(['./' par.subject '/KRpres.mat']);
%     
%     experiment(nextBlock).story = cell(1, numNew);
%     experiment(nextBlock).answer = cell(1, numNew);
%     experiment(nextBlock).storyLength = zeros(1,numNew);
%     experiment(nextBlock).parPort = zeros(1,numNew);
%     kInd = [1 ind(randperm(numNew/2-1))];
%     l = 1;
%     for i = 1:length(corr)
%         if ~corr(i)
%             for j = 1:length(par.parPort)
%                 if experiment(par.blockNum).parPort(j) == par.parPort(i);
%                     k = kInd(l);
%             experiment(nextBlock).story{k} = experiment(par.blockNum).story{j}; %#ok<*SAGROW>
%             experiment(nextBlock).story{k+1} = experiment(par.blockNum).story{j+1};
%             experiment(nextBlock).storyLength(k) = experiment(par.blockNum).storyLength(j);
%             experiment(nextBlock).storyLength(k+1) = experiment(par.blockNum).storyLength(j+1);
%             experiment(nextBlock).parPort(k) = experiment(par.blockNum).parPort(j);
%             experiment(nextBlock).parPort(k+1) = experiment(par.blockNum).parPort(j+1);
%             l = l+1;
%                 end
%             end
%         end
%     end
%     experiment(nextBlock).storyTime = cumsum([0, experiment(nextBlock).storyLength]);
%     save(['./' par.subject '/KRpres.mat'], 'experiment');
end


