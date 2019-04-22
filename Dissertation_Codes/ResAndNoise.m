function [] = ResAndNoise()
% dir ~ Directory
% rp ~ Reference pipeline
% dsp ~ Downsampling pipelines
% d3d ~ Downsampled dream.3d file names
% Notes:
% Average grain size is hard-coded
% File Directory is hard-coded
% Resolution and Noise levels should be modified based on the experiment
% being performed
 
 
clear all; close all; clc; format compact
dir = sprintf('%s', 'cd C:\Users\Greg\Documents\__rx114data\2015_MarChar_Noise_Random');
 
% % Build Phantom & Compute Stats 
rp= sprintf('%s', 'Phantom_Build_Stats.json');
PipelineRunner(dir,rp)
 
%% Create and Run Pipelines for Error Analysis
res3 = linspace(1.0,1.0,4); %3 VRAD Results
res5 = linspace(0.6,0.6,4); %5 VRAD Results
res10 = linspace(0.3,0.3,4); %10 VRAD Results
res15 = linspace(0.2,0.2,4); %15 VRAD Results
res = [res3 res5 res10 res15];
 
%NOISE
pn_repeat = 0.01*[10 25 50 75]; % Random Noise Levels to be investigated
bn_repeat = linspace(0.0,0.0,4); % Boundary Noise Levels to be investigated
pn = [pn_repeat pn_repeat pn_repeat pn_repeat]; %Repeat for each res
bn = [bn_repeat bn_repeat bn_repeat bn_repeat]; %Repeat for each res
 
k = 20; % Instantiation Number (# structures at each noise level)
for i = 1:k
    instant_num = num2str(i);
    
    for j = 1:length(res)
        
        AvgGrainSize = 3;
        A = sprintf('%i',(AvgGrainSize/res(j)));
        B = sprintf('%i',(pn(j)*100));
        C = sprintf('%i',(bn(j)*100));
        
        % New Pipeline .txt File Names
        dsp{j,1} = strcat(sprintf('%s', '_'), num2str(A), sprintf('%s', 'VRAD'), ...
                          sprintf('%s', '_'), num2str(B), sprintf('%s', 'pn'), ...
                          sprintf('%s', '_'), num2str(C), sprintf('%s', 'bn'), ...
                          sprintf('%s', '_'), instant_num, sprintf('%s', '.json')); 
 
        % DREAM.3D File Names
        d3d{j,1} = strcat(sprintf('%s', '_'), num2str(A), sprintf('%s', 'VRAD'), ...
                          sprintf('%s', '_'), num2str(B), sprintf('%s', 'pn'), ...
                          sprintf('%s', '_'), num2str(C), sprintf('%s', 'bn'), ...
                          sprintf('%s', '_'), instant_num, sprintf('%s', '.dream3d')); 
 
        % CSV File Names
        csv{j,1} = strcat(sprintf('%s', '_'), num2str(A), sprintf('%s', 'VRAD'), ...
                          sprintf('%s', '_'), num2str(B), sprintf('%s', 'pn'), ...
                          sprintf('%s', '_'), num2str(C), sprintf('%s', 'bn'), ...
                          sprintf('%s', '_'), instant_num, sprintf('%s', '.csv')); 
 
        % Run New Pipeline
        NewRes = sprintf('%2.2f',(res(j)));
        NewPN = sprintf('% 2.2f',(pn(j)));
        NewBN = sprintf('% 2.2f',(bn(j)));
        PipelineCreator('DownSample.json',dsp{j},d3d{j},csv{j},NewRes,NewPN,NewBN);
        PipelineRunner(dir,dsp{j});
    
    end
    
end
