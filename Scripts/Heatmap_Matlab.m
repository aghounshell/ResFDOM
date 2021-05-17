% Script to plot heatmaps of Temp and DO for 2015-2020
% 17 May 2021  A Hounshell

% Load in Temp and DO data (merged YSI and CTD files from Heatmaps.R)
temp = importdata('Heatmap_temp.csv');
time = temp.textdata(2:358, 1);
depth = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.3 1.5 1.6 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.2 3.3 3.5 3.8 4 4.1 4.3 4.4 4.5 4.7 5 5.3 5.5 5.6 5.7 5.9 6 6.2 6.3 6.5 6.8 7 7.1 7.3 7.4 7.5 7.7 7.8 8 8.3 8.5 8.6 8.7 8.8 8.9 9 9.2 9.3 9.4 9.5 9.6 10]';

do = importdata('Heatmap_DO.csv');

% Try to plot using contourf
figure
contourf(temp.data)
colorbar()