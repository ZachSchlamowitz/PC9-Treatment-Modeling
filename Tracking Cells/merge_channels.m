%% merge_channels.m
% Script to sum signal from G1 (cdt1) reporter and S/G2 (geminin) reporter
% to reconstrue a nuclear marker channel to aid in tracking Fucci clonal
% cells

% Author: Zach Schlamowitz (8/10/20)

%% Main
cd 'D:\Zach\2021-07-26\RAW_DATA'  % point MATLAB to the RAW_DATA folder contained desired tifs
    % From flashdrive (Seagull): 'F:\2021-07-26 Partial Copy\RAW_DATA'
    % From Lisa Shanks' Harddrive (Seagate Expansion Desk): 'D:\Zach\2021-07-26\RAW_DATA'
date = '2021-07-26-';
channel_a = 1; % channel numbers for desired channels
channel_b = 3;
num_positions = 120;
num_timepoints = 289; % number of timepoints for one channel at one position

for xy = 1:num_positions
    for t=1:num_timepoints
        filename_a = horzcat(date,'xy',sprintf('%03d',xy),'c',num2str(channel_a),'t',sprintf('%03d',t));
        filename_b = horzcat(date,'xy',sprintf('%03d',xy),'c',num2str(channel_b),'t',sprintf('%03d',t));
        data_a = imread(filename_a, 'tif');
        data_b = imread(filename_b, 'tif');
        sum = data_a + data_b;
        
        out_file_name = horzcat(date,'xy',sprintf('%03d',xy),'c4','t',sprintf('%03d',t),'.tif');
        imwrite(sum,out_file_name,'tif');
    end
end


% Andrew's commands to view images side-by-side
% close all
% figure, imshow(sum,[])
% title('sum')
% 
% figure, imshow(data_a,[])
% title('data_a')
% 
% figure, imshow(data_b,[])
% title('data_b')
% 
% imtool(sum) % allows contrast adjustment
