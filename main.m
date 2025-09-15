%% Pixel-Level Fusion of SAR Image Despeckling Algorithms
% 
% Dr. Emrah Onat
% 16.09.2025


% Despeckling Algorithms
% 01 - Mean Filter
% 02 - Median Filter
% 03 - Frost Filter
% 04 - Mod. Frost Filter
% 05 - Lee Filter
% 06 - Kuan Filter
% 07 - Kuwahara Filter
% 08 - SDD Filter
% 09 - SDD-QL Filter
% 10 - SARBM3D Filter
% 11 - NCDF Filter
% 12 - SRAD Filter
% ...
% XX - FUSION Filter

% Metrics
% 01 - PSNR
% 02 - SSIM
% 03 - EQP
% 04 - Entrophy
% 05 - 

close all; clearvars; clc;

number_of_filter = 12;
colormap(1,:) = [1 0 0]; 
colormap(2,:) = [0 1 0];
colormap(3,:) = [0 0 1];
colormap(4,:) = [0 1 1];
colormap(5,:) = [1 0 1];
colormap(6,:) = [1 1 0];
colormap(7,:) = [0 0 0];
colormap(8,:) = [1 1 1];
colormap(9,:) = [0 0.4470 0.7410];
colormap(10,:) = [0.8500 0.3250 0.0980];
colormap(11,:) = [0.9290 0.6940 0.1250];
colormap(12,:) = [0.4940 0.1840 0.5560];
colormap(13,:) = [0.4660 0.6740 0.1880];
colormap(14,:) = [0.3010 0.7450 0.9330];
colormap(15,:) = [0.6350 0.0780 0.1840];

% Load Test Image
% imagePath = 'test_image.tif';
% sar_img = imread(imagePath);
% sar_img = im2double(sar_img);

fp = fopen('napoli.raw','rb'); 
X  = fread(fp,[256,256],'float')'; 
fclose(fp);

fp = fopen('napoli_noisy.raw','rb'); 
Z  = fread(fp,[256,256],'float')'; 
fclose(fp);

sar_img = Z;

% SSIM parameters
K = [0.05 0.05];
window = ones(9);
% window = fspecial('gaussian', 11, 1.5);
L = 100;
L = 255;
[mssim_noisy, ssim_map_noisy] = ssim_index(X, sar_img, K, window, L);

%% DESPECKLING

mean_window = 3;
tic
filtered_image(1,:,:) = localMean(sar_img, mean_window, 'off');
duration_desp(1) = toc;

median_window = [3 3];
tic
filtered_image(2,:,:) = medfilt2(sar_img, median_window);
duration_desp(2) = toc;

tic
filtered_image(3,:,:) = FrostFilter(sar_img, 4, 2.2);
duration_desp(3) = toc;

tic
filtered_image(4,:,:) = modifiedFrost(sar_img, 4, 5);
duration_desp(4) = toc;

tic
filtered_image(5,:,:) = LeeFilter(sar_img, 4, 4.5, 'r');
duration_desp(5) = toc;

tic
filtered_image(6,:,:) = double(kuan(sar_img, 1));
duration_desp(6) = toc;

tic
filtered_image(7,:,:) = double(kuwahara(sar_img,5,1));
duration_desp(7) = toc;

lambda = 25;
targetEpsilon = 0.1;
targetNorm = 0.5;
maxIter = 5;
solverTolerance = 1e-2;
solverMaxIter = 10000;
tic
filtered_image(8,:,:) = double(SDD(sar_img, lambda, targetEpsilon, targetNorm, maxIter, solverTolerance, solverMaxIter));
duration_desp(8) = toc;

lambda = 50;
epsilon = 1e-5;
adaptiveEpsilon = 1;
alpha = 0.5;
gamma = 1.0;
maximumIteration = 5;
solverMaximumIteration = 1000;
solverTolerance = 1e-2;
%preconditioner = 'JACOBI';
preconditioner = 'IC';
tic
[filtered_image(9,:,:), errList, errListCR, totalCGIteration] = SDD_QL(sar_img, lambda, epsilon, adaptiveEpsilon, alpha, gamma, maximumIteration, solverMaximumIteration, solverTolerance, preconditioner);
duration_desp(9) = toc;

addpath '...\SARBM3D_win64'
L_SARBM3D = 1;
tic
filtered_image(10,:,:) = SARBM3D(sar_img,L_SARBM3D);
duration_desp(10) = toc;

TMAX = .75;
tic
filtered_image(11,:,:) = twodncdf(sar_img, TMAX);
duration_desp(11) = toc;

tic
filtered_image(12,:,:) = SRAD(uint8(sar_img),300,0.05);
duration_desp(12) = toc;

addpath '...\IDF'
tic
filtered_image(13,:,:) = addnoisefilter(sar_img,[3 3], 3);
duration_desp(13) = toc;

addpath '...\IDF'
tic
filtered_image(14,:,:) = aditiv (uint8(sar_img), [3 3], 3);
duration_desp(14) = toc;

addpath '...\IDF'
tic
filtered_image(15,:,:) = am_noisenew(uint8(sar_img), [3 3], 4); 
duration_desp(15) = toc;

tic
filtered_image(16,:,:) = sar_img;
duration_desp(16) = toc;

%% Parameters
for i = 1:1:16
    Y = squeeze(filtered_image(i,:,:));
    [mssim_vals(i), ssim_maps(i,:,:)] = ssim_index(X, Y, K, window, L);
    [mdssim_vals(i), dssim_maps(i,:,:)] = dssim_index(X, Y, K, window, L);
    [psnr_vals(i), psnr_maps(i,:,:)] = psnr_index(X, Y);
    [mse_vals(i), mse_maps(i,:,:)] = mse_index(X, Y);
    [enl_vals(i), enl_maps(i,:,:)] = enl_index(X, Y, 9); 
    [ent_vals(i), ent_maps(i,:,:)] = entropy_index(X, Y, 9); 
    [epi_vals(i), epi_maps(i,:,:)] = epi_index(X, Y); 
    [q_vals(i), q_maps(i,:,:)] = qindex_sar(X, Y); 
    [uqi_vals(i), uqi_maps(i,:,:)] = uqi_index(X, Y); 
    [mEQP_val(i)] = eqp_cal(Y);
end

%% Write in a .txt file

mSSIMs = [round(mssim_vals,3)]';
mDSSIMs = [round(mssim_vals,3)]';
PSNRs = [round(psnr_vals,3)]';
MSEs = [round(mse_vals,1)]';
ENLs = [round(enl_vals,1)]';
ENTs = [round(ent_vals,3)]';
EPIs = [round(epi_vals,3)]';
Q_index = [round(q_vals,3)]';
UQIs = [round(uqi_vals,3)]';
mEQPs = [round(mEQP_val,3)]';
Durations = [round(duration_desp,3)]';

rowNames = {'Mean';'Median';'Frost';'M.Frost';'Lee';'Kuan';'Kuwahar';'SDD';'SDD-QL';'SARBM3D';'NCDF';'SRAD';'adnoise';'aditiv';'amnoise';'noisy'};
T2 = table(mSSIMs,mDSSIMs,PSNRs,MSEs,ENLs,ENTs,EPIs,Q_index,UQIs,mEQPs,Durations,'RowNames',rowNames);
writetable(T2,'tabledata.txt','Delimiter','\t','WriteRowNames',true);
type tabledata.txt