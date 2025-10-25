%% Pixel-Level Fusion of SAR Image Despeckling Algorithms
% 
% Dr. Emrah Onat
% 25.10.2025


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
% XX - SSIM-Based FUSION Filter
% XX - PSNR-Based FUSION Filter
% XX - ENT-Based FUSION Filter
% XX - UQI-Based FUSION Filter


% Metrics
% 01 - SSIM
% 02 - MSSIM
% 03 - PSNR
% 04 - MSE
% 05 - ENL
% 06 - Entrophy
% 07 - EPi
% 08 - Q index
% 09 - UQI
% 20 - EQP


Entrophy

clear all;
close all; clearvars; clc;

colormap(1,:) = [1 0 0]; 
colormap(2,:) = [0 1 0];
colormap(3,:) = [0 0 1];
colormap(4,:) = [0 1 1];
colormap(5,:) = [1 0 1];
colormap(6,:) = [1 1 0];
colormap(7,:) = [0.2950 0.1980 0.2240];
colormap(8,:) = [1 1 1];
colormap(9,:) = [0 0.4470 0.7410];
colormap(10,:) = [0.8500 0.3250 0.0980];
colormap(11,:) = [0.9290 0.6940 0.1250];
colormap(12,:) = [0.4940 0.1840 0.5560];
colormap(13,:) = [0.4660 0.6740 0.1880];
colormap(14,:) = [0.3010 0.7450 0.9330];
colormap(15,:) = [0.6350 0.0780 0.1840];
colormap(16,:) = [0.6150 0.0880 0.1540];
colormap(17,:) = [0.5950 0.0980 0.1240];
colormap(30,:) = [0 0 0];
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
windowsize = 5;
WS = (windowsize-1)/2;

% Mean Filter
tic
filtered_image(1,:,:) = localMean(sar_img, windowsize, 'off');
duration_desp(1) = toc;

% Median Filter
median_window = [windowsize windowsize];
tic
filtered_image(2,:,:) = medfilt2(sar_img, median_window);
duration_desp(2) = toc;

% Frost Filter
tic
filtered_image(3,:,:) = FrostFilter(sar_img, WS, 2.2);
duration_desp(3) = toc;

% Modified Frost Filter
tic
filtered_image(4,:,:) = modifiedFrost(sar_img, WS, 5);
duration_desp(4) = toc;

% Lee Filter
tic
filtered_image(5,:,:) = LeeFilter(sar_img, WS, 4.5, 'r');
duration_desp(5) = toc;

% Kuan Filter
tic
filtered_image(6,:,:) = double(kuan(sar_img, 1));
duration_desp(6) = toc;

% Kuwahara Filter
tic
filtered_image(7,:,:) = double(kuwahara(sar_img,windowsize,1));
duration_desp(7) = toc;

% SDD Filter
lambda = 25;
targetEpsilon = 0.1;
targetNorm = 0.5;
maxIter = 5;
solverTolerance = 1e-2;
solverMaxIter = 10000;
tic
filtered_image(8,:,:) = double(SDD(sar_img, lambda, targetEpsilon, targetNorm, maxIter, solverTolerance, solverMaxIter));
duration_desp(8) = toc;

% SDD-QL Filter
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

% SARBM3D Filter
addpath 'C:\Users\emrah\OneDrive\Documents\Academia\Despeckling\MATLAB_Despeckling\SARBM3D_win64'
L_SARBM3D = 1;
tic
filtered_image(10,:,:) = SARBM3D(sar_img,L_SARBM3D);
duration_desp(10) = toc;

% 2DNCDF Filter
TMAX = .75;
tic
filtered_image(11,:,:) = twodncdf(sar_img, TMAX);
duration_desp(11) = toc;

% SRAD Filter
tic
filtered_image(12,:,:) = SRAD(uint8(sar_img),300,0.05);
duration_desp(12) = toc;

% Add Noise Filter
addpath 'C:\Users\emrah\OneDrive\Documents\Academia\Despeckling\MATLAB_Despeckling\IDF'
tic
filtered_image(13,:,:) = addnoisefilter(sar_img,[3 3], 3);
duration_desp(13) = toc;

% Aditiv Filter
addpath 'C:\Users\emrah\OneDrive\Documents\Academia\Despeckling\MATLAB_Despeckling\IDF'
tic
filtered_image(14,:,:) = aditiv (uint8(sar_img), [3 3], 3);
duration_desp(14) = toc;

% AM Noise Filter
addpath 'C:\Users\emrah\OneDrive\Documents\Academia\Despeckling\MATLAB_Despeckling\IDF'
tic
filtered_image(15,:,:) = am_noisenew(uint8(sar_img), [3 3], 4); 
duration_desp(15) = toc;

% Noisy Image
tic
filtered_image(16,:,:) = sar_img;
duration_desp(16) = toc;

number_of_filter = 12;

%% Figures

[output] = fig_filters(X, filtered_image, colormap);

%% Parameters
for i = 1:1:number_of_filter
    Y = squeeze(filtered_image(i,:,:));
    X1 = pad_with_zeros(X, 9);
    Y1 = pad_with_zeros(Y, 9);
    [mssim_vals(i), ssim_maps(i,:,:)] = ssim_index(X1, Y1, K, window, L);   % 1
    [mdssim_vals(i), dssim_maps(i,:,:)] = dssim_index(X1, Y1, K, window, L);% 2
    [psnr_vals(i), psnr_maps(i,:,:)] = psnr_index(X, Y);                    % 3
    [mse_vals(i), mse_maps(i,:,:)] = mse_index(X, Y);                       % 4
    [enl_vals(i), enl_maps(i,:,:)] = enl_index(X1, Y1, 9);                  % 5
    [ent_vals(i), ent_maps(i,:,:)] = entropy_index(X, Y, 9);                % 6
    [epi_vals(i), epi_maps(i,:,:)] = epi_index(X, Y);                       % 7
    [q_vals(i), q_maps(i,:,:)] = qindex_sar(X, Y);                          % 8
    [uqi_vals(i), uqi_maps(i,:,:)] = uqi_index(X, Y);                       % 9
    [mEQP_val(i)] = eqp_cal(Y);
end

%% FUSION
tic
[fused_despekled_images, type_desp_maps, fused_color_maps fused_param_map_values] = fusions_despekling(filtered_image, ssim_maps, psnr_maps, mse_maps, enl_maps, ent_maps, epi_maps, q_maps, uqi_maps, colormap);
duration_desp(number_of_filter+1) = toc;
duration_desp(number_of_filter+2) = toc;
duration_desp(number_of_filter+3) = toc;
duration_desp(number_of_filter+4) = toc;
duration_desp(number_of_filter+5) = toc;
duration_desp(number_of_filter+6) = toc;
duration_desp(number_of_filter+7) = toc;

% tic
% [fused_image_Pareto, type_map_Pareto, fused_color_maps_Pareto, score_map_Pareto] = pareto_fusion(filtered_image, ssim_maps, psnr_maps, mse_maps, ent_maps, uqi_maps, colormap);
% duration_desp(number_of_filter+6) = toc;
% 
% tic
% [fused_image_learned, type_map_learned, fused_color_maps_learned, score_map_learned, net_learned] = learned_fusion(filtered_image, ssim_maps, psnr_maps, mse_maps, ent_maps, uqi_maps, colormap, X);
% duration_desp(number_of_filter+7) = toc;

%% Calculations - Histograms - Figures

% SSIM
fused_despekled_image_SSIM = squeeze(fused_despekled_images(1,:,:));
type_desp_image_SSIM = squeeze(type_desp_maps(1,:,:));
fused_color_maps_SSIM = squeeze(fused_color_maps(1,:,:,:));
fused_param_map_values_SSIM = squeeze(fused_param_map_values(1,:,:));

[mean_param_SSIM param_SSIM] = calc_hist_fig(X, fused_despekled_image_SSIM, type_desp_image_SSIM, fused_color_maps_SSIM, fused_param_map_values_SSIM, colormap, number_of_filter, K, window, L);

% PSNR
fused_despekled_image_PSNR = squeeze(fused_despekled_images(3,:,:));
type_desp_image_PSNR = squeeze(type_desp_maps(3,:,:));
fused_color_maps_PSNR = squeeze(fused_color_maps(3,:,:,:));
fused_param_map_values_PSNR = squeeze(fused_param_map_values(3,:,:));

[mean_param_PSNR param_PSNR] = calc_hist_fig(X, fused_despekled_image_PSNR, type_desp_image_PSNR, fused_color_maps_PSNR, fused_param_map_values_PSNR, colormap, number_of_filter, K, window, L);

% MSE
fused_despekled_image_MSE = squeeze(fused_despekled_images(4,:,:));
type_desp_image_MSE = squeeze(type_desp_maps(4,:,:));
fused_color_maps_MSE = squeeze(fused_color_maps(4,:,:,:));
fused_param_map_values_MSE = squeeze(fused_param_map_values(4,:,:));

[mean_param_MSE param_MSE] = calc_hist_fig(X, fused_despekled_image_MSE, type_desp_image_MSE, fused_color_maps_MSE, fused_param_map_values_MSE, colormap, number_of_filter, K, window, L);

% ENL
fused_despekled_image_ENL = squeeze(fused_despekled_images(5,:,:));
type_desp_image_ENL = squeeze(type_desp_maps(5,:,:));
fused_color_maps_ENL = squeeze(fused_color_maps(5,:,:,:));
fused_param_map_values_ENL = squeeze(fused_param_map_values(5,:,:));

[mean_param_ENL param_ENL] = calc_hist_fig(X, fused_despekled_image_ENL, type_desp_image_ENL, fused_color_maps_ENL, fused_param_map_values_ENL, colormap, number_of_filter, K, window, L);

% ENT
fused_despekled_image_ENT = squeeze(fused_despekled_images(6,:,:));
type_desp_image_ENT = squeeze(type_desp_maps(6,:,:));
fused_color_maps_ENT = squeeze(fused_color_maps(6,:,:,:));
fused_param_map_values_ENT = squeeze(fused_param_map_values(6,:,:));

[mean_param_ENT param_ENT] = calc_hist_fig(X, fused_despekled_image_ENT, type_desp_image_ENT, fused_color_maps_ENT, fused_param_map_values_ENT, colormap, number_of_filter, K, window, L);

% EPI
fused_despekled_image_EPI = squeeze(fused_despekled_images(7,:,:));
type_desp_image_EPI = squeeze(type_desp_maps(7,:,:));
fused_color_maps_EPI = squeeze(fused_color_maps(7,:,:,:));
fused_param_map_values_EPI = squeeze(fused_param_map_values(7,:,:));

[mean_param_EPI param_EPI] = calc_hist_fig(X, fused_despekled_image_EPI, type_desp_image_EPI, fused_color_maps_EPI, fused_param_map_values_EPI, colormap, number_of_filter, K, window, L);

% ENT
fused_despekled_image_Q = squeeze(fused_despekled_images(8,:,:));
type_desp_image_Q = squeeze(type_desp_maps(8,:,:));
fused_color_maps_Q = squeeze(fused_color_maps(8,:,:,:));
fused_param_map_values_Q = squeeze(fused_param_map_values(8,:,:));

[mean_param_Q param_Q] = calc_hist_fig(X, fused_despekled_image_Q, type_desp_image_Q, fused_color_maps_Q, fused_param_map_values_Q, colormap, number_of_filter, K, window, L);

% UQI
fused_despekled_image_UQI = squeeze(fused_despekled_images(9,:,:));
type_desp_image_UQI = squeeze(type_desp_maps(9,:,:));
fused_color_maps_UQI = squeeze(fused_color_maps(9,:,:,:));
fused_param_map_values_UQI = squeeze(fused_param_map_values(9,:,:));

[mean_param_UQI param_UQI] = calc_hist_fig(X, fused_despekled_image_UQI, type_desp_image_UQI, fused_color_maps_UQI, fused_param_map_values_UQI, colormap, number_of_filter, K, window, L);

% % Pareto
% fused_despekled_image_Pareto = squeeze(fused_image_Pareto(1,:,:));
% type_desp_image_Pareto = squeeze(type_map_Pareto(1,:,:));
% fused_color_maps_Pareto = fused_color_maps_Pareto;
% fused_param_map_values_Pareto = squeeze(score_map_Pareto(1,:,:));
% 
% [mean_param_Pareto param_Pareto] = calc_hist_fig(X, fused_despekled_image_Pareto, type_desp_image_Pareto, fused_color_maps_Pareto, fused_param_map_values_Pareto, colormap, number_of_filter, K, window, L);
% 
% % Learned
% fused_despekled_image_learned = squeeze(fused_image_learned(1,:,:));
% type_desp_image_learned = squeeze(type_map_learned(1,:,:));
% fused_color_maps_learned = fused_color_maps_learned;
% fused_param_map_values_learned = squeeze(score_map_learned(1,:,:));
% 
% [mean_param_learned param_learned] = calc_hist_fig(X, fused_despekled_image_learned, type_desp_image_learned, fused_color_maps_learned, fused_param_map_values_learned, colormap, number_of_filter, K, window, L);


% Calculated Values
mssim_vals  = [mssim_vals param_SSIM(1) param_PSNR(1) param_ENL(1) param_ENT(1) param_EPI(1) param_Q(1) param_UQI(1)];
mdssim_vals = [mdssim_vals param_SSIM(2) param_PSNR(2) param_ENL(2) param_ENT(2) param_EPI(2) param_Q(2) param_UQI(2)];
psnr_vals   = [psnr_vals param_SSIM(3) param_PSNR(3) param_ENL(3) param_ENT(3) param_EPI(3) param_Q(3) param_UQI(3)];
mse_vals    = [mse_vals param_SSIM(4) param_PSNR(4) param_ENL(4) param_ENT(4) param_EPI(4) param_Q(4) param_UQI(4)];
enl_vals    = [enl_vals param_SSIM(5) param_PSNR(5) param_ENL(5) param_ENT(5) param_EPI(5) param_Q(5) param_UQI(5)];
ent_vals    = [ent_vals param_SSIM(6) param_PSNR(6) param_ENL(6) param_ENT(6) param_EPI(6) param_Q(6) param_UQI(6)];
epi_vals    = [epi_vals param_SSIM(7) param_PSNR(7) param_ENL(7) param_ENT(7) param_EPI(7) param_Q(7) param_UQI(7)];
q_vals      = [q_vals param_SSIM(8) param_PSNR(8) param_ENL(8) param_ENT(8) param_EPI(8) param_Q(8) param_UQI(8)];
uqi_vals    = [uqi_vals param_SSIM(9) param_PSNR(9) param_ENL(9) param_ENT(9) param_EPI(9) param_Q(9) param_UQI(9)];
mEQP_val    = [mEQP_val param_SSIM(10) param_PSNR(10) param_ENL(10) param_ENT(10) param_EPI(10) param_Q(10) param_UQI(10)];

%% Percentages
centers = 0:number_of_filter;

[N,binCenters] = hist(type_desp_image_SSIM(:),centers);
pSSIM = round((100*N./sum(N))',4);

[N,binCenters] = hist(type_desp_image_PSNR(:),centers);
pPSNR = round((100*N./sum(N))',4);

[N,binCenters] = hist(type_desp_image_ENL(:),centers);
pENL = round((100*N./sum(N))',4);

[N,binCenters] = hist(type_desp_image_ENT(:),centers);
pENT = round((100*N./sum(N))',4);

[N,binCenters] = hist(type_desp_image_EPI(:),centers);
pEPI = round((100*N./sum(N))',4);

[N,binCenters] = hist(type_desp_image_Q(:),centers);
pQ = round((100*N./sum(N))',4);

[N,binCenters] = hist(type_desp_image_UQI(:),centers);
pUQI = round((100*N./sum(N))',4);

rowNames = {'Zeros';'Mean';'Median';'Frost';'M.Frost';'Lee';'Kuan';'Kuwahar';'SDD';'SDD-QL';'SARBM3D';'NCDF';'SRAD';'adnoise';'aditiv';'amnoise';'noisy'};
rowNames = rowNames(1:number_of_filter+1,1);
T3 = table(pSSIM,pPSNR,pENL,pENT,pEPI,pQ,pUQI,'RowNames',rowNames);
writetable(T3,'tablehist.txt','Delimiter','\t','WriteRowNames',true);
type tablehist.txt

%% DURATIONS

Durations = [round(duration_desp,3)]';
Durations = Durations(1:number_of_filter+7,1);
centers = 1:number_of_filter+7;
figure;

hBar = bar(centers,Durations,'hist');
set(hBar,'FaceVertexCData',colormap(1:number_of_filter+7,:));
labels = [string('Mean');string('Median');string('Frost');string('M. Frost');string('Lee');string('Kuan');string('Kuwahara');string('SDD');string('SDD-QL');string('SARBM3D');string('NCDF');string('SRAD');string('AddNoise');string('Aditiv');string('AMNoise');string('NoisyImage');string('SSIM');string('PSNR');string('ENL');string('ENT');string('EPI');string('Q');string('UQI')];
labels = [labels(1:number_of_filter,1); labels(17:end,1)];

x = get(hBar,'XData');
y = get(hBar,'YData');
ygap = 0.1;  %// Specify vertical gap between the bar and label
ylimits = get(gca,'YLim');

for i = 1:1:length(x) %// Loop over each bar
    xpos = x(2,i)+0.5;        %// Set x position for the text label
    ypos = y(2,i) + ygap; %// Set y position, including gap
    htext = text(xpos,ypos,labels{i});          %// Add text label
    set(htext,'VerticalAlignment','bottom', 'HorizontalAlignment','center')
end

title('Durations of Algorithms (second)');

% bar(centers,duration_vals);

%% Write in a .txt file

mSSIMs = [round(mssim_vals,3)]';
mDSSIMs = [round(mdssim_vals,3)]';
PSNRs = [round(psnr_vals,3)]';
MSEs = [round(mse_vals,1)]';
ENLs = [round(enl_vals,1)]';
ENTs = [round(ent_vals,3)]';
EPIs = [round(epi_vals,3)]';
Q_index = [round(q_vals,3)]';
UQIs = [round(uqi_vals,3)]';
mEQPs = [round(mEQP_val,3)]';

rowNames = {'Mean';'Median';'Frost';'M.Frost';'Lee';'Kuan';'Kuwahar';'SDD';'SDD-QL';'SARBM3D';'NCDF';'SRAD';'adnoise';'aditiv';'amnoise';'noisy';'SSIM';'PSNR';'ENL';'ENT';'EPI';'Q';'UQI'};
rowNames = [rowNames(1:number_of_filter,1); rowNames(17:end,1)];

T2 = table(mSSIMs,mDSSIMs,PSNRs,MSEs,ENLs,ENTs,EPIs,Q_index,UQIs,mEQPs,Durations,'RowNames',rowNames);
writetable(T2,'tabledata.txt','Delimiter','\t','WriteRowNames',true);
type tabledata.txt

open('tabledata.txt');

open('tablehist.txt');