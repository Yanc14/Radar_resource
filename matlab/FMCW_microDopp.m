clear all; 
close all;
clc;

GHz=1e+9;
MHz=1e+6;
us=1e-6;
c = physconst('LightSpeed'); %speed of light
load('FMCW-MIMO-Radar-Simulation/20220817_191142_00.mat');
%%===========================雷达参数设置===========================%%
BW = 4*GHz;                     %调频带宽
start_freq = 77*GHz;            %起始频率
end_freq = start_freq + BW;     %结束频率
fc = (start_freq + end_freq)/2; %载波频率
slope = BW/T;                   %调频系数
lambda = c/fc;
numChirps = 256;                %chirp数
frameDuration = 40e-3;          %帧周期
T = frameDuration/numChirps;    %PRI脉冲重复周期T_r
PRF = 1/T;                      %脉冲重复频率

%采样参数设置
numADC = 256; % 采样点数
F = numADC/T; %采样频率
dt = 1/F;     %采样间隔
t_onePulse = 0:dt:dt*numADC-dt;
%%=================================雷达性能设置========================%%
numTX = 1;
numRX = 1;


Vmax = lambda/(T*4);  %最大不模糊速度
Rmax = F*c/(2*slope); %最大不模糊距离
dR = c/(2*BW);        %距离分辨率

DFmax = 1/2*PRF;      %最大不模糊多普勒频率
Rmax2 = c/2/PRF; % lecture 2.3

d_rx = lambda/2; % dist. between rxs
d_tx = 4*d_rx; % dist. between txs

N_range = numADC; % length of range FFT
N_azimuth = numTX*numRX;
R = 0:dR:Rmax-dR; % range axis
ang_ax = -90:90; % angle axis

%%===================================空间分布=======================%%
%%1、雷达空间分布
tx_loc = cell(1,numTX);
for i = 1:numTX
   tx_loc{i} = [(i-1)*d_tx 0 0];%x轴等间隔排列
end

rx_loc = cell(1,numRX);
for i = 1:numRX
   rx_loc{i} = [tx_loc{numTX}(1)+d_tx+(i-1)*d_rx 0 0];%x轴等间隔排列
end

%%2、目标空间分布，
%skel_hist为3D coordinates of 关节 captured by Kinect v2 device
fps_skel = 30;%帧率
num_tar = size(skel_hist,1);%目标数
durationx = size(skel_hist,3)/fps_skel;%处理skel数据需要的时间
%同步处理skel和chirp的时间->需要的chirp数量，即总帧数
sumChirps = floor(numChirps*(durationx/frameDuration));

tar_loc = zeros(num_tar, size(skel_hist,2), sumChirps*numADC);
for t = 1:num_tar
    for i = 1:size(skel_hist,2)
        %在1:size(skel_hist,3)之间插入linspace(1,size(skel_hist,3),numChirps*numADC)x
        tar_loc(t,i,:) = interp1(1:size(skel_hist,3), squeeze(skel_hist(t,i,:)), linspace(1,size(skel_hist,3),sumChirps*numADC));      
    end
end

temp = linspace(2.7, 0.5, sumChirps*numADC);
v_avg = (max(tar_loc(1,2,:)) - min(tar_loc(1,2,:))) * sqrt(3) / durationx;

%%3、计算时延
delays_targets = cell(numTX,numRX,num_tar);
for t = 1:num_tar
    for i = 1:numTX
        for j = 1:numRX
            delays_targets{i,j,t} = (vecnorm(squeeze(tar_loc(t,:,:))-repmat(rx_loc{j},size(tar_loc,3),1).',2,1)+vecnorm(squeeze(tar_loc(t,:,:))-repmat(tx_loc{i},size(tar_loc,3),1).',2,1))/c; 
        end
    end
end
%%4、混频
phase = @(tx,fx) 2*pi*(fx.*tx+slope/2*tx.^2); 
phase_t = phase(t_onePulse,fc);%发射信号的瞬时相位

mixed = zeros(numTX,numRX,size(tar_loc,3));
for i = 1:numTX
    for j = 1:numRX
        disp(['处理关节数: ' num2str(j) '/' num2str(numRX)]);
        for t = 1:num_tar
            disp([int2str(t) '/' int2str(num_tar)]);
            for k = 1:sumChirps
                phase_tar = phase(t_onePulse-delays_targets{i,j,t}(k*numADC),fc); % received
                signal_tar((k-1)*numADC+1:k*numADC) = exp(1j*(phase_t - phase_tar));
            end
            %混合目标回波
            mixed(i,j,:) = squeeze(mixed(i,j,:)) + signal_tar.';
        end
    end
end
%%==========================距离-速度 2D-FFT========================%%
% RDC = reshape(cat(3,mixed{:}),numADC,numChirps*numCPI,numRX*numTX); 
% radar data cube
RDC = reshape(mixed,numADC,sumChirps,numRX*numTX);
numCPI = floor(sumChirps/numChirps);
RDMs = zeros(numADC,numChirps,numTX*numRX,numCPI);
for i = 1:numCPI
    %截取一帧（256chirps）做range-dop FFT
    RD_frame = RDC(:,(i-1)*numChirps+1:i*numChirps,:);
    RDMs(:,:,:,i) = fftshift(fft2(RD_frame,[],[]),2);
end

figure
colormap(jet(256))
for f = 1:numCPI
    imagesc([-Vmax Vmax], [0 Rmax], 20*log10(abs(RDMs(:,:,1,f))/max(max(abs(RDMs(:,:,1,f))))));
    clim = get(gca,'clim');
    caxis([clim(1)/2 0])
    xlabel('Velocity (m/s)');
    ylabel('Range (m)');
    title(['Range-Doppler Map, Frame: ' int2str(f) '/' int2str(numCPI)]);
    drawnow;
    F2(f) = getframe(gcf); % gcf returns the current figure handle
    pause(frameDuration);
end

%保存帧播放文件
writerObj = VideoWriter('FMCW-MIMO-Radar-Simulation/test.avi');
writerObj.FrameRate = floor(1/frameDuration);
open(writerObj);
for i=1:length(F2)
    frame = F2(i);
    writeVideo(writerObj, frame);
end
close(writerObj);
%%=================================微多普勒谱图==============================%%
rBin = 1:256;
nfft = 2^12;
window = 2^7;
noverlap = 100;
shift = window - noverlap;%窗口移动步长，窗口必须大于移动步长
sx = myspecgramnew(sum(RDC(rBin,:,:),1),window,nfft,shift); %对慢时域进行时频分析
sx2 = abs(flipud(fftshift(sx,1)));%上下翻转sx

timeAxis = [1:numCPI]*frameDuration; % Time
freqAxis = linspace(-PRF/2,PRF/2,nfft); % Frequency Axis

fig = figure('visible','on');
colormap(jet(256));
imagesc(timeAxis,[-PRF/2 PRF/2],20*log10(abs(sx2/max(sx2(:)))));
title('micro-Doppler Spectrogram');
xlabel('Time (sec)');ylabel('Frequency (Hz)');
caxis([-45 0]);
set(gca, 'YDir','normal');
set(gcf,'color','w');
