clear all; 
close all;
clc;

GHz=1e+9;
MHz=1e+6;
us=1e-6;
c = physconst('LightSpeed'); %speed of light
load('FMCW-MIMO-Radar-Simulation/20220817_191142_00.mat');
%%===========================�״��������===========================%%
BW = 4*GHz;                     %��Ƶ����
start_freq = 77*GHz;            %��ʼƵ��
end_freq = start_freq + BW;     %����Ƶ��
fc = (start_freq + end_freq)/2; %�ز�Ƶ��
slope = BW/T;                   %��Ƶϵ��
lambda = c/fc;
numChirps = 256;                %chirp��
frameDuration = 40e-3;          %֡����
T = frameDuration/numChirps;    %PRI�����ظ�����T_r
PRF = 1/T;                      %�����ظ�Ƶ��

%������������
numADC = 256; % ��������
F = numADC/T; %����Ƶ��
dt = 1/F;     %�������
t_onePulse = 0:dt:dt*numADC-dt;
%%=================================�״���������========================%%
numTX = 1;
numRX = 1;


Vmax = lambda/(T*4);  %���ģ���ٶ�
Rmax = F*c/(2*slope); %���ģ������
dR = c/(2*BW);        %����ֱ���

DFmax = 1/2*PRF;      %���ģ��������Ƶ��
Rmax2 = c/2/PRF; % lecture 2.3

d_rx = lambda/2; % dist. between rxs
d_tx = 4*d_rx; % dist. between txs

N_range = numADC; % length of range FFT
N_azimuth = numTX*numRX;
R = 0:dR:Rmax-dR; % range axis
ang_ax = -90:90; % angle axis

%%===================================�ռ�ֲ�=======================%%
%%1���״�ռ�ֲ�
tx_loc = cell(1,numTX);
for i = 1:numTX
   tx_loc{i} = [(i-1)*d_tx 0 0];%x��ȼ������
end

rx_loc = cell(1,numRX);
for i = 1:numRX
   rx_loc{i} = [tx_loc{numTX}(1)+d_tx+(i-1)*d_rx 0 0];%x��ȼ������
end

%%2��Ŀ��ռ�ֲ���
%skel_histΪ3D coordinates of �ؽ� captured by Kinect v2 device
fps_skel = 30;%֡��
num_tar = size(skel_hist,1);%Ŀ����
durationx = size(skel_hist,3)/fps_skel;%����skel������Ҫ��ʱ��
%ͬ������skel��chirp��ʱ��->��Ҫ��chirp����������֡��
sumChirps = floor(numChirps*(durationx/frameDuration));

tar_loc = zeros(num_tar, size(skel_hist,2), sumChirps*numADC);
for t = 1:num_tar
    for i = 1:size(skel_hist,2)
        %��1:size(skel_hist,3)֮�����linspace(1,size(skel_hist,3),numChirps*numADC)x
        tar_loc(t,i,:) = interp1(1:size(skel_hist,3), squeeze(skel_hist(t,i,:)), linspace(1,size(skel_hist,3),sumChirps*numADC));      
    end
end

temp = linspace(2.7, 0.5, sumChirps*numADC);
v_avg = (max(tar_loc(1,2,:)) - min(tar_loc(1,2,:))) * sqrt(3) / durationx;

%%3������ʱ��
delays_targets = cell(numTX,numRX,num_tar);
for t = 1:num_tar
    for i = 1:numTX
        for j = 1:numRX
            delays_targets{i,j,t} = (vecnorm(squeeze(tar_loc(t,:,:))-repmat(rx_loc{j},size(tar_loc,3),1).',2,1)+vecnorm(squeeze(tar_loc(t,:,:))-repmat(tx_loc{i},size(tar_loc,3),1).',2,1))/c; 
        end
    end
end
%%4����Ƶ
phase = @(tx,fx) 2*pi*(fx.*tx+slope/2*tx.^2); 
phase_t = phase(t_onePulse,fc);%�����źŵ�˲ʱ��λ

mixed = zeros(numTX,numRX,size(tar_loc,3));
for i = 1:numTX
    for j = 1:numRX
        disp(['����ؽ���: ' num2str(j) '/' num2str(numRX)]);
        for t = 1:num_tar
            disp([int2str(t) '/' int2str(num_tar)]);
            for k = 1:sumChirps
                phase_tar = phase(t_onePulse-delays_targets{i,j,t}(k*numADC),fc); % received
                signal_tar((k-1)*numADC+1:k*numADC) = exp(1j*(phase_t - phase_tar));
            end
            %���Ŀ��ز�
            mixed(i,j,:) = squeeze(mixed(i,j,:)) + signal_tar.';
        end
    end
end
%%==========================����-�ٶ� 2D-FFT========================%%
% RDC = reshape(cat(3,mixed{:}),numADC,numChirps*numCPI,numRX*numTX); 
% radar data cube
RDC = reshape(mixed,numADC,sumChirps,numRX*numTX);
numCPI = floor(sumChirps/numChirps);
RDMs = zeros(numADC,numChirps,numTX*numRX,numCPI);
for i = 1:numCPI
    %��ȡһ֡��256chirps����range-dop FFT
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

%����֡�����ļ�
writerObj = VideoWriter('FMCW-MIMO-Radar-Simulation/test.avi');
writerObj.FrameRate = floor(1/frameDuration);
open(writerObj);
for i=1:length(F2)
    frame = F2(i);
    writeVideo(writerObj, frame);
end
close(writerObj);
%%=================================΢��������ͼ==============================%%
rBin = 1:256;
nfft = 2^12;
window = 2^7;
noverlap = 100;
shift = window - noverlap;%�����ƶ����������ڱ�������ƶ�����
sx = myspecgramnew(sum(RDC(rBin,:,:),1),window,nfft,shift); %����ʱ�����ʱƵ����
sx2 = abs(flipud(fftshift(sx,1)));%���·�תsx

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
