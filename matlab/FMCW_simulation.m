clear all;
close all;
clc;
GHz=1e+9;
MHz=1e+6;
us=1e-6;
c = physconst('LightSpeed'); %����
%%======================================�״������������==================================%%
BW = 150*MHz;    %��Ƶ����
fc = 77*GHz;     %����Ƶ�� 
T = 10*us;       %PRI�����ظ�ʱ�䣬�˴������Ȧ�=�����ظ�ʱ��T_r
PRF = 1/T;       %�����ظ�Ƶ��
slope = BW/T;    %��Ƶϵ��
numChirps = 256; %֡chirp����
%CPI Coherent Processing Interval 
%��ɴ�������CPI��ͨ���й̶���PRI���״�Ƶ�ʣ����ҷ��䲨����ͬ���Ա���ж����ղ�����
numCPI = 10;     %256chirp����һ֡����10֡


%������������
numADC = 256; % # of adc samples
F = numADC/T;   % ����Ƶ��
dt = 1/F;       % �������
lambda = c/fc;  
N = numChirps*numADC*numCPI; % �ܲ�������
t = linspace(0,T*numChirps*numCPI,N); %��ʱ��ʱ����Ϊ10us/256��
t_onePulse = 0:dt:dt*numADC-dt;%��ʱ�򣬼��Ϊ10us/256
%%=================================�״���������==============================================%%
numTX = 1;
numRX = 8;

%���ģ���ٶ�
Vmax = lambda/(T*4); % Max Unamb velocity m/s
%�ٶȷֱ���
dV = lambda/(2*numChirps*T); % velocity resol, lambda/(2*framePeriod)
%���ģ�����루����Ƶ�ʣ�
Rmax = F*c/(2*slope); % TI's MIMO Radar doc
%����ֱ���
dR = c/(2*BW); % range resol

%Rmax2 = c/2/PRF; % lecture 2.3
%DFmax = 1/2*PRF; % = Vmax/(c/fc/2); % Max Unamb Dopp Freq

%�������߼��Ϊ���߳���lambda/2 (meter)
d_rx = lambda/2; % dist. between rxs
%�������߼��
d_tx = 4*d_rx; % dist. between txs

%%====================================�ռ�ֲ�===============================================%%
%���״�Ϊ���ģ��״ﲨ��Ϊ������������άƽ������ϵ
%�״�ռ�ֲ�

figure;
tx_loc = cell(1,numTX);
for i = 1:numTX
   tx_loc{i} = [(i-1)*d_tx 0 0];
   scatter3(tx_loc{i}(1),tx_loc{i}(2),tx_loc{i}(3),'blue','filled');
   hold on;
end
rx_loc = cell(1,numRX);
for i = 1:numRX
   rx_loc{i} = [tx_loc{numTX}(1)+d_tx+(i-1)*d_rx 0 0];
   scatter3(rx_loc{i}(1),rx_loc{i}(2),rx_loc{i}(3),'red','filled');
end
title("�״�ռ�ֲ�");grid on;
hold off;
xlabel("x(m)");ylabel("y(m)"),zlabel("z(m)");
legend("��������","��������");


%Ŀ��Ŀռ�ֲ�[x(t),y(t),z(t)]
%Ŀ��1
r1_radial = 50;     %�뾶��m��
v1_radial = 10;     %�����ٶȣ�m/s��
tar1_angle = -15;   %�Ƕȣ����ȣ�
r1_y = cosd(tar1_angle)*r1_radial;
r1_x = sind(tar1_angle)*r1_radial;
v1_y = cosd(tar1_angle)*v1_radial;
v1_x = sind(tar1_angle)*v1_radial;
r1 = [r1_x r1_y 0];
%�ռ�����
tar1_loc = zeros(length(t),3);
tar1_loc(:,1) = r1(1) + v1_x*t;
tar1_loc(:,2) = r1(2) + v1_y*t;
%Ŀ��2
r2_radial = 100; %����
v2_radial = -15; %�ٶ�
tar2_angle = 10; %�Ƕ�
r2_y = cosd(tar2_angle)*r2_radial;
r2_x = sind(tar2_angle)*r2_radial;
v2_y = cosd(tar2_angle)*v2_radial;
v2_x = sind(tar2_angle)*v2_radial;
r2 = [r2_x r2_y 0];
%�ռ�����
tar2_loc = zeros(length(t),3);
tar2_loc(:,1) = r2(1) + v2_x*t;
tar2_loc(:,2) = r2(2) + v2_y*t;

%ʱ�Ӳ���
delays_tar1 = cell(numTX,numRX);
delays_tar2 = cell(numTX,numRX);
r1_at_t = cell(numTX,numRX);
r2_at_t = cell(numTX,numRX);
%��λ����
tar1_angles = cell(numTX,numRX);
tar2_angles = cell(numTX,numRX);
%�ٶȲ���
tar1_velocities = cell(numTX,numRX);
tar2_velocities = cell(numTX,numRX);

%����ʱ��t_d=(R_tx+R_rx)/c
for i = 1:numTX
    for j = 1:numRX
        delays_tar1{i,j} = (vecnorm(tar1_loc-repmat(rx_loc{j},N,1),2,2)+vecnorm(tar1_loc-repmat(tx_loc{i},N,1),2,2))/c; 
        delays_tar2{i,j} = (vecnorm(tar2_loc-repmat(rx_loc{j},N,1),2,2)+vecnorm(tar2_loc-repmat(tx_loc{i},N,1),2,2))/c;
    end
end

%����Ƶ��=��ֹ��ƵƵ��+������Ƶ��
fr1 = 2*r1(2)*slope/c; 
fr2 = 2*r2(2)*slope/c;

fd1 = 2*v1_radial*fc/c; % doppler freq
fd2 = 2*v2_radial*fc/c;

f_if1 = fr1 + fd1; % beat or IF freq
f_if2 = fr2 + fd2;


%��Ƶ
phase = @(tx,fx) 2*pi*(fx.*tx+slope/2*tx.^2); % �����ź�˲ʱ��λ
mixed = cell(numTX,numRX);
for i = 1:numTX
    for j = 1:numRX
        disp(['Processing Channel: ' num2str(j) '/' num2str(numRX)]);
        for k = 1:numChirps*numCPI
            phase_t = phase(t_onePulse,fc);%�����Ƶ�ź���λ
            phase_1 = phase(t_onePulse-delays_tar1{i,j}(k*numADC),fc); % Ŀ��1�����ź���λ
            phase_2 = phase(t_onePulse-delays_tar2{i,j}(k*numADC),fc); % Ŀ��2�����ź���λ
            
            signal_t((k-1)*numADC+1:k*numADC) = exp(1j*phase_t);
            signal_1((k-1)*numADC+1:k*numADC) = exp(1j*(phase_t - phase_1));
            signal_2((k-1)*numADC+1:k*numADC) = exp(1j*(phase_t - phase_2));
        end
        mixed{i,j} = signal_1 + signal_2;%ͬʱ���ն��Ŀ��
    end
end

figure
subplot(4,1,1)
p1 = plot(t, real(signal_t));
title('�����ź�');
xlim([0 0.1e-4]);xlabel('ʱ��t/(10us)');ylabel('tx(t)');
subplot(4,1,2)
p2 = plot(t, real(signal_1));
title('Ŀ��1�ز��ź�');
xlim([0 0.1e-4]);xlabel('ʱ��t/(10us)');ylabel('rx1_1(t)');
subplot(4,1,3)
p3 = plot(t, real(signal_2));
title('Ŀ��1�ز��ź�');
xlim([0 0.1e-4]);xlabel('ʱ��t/(10us)');ylabel('rx1_2(t)');
subplot(4,1,4)
p4 = plot(t,real(mixed{i,j}));
title('�������߽����ź�');
xlim([0 0.1e-4]);xlabel('ʱ��t/(10us)');ylabel('rx1_1(t)+rx1_2(t)');
%%=====================================����-�ٶȶ�άFFT�任===========================%%
%N_range��FFT�任��N_dopp��FFT�任��N_azimuth��λ��
N_Dopp = numChirps; % length of doppler FFT
N_range = numADC; % length of range FFT
N_azimuth = numTX*numRX;
R = 0:dR:Rmax-dR; % range axis
V = linspace(-Vmax, Vmax, numChirps); % Velocity axis
ang_ax = -90:90; % angle axis

RDC = reshape(cat(3,mixed{:}),numADC,numChirps*numCPI,numRX*numTX); % radar data cube,�ںϽ��յĶ���ź�
RDMs = zeros(numADC,numChirps,numTX*numRX,numCPI);
for i = 1:numCPI
    RD_frame = RDC(:,(i-1)*numChirps+1:i*numChirps,:);%numADC*frame*numRx��rxͨ��
    RDMs(:,:,:,i) = fftshift(fft2(RD_frame,N_range,N_Dopp),2);%��ǰ����ά��[numchirp,numADC]��FFT
end

figure
%��256*256������Ԫ��ӳ��Ϊһ������ͬһλ�õ�һ������
imagesc(V,R,20*log10(abs(RDMs(:,:,1,1))/max(max(abs(RDMs(:,:,1,1))))));
colormap(gca,jet(256));
clim = get(gca,'clim');
caxis([clim(1)/2 0]);
title("����-�ٶ�FFT");
xlabel('Velocity (m/s)');ylabel('Range (m)');
%%==========================================CA-CFAR==================================%%
numGuard = 2; % ������Ԫ
numTrain = numGuard*2; % ѵ����Ԫ
P_fa = 1e-5; % �龯���� 
SNR_OFFSET = -5; % dB
RDM_dB = 10*log10(abs(RDMs(:,:,1,1))/max(max(abs(RDMs(:,:,1,1)))));

[RDM_mask, cfar_ranges, cfar_dopps, K] = ca_cfar(RDM_dB, numGuard, numTrain, P_fa, SNR_OFFSET);%����ȥ����֡�źţ�Ŀ��׼ȷĿ��[cfar_ranges, cfar_dopps]��Ŀ����K

figure
h=imagesc(V,R,RDM_mask);
title('CA-CFAR');
xlabel('Velocity (m/s)');ylabel('Range (m)')
%%====================================����-��λ�Ƕ�άFFT=====================================%%
rangeFFT = fft(RDC(:,1:numChirps,:),N_range);%��һ��chirp����numADC��FFT
angleFFT = fftshift(fft(rangeFFT,length(ang_ax),3),3);%��һ��chirp����rx��FFT
range_az = squeeze(sum(angleFFT,2)/numChirps); % ����-��λFFTͼ

figure
imagesc(ang_ax,R,20*log10(abs(range_az)./max(abs(range_az(:))))); 
colormap(jet);
set(gca,'clim', [-35, 0]);
title('����-��λ��FFT');
xlabel('Azimuth Angle');ylabel('Range (m)');

%��CA-CFAR���Ŀ��
doas = zeros(K,181); 
figure
hold on; grid on;
for i = 1:K
    doas(i,:) = fftshift(fft(rangeFFT(cfar_ranges(i),cfar_dopps(i),:),181));
    plot(ang_ax,10*log10(abs(doas(i,:))));
end
xlabel('Azimuth Angle');ylabel('dB');
%%======================================AOA ���� - MUSIC�㷨===============================%%
d = 0.5;
M = numCPI; % 

%���������������
%����Ŀ������ڵ��������ĵ��еĽǶ���
for k=1:length(ang_ax)
        a1(:,k)=exp(-1i*2*pi*(d*(0:numTX*numRX-1)'*sin(ang_ax(k).'*pi/180)));
end
 
for i = 1:K
    Rxx = zeros(numTX*numRX,numTX*numRX);
    for m = 1:M
       %��������źŵ�����ؾ��󣨶�numCPI��ƽ������������[cfar_ranges(i),cfar_dopps(i)]λ�õĲ�ͬrx��������
       A = squeeze(RDMs(cfar_ranges(i),cfar_dopps(i),:,m));
       Rxx = Rxx + 1/M * (A*A');
    end
    %���ź�����ؾ��󣬻����������������
    [Q,D] = eig(Rxx); % Q:��������, D: �Խ�Ԫ��Ϊ����ֵ
    [D, I] = sort(diag(D),'descend');
    Q = Q(:,I); 
    Qs = Q(:,1:2); % ����ź�����������Ϊʲô��1����ʱ�����ز��ز�Ƶ�ʲ�ͬ��������ؾ���rank=2
    Qn = Q(:,3:end); % ������������
    
    %MUSIC�ռ���
    for k=1:length(ang_ax)
        music_spectrum(i,k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
    end
end

figure
hold on 
grid on
title('MUSIC Spectrum')
xlabel('Angle in degrees')
for k = 1:K
    plot(ang_ax,log10(abs(music_spectrum(k,:))));
end
%%============================================���ɵ���=================================%%
[~, I] = max(music_spectrum(2,:));
angle1 = ang_ax(I);
[~, I] = max(music_spectrum(1,:));
angle2 = ang_ax(I);

coor1 = [cfar_ranges(2)*sind(angle1) cfar_ranges(2)*cosd(angle1) 0];
coor2 = [cfar_ranges(1)*sind(angle2) cfar_ranges(1)*cosd(angle2) 0];

figure
mesh(-100:1:100,-100:1:100,zeros(201));
hold on;
scatter3(coor1(1),coor1(2),coor1(3),100,'m','filled','linewidth',9)
scatter3(coor2(1),coor2(2),coor2(3),100,'b','filled','linewidth',9)
title('3D Coordinates (Point Cloud) of the targets');
xlabel('Range (m) X');ylabel('Range (m) Y');zlabel('Range (m) Z');
view(60,60);

%%=====================================MUSIC ����-��λ��ͼ=================================%%
rangeFFT = fft(RDC);
range_az_music=zeros(N_range,length(ang_ax));
for i = 1:N_range
    Rxx = zeros(numTX*numRX,numTX*numRX);
    for m = 1:M
       %��ÿһ�У�ADC��,����ѹ����ÿһ������[numRx]��Ϊ����ź�
       A = squeeze(sum(RDMs(i,:,:,m),2));
       %A = squeeze(sum(rangeFFT(i,(m-1)*numChirps+1:m*numChirps,:),2));
       Rxx = Rxx + 1/M * (A*A');
    end
    [Q,D] = eig(Rxx);
    [D, I] = sort(diag(D),'descend');
    Q = Q(:,I);
    index=0;
    for j=1:length(D)-1
        if((D(j)-D(j+1))/D(j)>0.999 && D(j)>=1)
            index=j;
            break
        end
    end
    if(index==0)
        Qn=Q;
    end
    Qn = Q(:,index+1:end);

    for k=1:length(ang_ax)
        range_az_music(i,k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
    end
end

figure
imagesc(ang_ax,R,20*log10(abs(range_az_music)./max(abs(range_az_music(:)))));
colormap(jet);colorbar;
clim = get(gca,'clim');
title('MUSIC Range-Angle Map');
xlabel('Azimuth');ylabel('Range (m)');
%%=====================================Angle Estimation -ѹ����֪==============================%%
numTheta = length(ang_ax); % divide FOV into fine grid
B = a1; % �������
figure
hold on; grid on;
title('Angle Estimation with Compressed Sensing')
xlabel('Azimuth')
ylabel('dB')
for i = 1:K
    A = squeeze(RDMs(cfar_ranges(i),cfar_dopps(i),:,1));%�Ե�һ��CPI�����γ���Ԫ�յ����ź�A
    cvx_begin
        %�ź�A��������������������B��Ϊ���ı�ʾ��A=Bs,��������Ŀ�꣬
        %Ŀ��ֻ���������Ƕȣ�����s��ϡ��ģ�ֻ��Ҫ��ϡ���ʾ���ɣ���ϡ���K�����ϡ�裩
        variable s(numTheta) complex; 
        minimize(norm(s,1))
        norm(A-B*s,2)<=1;
    cvx_end
    cvx_status
    cvx_optval
    plot(ang_ax,10*log10(abs(s)))
end
