clear all;close all;clc;
%%  inisiasi parameter
% VLC Parameter
theta = 50;                      % semi-angle at half power
ml=-log10(2)/log10(cosd(theta)); % Lambertian order of emission
P_LED=7000;                     % Total transmitted power (mW)
Adet=1e-4;                       % detector physical area of a PD (M)

Ib=202e-6;
q = 1.6e-19;                       % Charge of Electron
R = 0.55;                          % Photodetector responsivity
% be = 2e9;                          % Bandwidth
% rl = 30;                           % impedence
% fm = 0.7;                          % 
% T = 300;                           % temperature (kelvin)
% kb = 1.380658e-23;                 % konstanta boltzmann
%   N0 = (2*q*fm*be)+((4*kb*T*be)/rl); % Noise Density
N0 = 2*q*Ib; % white noise direceiver (hal 99)

Rb = 1e9;         % Bit rate (Hz) kecepatan
Tb = 1/Rb;        % bit duration (makin tinggi Rb, maka makin cepet pengiriman, 
nsamp = 10;       % samples per symbols
Tsamp = Tb/nsamp; % sampling time

EbN0 = [0:1:19];          % signal-to-noise ratio in dB.
SNR = 10.^(EbN0./10); % signal-to-noise ratio
konv=qfunc(sqrt(10.^(EbN0/10))); %konversi dr SNR ke numerik

%% Manggil file

 namafile='PrxFix';
 sheetQ='BER';
 data=xlsread(namafile,sheetQ);
 prx =data(1,1)/1000;  
 
%% OFDM parameter
M = 16; % input('masukan nilai M-array = '); 
NFFT = 128;
subcarrier = (NFFT/2)-1; % input('Jumlah subcarriers = ');
NCP = 16; % input('CP length NCP = ');
k = log2(M);
DC_bias = 16; % DC bias
bit = 1984248/2; % 61.440 for subcarrier = 64,128,256 ; Mapper = 16QAM,64QAM,QPSK ; NFFT = 128,256,512

iterasi = 1;

% %%
% lx=5; ly=5; lz=3; % room dimension in meter
% 
% h=2.15;           %the distance between source and receiver plane
% 
% % for one LED simulation located at the central of the room,
% XT = 0;
% YT = 0;
% XR = 0;
% YR = 0;
% 
% D1=sqrt((XR-XT).^2+(YR-YT).^2+h^2); % distance vector from source 1
% cosphi_A1=h./D1;                              % angle vector
% 
% H_A1=(ml+1)*Adet.*cosphi_A1.^(ml+1)./(2*pi.*D1.^2); % channel DC gain for source 1
% % H_A1=1e-9;    
    %% Transmitter %%
% iterasi monte carlo   
for mm = 1:iterasi
    %% binary data
    t_data = randi([0 1],1,bit);

    %% bit to decimal
    for i=1:length(t_data)/k
        t_bi2de(1,i)=bi2de(t_data(1,k*i+1-k:k*i),'left-msb');
    end
    
    %% modulasi QAM
    t_mod=qammod(t_bi2de,M);
 
 %% serial to pararel
    t_s2p= [];
    for j = 1:length(t_mod)/(subcarrier)
        t_s2p(j,1:subcarrier) = t_mod(((subcarrier)*j+1)-(subcarrier):(subcarrier)*j);
    end;
    
    % Hermitian symmetry
        HS_signal = zeros ((size(t_s2p,1)), NFFT);
        HS_signal(:,2:(NFFT/2)) = t_s2p(:,1:subcarrier);
        HS_signal(:,(NFFT/2)+2:end)=fliplr(conj(t_s2p(:,1:(subcarrier))));
    
    % IFFT 
        t_ifft_signal = (sqrt(NFFT))*ifft(fftshift(HS_signal.')).';
        %t_ifft_signal = ifft(HS_signal);
        
    %% Add CP
    t_cp = [t_ifft_signal(:,(NFFT-NCP)+1:NFFT) t_ifft_signal];
    
    %% pararel to serial
    t_p2s=[];
    for j=1:size(t_cp,1);
        t_p2s(1,((NFFT+NCP)*j+1)-(NFFT+NCP):(NFFT+NCP)*j)=t_cp(j,:);
    end;
    
    %% Add DC bias and Clipping
    % Add DC bias
    t_dc = t_p2s+DC_bias;

    %clipping
    t_clip = 0.5*(t_dc + abs(t_dc));
    
    %% amplifier
    %BO=40; %input backoff
    %p=2; %mengatur kehalusan transisi dari daerah linier ke daerah batas saturasi
    %rata=mean(abs(t_dc).^2);
    %a=(10^(BO/10)*rata)^0.5; %saturasi
    %t_amp=t_dc./((1+abs(t_dc)/a).^(2*p)).^(1/(2*p)); %AM/PM dianggap kecil (diabaikan)
    %t_amp = t_dc./((1+((abs(t_dc)/5)).^2).^0.5);
    
    %% kanal
    for ebno = 1:length(EbN0)
    %Snr = EbN0(ebno) 
    
        P_avg(ebno) = prx(ebno);          % average transmitted optical power
        i_peak(ebno) = 2*R*P_avg(ebno);              % Peak Electrical amplitude
        Ep(ebno) = i_peak(ebno).^2*Tb;                % Peak energy (Energy per bit is Ep/2)
        
        sgma(ebno) = sqrt(N0/2/Tsamp);            % noise variance
        %sgma(ebno) = sqrt(N0*Ep(ebno)/2);
%         sgma(l) = i_peak(l)/sqrt(2)*sqrt(nsamp/(2*SNR(l)));
    
        %pt = ones(1,nsamp)*i_peak(ebno); % tranmitter filter
        %rt = pt;                      % receiver filter matched to pt
    
        Tx_signal = t_clip.*i_peak(ebno);
        Rx_signal = R*Tx_signal+sgma(ebno)*randn(1,length(Tx_signal)); 
        Rx_pd = Rx_signal./i_peak(ebno);
        
            
%         sgma(ebno) = sqrt(N0/2/Tsamp);            % noise variance
%          sgma(ebno) = sqrt(N0*Ep(ebno)/2); %choice 1 (dayanya kecilin, sampe 1e-4 paling kecil 1watt
%          sgma(ebno) = i_peak(ebno)/sqrt(2)*sqrt(nsamp/(2*SNR(ebno)));
%         choice 2 (dayanya gedein, sampe 1e-4 paling max 25watt
            
        %% RECEIVER
        % Removing DC bias
        r_dc = Rx_pd - DC_bias;
    
        %% serial to pararel
        r_s2p = [];
        for j = 1:length(r_dc)/(NFFT+NCP)
            r_s2p(j,1:NFFT+NCP) = r_dc(((NFFT+NCP)*j+1)-(NFFT+NCP):(NFFT+NCP)*j);
        end;
        %r_s2p = reshape(outawgn,(NFFT+NCP),length(outawgn)/(NFFT+NCP));
          
        %% Removing CP
        %x_p_cpr = (x_p(:,[NCP+1:NFFT]));
        r_rem_cp = r_s2p(:,NCP+1:end);
    
       %% FFT
        r_fft = (sqrt(NFFT)/NFFT)*fftshift(fft(r_rem_cp.')).';
        r_fft_signal = [];
        r_fft_signal = [r_fft_signal r_fft(:,2:NFFT/2)];
    
        %% pararel to serial
        %[c d] = size (r_fft_signal); 
        %r_p2s = reshape(r_fft_signal,1,c*d);
        r_p2s=[];
        for j=1:size(r_fft_signal,1)
            r_p2s(1,((subcarrier)*j+1)-(subcarrier):(subcarrier)*j)=r_fft_signal(j,:);
        end;

        %% Demodulator
        r_demod = qamdemod(r_p2s,M);

        %% Decimal to Bit
        r_bit=de2bi(r_demod,'left-msb');
        r_bit_ser=[];
        for i=1:size(r_bit,1)
            r_bit_ser(1,k*i+1-k:k*i)=r_bit(i,:);
        end;
    
        % Calculating BER
        [nErr bErr(ebno,mm)] = symerr(t_data,r_bit_ser); 
   
    end; 
end;
semilogy(EbN0,bErr','b-*');
hold on;
% grid on;
% semilogy(EbN0,qfunc(sqrt(10.^(EbN0/10))),'r-X'); % theoretical ber, 'mx-');
% legend('DCO-OFDM','theory');
title('BER');
xlabel('Eb/No');
ylabel('Bit Error Rate');