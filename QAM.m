clear all;close all;clc;
%%  inisiasi parameter
% P_LED=7000;                       % Total transmitted power (mW)
Adet=1e-4;                          % detector physical area of a PD (M)

q = 1.6e-19;                       % Charge of Electron
R = 0.55;                          % Photodetector responsivity
be = 2e9;                          % Bandwidth 2mhz normalnya
rl = 30;                           % impedence/resistansi
fm = 0.7;                          % noise figure
T = 300;                           % temperature (kelvin) ngaruh ke perangkat receiver (makin panas makin bahaya)
kb = 1.380658e-23;                 % konstanta boltzmann
Ib=202e-6;                         % background noise
N0 = 2*q*Ib;
% N0 = (2*q*fm*be)+((4*kb*T*be)/rl); % Noise Density

Rb = 1e9;         % Bit rate (Hz)
Tb = 1/Rb;        % bit duration
nsamp = 10;       % samples per symbols
Tsamp = Tb/nsamp; % sampling time

EbN0 = [0:1:14];          % signal-to-noise ratio in dB.
SNR = 10.^(EbN0./10); % signal-to-noise ratio
% konv=qfunc(sqrt(10.^(EbN0/10)));

%% Manggil file

 namafile='PrxFix';
 sheetQ='BER';
 data=xlsread(namafile,sheetQ);
 prx =data(1,1)/1000;  
 
%% OFDM parameter
M = 4; % input('masukan nilai M-array = '); 
NFFT = 128;
subcarrier = (NFFT/2)-1; % input('Jumlah subcarriers = ');
NCP = 16; % input('CP length NCP = ');
k = log2(M); % 2 pangkat
DC_bias = 16 ; % DC bias
bit = 1984248/2; %1960560 (16psk); %1984248/2 (rana); %1984122 (8psk);% 61.440 for subcarrier = 64,128,256 ; Mapper = 16QAM,64QAM,QPSK ; NFFT = 128,256,512

literasi = 1 ;
    %% Transmitter %%
    
for mm = 1:literasi
    %% binary data
    t_data = randi([0 1],1,bit);

    %% bit to decimal
    for i=1:length(t_data)/k
        t_bi2de(1,i)=bi2de(t_data(1,k*i+1-k:k*i),'left-msb');
    end
    
    %% modulasi QAM
    t_mod=pskmod(t_bi2de,M);
 
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
    
        %% kanal
   
    for l = 1:length(prx)
    %Snr = EbN0(ebno) 
        %p_vlc = t_dc*P_LED;
        
        P_avg(l) = prx(l);          % average transmitted optical power
        i_peak(l) = 2*R*P_avg(l);              % Peak Electrical amplitude
        Ep(l) = i_peak(l).^2*Tb;                % Peak energy (Energy per bit is Ep/2)
        
        sgma(l) = sqrt(N0/2/Tsamp);            % noise variance
        %sgma(ebno) = sqrt(N0*Ep(ebno)/2);
%         sgma(l) = i_peak(l)/sqrt(2)*sqrt(nsamp/(2*SNR(l)));
    
        %pt = ones(1,nsamp)*i_peak(ebno); % tranmitter filter
        %rt = pt;                      % receiver filter matched to pt
    
        Tx_signal = t_clip.*i_peak(l);
        Rx_signal = R*Tx_signal+sgma(l)*randn(1,length(Tx_signal)); 
        Rx_pd = Rx_signal./i_peak(l);
        
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
        for j=1:size(r_fft_signal,1);
            r_p2s(1,((subcarrier)*j+1)-(subcarrier):(subcarrier)*j)=r_fft_signal(j,:);
        end;

        %% Demodulator
        r_demod = pskdemod(r_p2s,M);

        %% Decimal to Bit
        r_bit=de2bi(r_demod,'left-msb');
        r_bit_ser=[];
        for i=1:size(r_bit,1)
            r_bit_ser(1,k*i+1-k:k*i)=r_bit(i,:);
        end;
    
        % Calculating BER
        [nErr bErr(l,mm)] = symerr(t_data,r_bit_ser); 
   end;
end;

% mean_ber_ori = mean(bErr,2);
% xlswrite('data_ber_QAM_ori', mean_ber_ori);
% 
%  semilogy(data(:,6),data(:,2),'b-*','linewidth',1),grid on; hold on;
%  title('BER');
%  xlabel('Eb/No');
%  ylabel('Bit Error Rate');