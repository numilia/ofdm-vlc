clear all;close all;clc;

M = 4; % Alphabet size
EbN0_min=0;EbN0_max=10;step=2;
SNR=[];SER=[];
for EbN0 = EbN0_min:step:EbN0_max
SNR_dB=EbN0 + 3; %for QPSK Eb/N0=0.5*Es/N0=0.5*SNR
x = randi(1000000,1,M);
y=modulate(modem.qammod(M),x);
ynoisy = awgn(y,SNR_dB,'measured');
z=demodulate(modem.qamdemod(M),ynoisy);
[num,rt]= symerr(x,z);
SNR=[SNR EbN0];
SER=[SER rt];
end;
semilogy(SNR,SER);grid;titel('Symbol error rate for QPSK over AWGN');
xlabel('E_b/N_0');ylabel('SER');