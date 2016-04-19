clc, clear, close all;

M = 2;
k = log2(M);
SNR = 20;
totPak = 64;
N = 2304;
EbNo_Vec = 0:2:10;
traceback=12;
sigpow=1;
numTx = 2;

berVec = zeros(totPak,size(EbNo_Vec,2));
trel = poly2trellis(7,[133,171]); %[1011011,1111001] in binary
spect = distspec(trel,2);

h11 = [1 .3 0 0 0 0];
h12 = [.2 .5 1 .3 0 0];
h21 = [1 .05 .3 .4 .2 .05];
h22 = [.3 1 0 0 0 0];

h11fft = fft(h11',64,1);
h12fft = fft(h12',64,1);
h21fft = fft(h21',64,1);
h22fft = fft(h22',64,1);

BigHFFT = zeros(numTx,numTx,64);
BigHFFT(1,1,:) = h11fft;
BigHFFT(1,2,:) = h12fft;
BigHFFT(2,1,:) = h21fft;
BigHFFT(2,2,:) = h22fft;



bits = randi([0,1],numTx,N);
msg = zeros(2,size(bits,2)/k);
for ii = 1:numTx
    msg(ii,:) = bi2de(reshape(bits(ii,:),k,size(bits,2)/k).','left-msb')';
end
xTilde = qammod(msg,M);

N0linear = std(xTilde(:))/10^(SNR/10);
N0db = 10*log10(N0linear);
noise = wgn(size(xTilde,1),64/48*size(xTilde,2),N0db,'complex');

%OFDM MOD
frameCount=numel(xTilde)/48/numTx;

OFDMsig1 = OFDMmod(xTilde(1,:),frameCount);
OFDMsig2 = OFDMmod(xTilde(2,:),frameCount);

Y1 = conv(h11,OFDMsig1) + conv(h21,OFDMsig2) ;
Y2 = conv(h12,OFDMsig1) + conv(h22,OFDMsig2) ;
Y1 = Y1(1:length(OFDMsig1));
Y2 = Y2(1:length(OFDMsig2));
Y1N = Y1.' + noise(1,:);
Y2N = Y2.' + noise(2,:);

%demodulate using OFDM
rxOFDMsym1 = reshape(Y1N,80,frameCount);
rxFFTin1 = rxOFDMsym1(17:80,:);
rxFFT1 = fft(rxFFTin1,64,1);
rxOFDMsym2 = reshape(Y2N,80,frameCount);
rxFFTin2 = rxOFDMsym2(17:80,:);
rxFFT2 = fft(rxFFTin2,64,1);


yZFN = zeros(numTx,frameCount,64);
yMMSEN = zeros(numTx,frameCount,64);
for k = 1:64
    Hk = BigHFFT(:,:,k);
    WZFk = (Hk'*Hk)\Hk;
    WMMSEk = (Hk'*Hk + N0linear)\Hk;
    yZFN(:,:,k) = WZFk*[rxFFT1(k,:);...
                        rxFFT2(k,:)];
    yMMSEN(:,:,k) = WMMSEk * [rxFFT1(k,:);...
                             rxFFT2(k,:)];
end
    
    

OFDMdemod1 = OFDMdemod(rxFFT1,frameCount);
OFDMdemod2 = OFDMdemod(rxFFT2,frameCount);


%rxFFTZF1 = rxFFT1./repmat(fft(H',64,1),1,frameCount);
%rxFFTMMSE1 = rxFFT1./(repmat(fft(H',64,1),1,frameCount)+N0linear);
rxFFT2 = fft(rxFFTin2,64,1);
%rxFFTZF2 = rxFFT2./repmat(fft(H',64,1),1,frameCount);
%rxFFTMMSE2 = rxFFT2./(repmat(fft(H',64,1),1,frameCount)+N0linear);






