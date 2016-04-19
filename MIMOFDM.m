clc, clear, close all;

M = 2;
k = log2(M);
SNR = 20;
totPak = 64;
N = 2304*2;
EbNo_Vec = 0:2:10;
traceback=12;
sigpow=1;
numTx = 2;

berVec = zeros(totPak,size(EbNo_Vec,2));
trel = poly2trellis(7,[133,171]); %[1011011,1111001] in binary
spect = distspec(trel,2);
R=1/2;

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
%txVitBits = convenc(bits,trel);
%msg = bi2de(reshape(txVitBits,k,length(txVitBits)/k).','left-msb')';
msg = zeros(2,size(bits,2)/k);
for ii = 1:numTx
    msg(ii,:) = bi2de(reshape(bits(ii,:),k,size(bits,2)/k).','left-msb')';
end
xTilde = qammod(msg,M);

%OFDM MOD
frameCount=numel(xTilde)/48/numTx;

OFDMsig1 = OFDMmod(xTilde(1,:),frameCount);
OFDMsig2 = OFDMmod(xTilde(2,:),frameCount);


%MIMO Channel
N0linear = std([OFDMsig1;OFDMsig2])/10^(SNR/10);
N0db = 10*log10(N0linear);
noise = wgn(numTx,length(OFDMsig1),N0db,'complex');
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
rxZF1 = squeeze(yZFN(1,:,:)).';    
rxZF2 = squeeze(yZFN(2,:,:)).';  
rxMMSE1 = squeeze(yMMSEN(1,:,:)).';    
rxMMSE2 = squeeze(yMMSEN(2,:,:)).'; 

OFDMdemodZF(1,:) = OFDMdemod(rxZF1,frameCount).';
OFDMdemodZF(2,:) = OFDMdemod(rxZF2,frameCount).';
OFDMdemodMMSE(1,:) = OFDMdemod(rxMMSE1,frameCount).';
OFDMdemodMMSE(2,:) = OFDMdemod(rxMMSE2,frameCount).';

for eq = 0:1
    if eq
        Yout = OFDMdemodZF;
    else
        Yout = OFDMdemodMMSE;
    end
    Yout = Yout * std(xTilde(1,:)) / std(Yout(1,:));
    rx = qamdemod(Yout,M,0,'gray');
    rx = de2bi(rx,'left-msb');
    rxBits = reshape(rx.',numel(rx),1);
    bits = reshape(bits,numel(rx),1);
    %rxVitBits = vitdec(rxBits,trel,traceback,'cont','hard'); % Decode

    if eq
        [~,berVecZF] = biterr(bits, rxBits,'overall')
    else
        [~,berVecMMSE] = biterr(bits, rxBits,'overall')
    end
end 






