clc;clear;close all;

M = 2;
k = log2(M);
N = 10E5;
SNR = 20;
N0=-SNR;%Transmit power is normalized to 0dBW
numTx = 2;
T = 1e-5;
dopshift = 0;

numIter = 100;
berPlain = zeros(1,numIter);
berPrecode = zeros(1,numIter);
berZF = zeros(1,numIter);
berMMSE = zeros(1,numIter);

for jj = 1:numIter
    bits = randi([0,1],numTx,N);
    msg = zeros(2,size(bits,2)/k);
    for ii = 1:numTx
        msg(ii,:) = bi2de(reshape(bits(ii,:),k,size(bits,2)/k).','left-msb')';
    end
    xTilde = qammod(msg,M);
    xTilde = xTilde./sqrt(var(xTilde(:))); %normalize power

    chan1 = rayleighchan(T,dopshift);
    chan2 = rayleighchan(T,dopshift);
    chan3 = rayleighchan(T,dopshift);
    chan4 = rayleighchan(T,dopshift);

    H = [chan1.PathGains,chan2.PathGains;chan3.PathGains,chan4.PathGains];
    [U,S,V] = svd(H);


    noise = wgn(size(xTilde,1),size(xTilde,2),N0,'complex');

    %PLAIN MIMO
    yPlain = H*xTilde;
    yPlainN = yPlain + noise;
    msgRxPlain = qamdemod(yPlainN,M);
    [~,berPlain(jj)] = biterr(msg,msgRxPlain);

    %PRECODING
    xPrecode = V'*xTilde;
    yPrecode = H*xPrecode;
    yPrecodeN = yPrecode + noise;
    scatterplot(yPrecodeN(:));
    title('Precode - Before');
    yTilde = U' * yPrecodeN;
    scatterplot(yTilde(:));
    title('Precode - After');
    msgRxPrecode = qamdemod(yTilde,M);
    [~,berPrecode(jj)] = biterr(msg,msgRxPrecode);
    %ZEROFORCING
    W = pinv(H);
    scatterplot(yPlainN(:));
    title('ZF - Before');
    yZFN = W*yPlainN;
    scatterplot(yZFN(:));
    title('ZF - After');
    msgRxZF = qamdemod(yZFN,M);
    [~,berZF(jj)] = biterr(msg,msgRxZF);

    %MMSE
    W = inv(H'*H + N0)*H';
    scatterplot(yPlainN(:));
    title('MMSE - Before');
    yMMSEN = W*yPlainN;
    scatterplot(yMMSEN(:));
    title('MMSE - After');
    msgRxMMSE = qamdemod(yMMSEN,M);
    [~,berMMSE(jj)] = biterr(msg,msgRxMMSE);
end
