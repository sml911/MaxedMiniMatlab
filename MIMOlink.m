clc;clear;close all;

M = 2;
k = log2(M);
N = 10E5;
SNR = 20;
N0=-SNR;%Transmit power is normalized to 0dBW
numTx = 2;
T = 1e-5;
dopshift = 0;

bits = randi([0,1],numTx,N);
msg = zeros(2,size(bits,2)/k);
for ii = 1:numTx
    msg(ii,:) = bi2de(reshape(bits(ii,:),k,size(bits,2)/k).','left-msb')';
end
xTilde = qammod(msg,M);

chan1 = rayleighchan(T,dopshift);
chan2 = rayleighchan(T,dopshift);
chan3 = rayleighchan(T,dopshift);
chan4 = rayleighchan(T,dopshift);

H = [chan1.PathGains,chan2.PathGains;chan3.PathGains,chan4.PathGains];
[U,S,V] = svd(H);


noise = wgn(size(xTilde,1),size(xTilde,2),N0);

%PLAIN MIMO
yPlain = H*xTilde;
yPlainN = yPlain + noise;
%PRECODING
xPrecode = V'*xTilde;
yPrecode = H*xPrecode;
yPrecodeN = yPrecode + noise;
yTilde = U' * yPrecodeN;
%ZEROFORCING
W = pinv(H);
yZF = W*yPlainN;
%MMSE
W = inv(H'*H + N0)*H'
