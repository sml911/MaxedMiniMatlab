clc;clear;close all;

M = 2;
k = log2(M);
N = 1E6; %number of bits per iteration
SNR = 1:20;
N0=-SNR;%Transmit power is normalized to 0dBW
numTx = 2;
T = 1e-5;
dopshift = 0;

numIter = 1;
berPlain   = zeros(length(SNR),numIter);
berPrecode = zeros(length(SNR),numIter);
berZF      = zeros(length(SNR),numIter);
berMMSE    = zeros(length(SNR),numIter);

h = waitbar(0,'...');
for kk = 1:length(SNR)
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


        noise = wgn(size(xTilde,1),size(xTilde,2),N0(kk),'complex');

        %PLAIN MIMO
        yPlain = H*xTilde;
        yPlainN = yPlain + noise;
        msgRxPlain = qamdemod(yPlainN,M);
        [~,berPlain(kk,jj)] = biterr(msg,msgRxPlain);

        %PRECODING
        xPrecode = V*xTilde;
        yPrecode = H*xPrecode;
        yPrecodeN = yPrecode + noise;
            
        yTilde = (U' * yPrecodeN);
        yTilde = S\yTilde;
        if numIter == 1
            scatterplotColorful(yPrecodeN,xTilde,yTilde,'Precoding');
        end
        msgRxPrecode = qamdemod(yTilde,M);
        [~,berPrecode(kk,jj)] = biterr(msg,msgRxPrecode);
        %ZEROFORCING
        W = (H'*H)\H';

        yZFN = W*yPlainN;
        if numIter == 1
            scatterplotColorful(yPlainN,xTilde,yZFN,'Zero Forcing');
        end
        msgRxZF = qamdemod(yZFN,M);
        [~,berZF(kk,jj)] = biterr(msg,msgRxZF);

        %MMSE
        W = (H'*H + 10^(N0(kk)/10))\H';
        yMMSEN = W*yPlainN;

        if (numIter == 1) 
            scatterplotColorful(yPlainN,xTilde,yMMSEN,'MMSE');
        end
        msgRxMMSE = qamdemod(yMMSEN,M);
        [~,berMMSE(kk,jj)] = biterr(msg,msgRxMMSE);

    end
    waitbar(kk/length(SNR));
end
close(h);

figure;
semilogy(SNR,mean(berPlain,2)','r'); hold on;
semilogy(SNR,mean(berPrecode,2)','b');
semilogy(SNR,mean(berZF,2)','g');
semilogy(SNR,mean(berMMSE,2)','m');
hold off;
legend('plain','precoding','zf','mmse');

mean(berPrecode,2)
