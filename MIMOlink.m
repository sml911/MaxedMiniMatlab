clc;clear;close all;
%% MIMO - Part 1 - Stephen Leone, Noah Santacruz
N = 24000; %number of bits per iteration

%N0=-SNR;%Transmit power is normalized to 0dBW
numTx = 2;
dopshift = 0;
bitsThresh = 0;

numIter = 1000;


berPlain   = zeros(3,numIter);
berPrecode = zeros(3,numIter);
berZF      = zeros(3,numIter);
berMMSE    = zeros(3,numIter);

bitsPlain   = zeros(3,numIter);
bitsPrecode = zeros(3,numIter);
bitsZF      = zeros(3,numIter);
bitsMMSE    = zeros(3,numIter);

H1 = [1   0;...
      0.01   0.8];
H2 = [1.2 i;...
      0.7 1.5];

H3 = [0.85 -0.4;...
      0.4   0.95];...

Hs(:,:,1) = H1;
Hs(:,:,2) = H2;
Hs(:,:,3) = H3;
Hs

disp('Channel BER Targets')
H_BERtargets = [1E-5 1E-5 1E-5]

SNR = [20 15 12]
Ms = [64,16,4] %modulation for each channel
precodeMs = [64,16,4];
zfMs = [64,16,4];
mmseMs = [64,16,4];


h = waitbar(0,'...');
for Htype = 1:3
    H = Hs(:,:,Htype);
    for Mtype = 1:3
        if Mtype == 1
            M = precodeMs(Htype);
        elseif Mtype == 2
            M = zfMs(Htype);
        elseif Mtype == 3
            M = mmseMs(Htype);
        end
        k = log2(M);
        
        for jj = 1:numIter
            bits = randi([0,1],numTx,N);
            msg = zeros(2,size(bits,2)/k);
            
            for ii = 1:numTx
                msg(ii,:) = bi2de(reshape(bits(ii,:),k,size(bits,2)/k).','left-msb')';
            end
            xTilde = qammod(msg,M);
            %xTilde = xTilde./std(xTilde(:)); %normalize power

            %chan1 = rayleighchan(T,dopshift);
            %chan2 = rayleighchan(T,dopshift);
            %chan3 = rayleighchan(T,dopshift);
            %chan4 = rayleighchan(T,dopshift);

            %H = [chan1.PathGains,chan2.PathGains;chan3.PathGains,chan4.PathGains];
            [U,S,V] = svd(H);

            N0linear = std(xTilde(:))/10^(SNR(Htype)/10);
            N0db = 10*log10(N0linear);
            noise = wgn(size(xTilde,1),size(xTilde,2),N0db,'complex');

            %PLAIN MIMO
            yPlain = H*xTilde;
            yPlainN = yPlain + noise;
%             msgRxPlain = qamdemod(yPlainN,M);
%             [bitsPlain(Htype,jj),berPlain(Htype,jj)] = biterr(msg,msgRxPlain);

            
            if Mtype == 1
                %PRECODING
                xPrecode = V*xTilde;
                yPrecode = H*xPrecode;
                yPrecodeN = yPrecode + noise;

                yTilde = (U' * yPrecodeN);
                yTilde = S\yTilde;
                if numIter == 1
                    scatterplotColorful(yPrecodeN,xTilde,yTilde,['Precoding - ' int2str(Htype)]);
                end
                msgRxPrecode = qamdemod(yTilde,M);
                %Why Not This?
                %bitsRxPrecode = de2bi(msgRxPrecode,'left-msb');
                %bitsRxPrecode = reshape(bitsRxPrecode.',numel(bitsRxPrecode)/2,2);
                %[bitsLost,berPrecode(Htype,jj)] = biterr(bits',bitsRxPrecode);
                [bitsLost,berPrecode(Htype,jj)] = biterr(msg,msgRxPrecode);
                if bitsLost <= bitsThresh
                    bitsPrecode(Htype,jj) = N*k - bitsLost;
                end
            elseif Mtype == 2
                %ZEROFORCING
                W = (H'*H)\H';

                yZFN = W*yPlainN;
                if numIter == 1
                    scatterplotColorful(yPlainN,xTilde,yZFN,['Zero Forcing - ' int2str(Htype)]);
                end
                msgRxZF = qamdemod(yZFN,M);
                [bitsLost,berZF(Htype,jj)] = biterr(msg,msgRxZF);
                if bitsLost <= bitsThresh
                    bitsZF(Htype,jj) = N*k - bitsLost;
                end
            elseif Mtype == 3
                %MMSE
                W = (H'*H - N0linear*eye(2))\H';
                yMMSEN = W*yPlainN;
                yMMSEN = yMMSEN * std(xTilde(:)) / std(yMMSEN(:));
                if (numIter == 1) 
                    scatterplotColorful(yPlainN,xTilde,yMMSEN,['MMSE - ' int2str(Htype)]);
                end
                msgRxMMSE = qamdemod(yMMSEN,M);
                [bitsLost,berMMSE(Htype,jj)] = biterr(msg,msgRxMMSE);
                if bitsLost <= bitsThresh
                    bitsMMSE(Htype,jj) = N*k - bitsLost;
                end
            end
        end
    end
    waitbar(Htype/length(1:3));
end
close(h);
%meanBerPlain = mean(berPlain,2)'
meanBerPrecode = mean(berPrecode,2)'
meanBerZF = mean(berZF,2)'
meanBerMMSE = mean(berMMSE,2)'

%meanBitsPrecode = sum(bitsPrecode,2)'
%meanBitsZF = sum(bitsZF,2)'
%meanBitsMMSE = sum(bitsMMSE,2)'

% figure; hold on;
% plot(1:3,meanBerPlain,'r');
% plot(1:3,meanBerPrecode,'b');
% plot(1:3,meanBerZF,'g');
% plot(1:3,meanBerMMSE,'m');
% legend('plain','precoding','zf','mmse');
% hold off;
