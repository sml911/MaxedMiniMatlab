clc, clear, close all;
%% 802.11a Spec - Stephen Leone, Noah Santacruz
totPak = 2560;
packetSize = 2304;
SNR = 20

H1 = [1 0.8]
H2 = [1 .5 1 .8]
H3 = [1 .05 .3 .4 .2 .05]

H_Bertargets = [1E-6,1E-5,1E-5];

zfMs = [64,8,16];
mmseMs = [64,8,16];

bitsVecZF = zeros(3,totPak);
bitsVecMMSE = zeros(3,totPak);
berVecZF = zeros(3,totPak);
berVecMMSE = zeros(3,totPak);

h = waitbar(0,'...');
for Htype = 1:3
    switch Htype
        
        case 1, 
            H=H1;
            titl='Channel1';
        case 2, 
            H=H2;
            titl='Channel2';
        case 3, 
            H=H3;
            titl='Channel3';
    end
    H = H/sqrt(var(H));
    
    
    for Mtype=1:2  % 2== number of equalization schemes
        if Mtype == 1
            M = zfMs(Htype);
        elseif Mtype == 2
            M = mmseMs(Htype);
        end
    
        k = log2(M);

        for pak=1:totPak
            bits = randi([0 1],k*packetSize,1);
            txVitBits = bits; %convenc(bits,trel);
            msg = bi2de(reshape(txVitBits,k,length(txVitBits)/k).','left-msb')';
            txMod = qammod(msg,M,0,'gray');
            
            N0linear = std(txMod(:))/10^(SNR/10);
            N0db = 10*log10(N0linear);
            frameCount=numel(txMod)/48;
            OFDMsig = OFDMmod(txMod,frameCount);
            %channel

            OFDMsig = conv(H,OFDMsig);

            noise = wgn(80*frameCount,1,N0db,'complex');
            noisyTx = OFDMsig(1:80*frameCount) + noise;

            %receive

            rxOFDMsym = reshape(noisyTx,80,frameCount);
            rxFFTin = rxOFDMsym(17:80,:);


            rxFFTZF = fft(rxFFTin,64,1)./repmat(fft(H',64,1),1,frameCount);


            rxFFTMMSE = fft(rxFFTin,64,1)./(repmat(fft(H',64,1),1,frameCount)+N0linear);
            if Mtype == 1
                rxFFTout = rxFFTZF;
            elseif Mtype == 2
                rxFFTout = rxFFTMMSE;
            end
            rxDemod = OFDMdemod(rxFFTout,frameCount);
            rxDemod = rxDemod * std(txMod) / std(rxDemod);
            rx = qamdemod(rxDemod,M,0,'gray');
            rx = de2bi(rx,'left-msb');
            rxBits = reshape(rx.',numel(rx),1);
            rxVitBits = rxBits; 

            if Mtype == 1
                [bitsLost,berVecZF(Htype,pak)] = biterr(bits, rxVitBits);
                if bitsLost == 0
                    bitsVecZF(Htype,pak) = packetSize * k - bitsLost;
                end
            elseif Mtype == 2
                [bitsLost,berVecMMSE(Htype,pak)] = biterr(bits, rxVitBits);
                if bitsLost == 0
                   bitsVecMMSE(Htype,pak) = packetSize * k - bitsLost;
                end
            end
        end
    end
    

    waitbar(Htype/3);
end
close(h);

berZF = mean(berVecZF,2)'
berMMSE = mean(berVecMMSE,2)'

bitsZF = sum(bitsVecZF,2)'
bitsMMSE = sum(bitsVecMMSE,2)'