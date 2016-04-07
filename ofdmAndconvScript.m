clc, clear, close all;
%% 802.11a Spec - Stephen Leone, Noah Santacruz
totPak = 2^4;
packetSize = 2304;
EbNo_Vec = 0:20;
traceback=12;
sigpow=-10;

berVec = zeros(totPak,size(EbNo_Vec,2));
trel = poly2trellis(7,[133,171]); %[1011011,1111001] in binary
spect = distspec(trel,2);
H1 = [1 .3];
H2 = [.2 .5 1 .3];
H3 = [1 .05 .3 .4 .2 .05];
R=1/2;
for H = 1:3
    figure;
    switch H
        
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
    for datarate=[6,12,24];
        switch datarate
            case 6,
                M=2;
                pointStyleZF='*b-';
                pointStyleMMSE='^b-';
            case 12,
                M=4;
                R=1/2;
                pointStyleZF='*r-';
                pointStyleMMSE='^r-';
            case 24,
                M=16;
                R=1/2;
                pointStyleZF='*g-';
                pointStyleMMSE='^g-';
            otherwise,
                disp(['Rate ' num2str(datarate) ' not supported']);
                exit(1);
        end

        disp(['Simulating Rate ' num2str(datarate) ' Mbps']);
        k = log2(M);
        snrs = EbNo_Vec + 10.*log10(k*R*48/64);

        for pak=1:totPak/k

            for ebno = 1:length(EbNo_Vec)
                bits = randi([0 1],k*packetSize,1);
                txVitBits = convenc(bits,trel);
                msg = bi2de(reshape(txVitBits,k,length(txVitBits)/k).','left-msb')';
                txMod = qammod(msg,M,0,'gray');
                frameCount=numel(txMod)/48;
                txtones = reshape(txMod,48,frameCount);
                txOFDM = [txtones(1:5,:);zeros(1,frameCount);
                          txtones(6:18,:);zeros(1,frameCount);
                          txtones(19:24,:);zeros(1,frameCount);
                          txtones(25:30,:);zeros(1,frameCount);
                          txtones(31:43,:);zeros(1,frameCount);
                          txtones(44:48,:)];
                txDTFTin = [txOFDM(27:53,:);
                          zeros(11,frameCount);
                          txOFDM(1:26,:)];
                FFTsym = ifft(txDTFTin,64,1);
                GI =FFTsym(49:64,:);
                OFDMsym = [GI;FFTsym];
                OFDMsig=reshape(OFDMsym,80*frameCount,1);

                %channel

                OFDMsig = conv(H,OFDMsig);


                noisyTx = awgn(OFDMsig(1:80*frameCount),snrs(ebno),sigpow);

                %receive

                rxOFDMsym = reshape(noisyTx,80,frameCount);
                rxFFTin = rxOFDMsym(17:80,:);

                rxFFTZF = fft(rxFFTin,64,1)./repmat(fft(H',64,1),1,frameCount);
                N0 = 10^(snrs(ebno)/10);
                rxFFTMMSE = fft(rxFFTin,64,1)./(repmat(fft(H',64,1),1,frameCount)+N0);

                for eq = 0:1
                    if eq
                        rxFFTout = rxFFTZF;
                    else
                        rxFFTout = rxFFTMMSE;
                    end
                    rxOFDM = [rxFFTout(39:64,:);rxFFTout(1:27,:)];
                    rxtones = [rxOFDM(1:5,:);
                               rxOFDM(7:19,:);
                               rxOFDM(21:26,:);
                               rxOFDM(28:33,:);
                               rxOFDM(35:47,:);
                               rxOFDM(49:53,:)];

                    rxDemod = reshape(rxtones,frameCount*48,1);
                    rx = qamdemod(rxDemod,M,0,'gray');
                    rx = de2bi(rx,'left-msb');
                    rxBits = reshape(rx.',numel(rx),1);
                    rxVitBits = vitdec(rxBits,trel,traceback,'cont','hard'); % Decode

                    if eq
                        [~,berVecZF(pak,ebno)] = biterr(bits(1:end-traceback), rxVitBits(traceback+1:end));
                    else
                        [~,berVecMMSE(pak,ebno)] = biterr(bits(1:end-traceback), rxVitBits(traceback+1:end));
                    end
                end
            end
        end

        berZF = mean(berVecZF(1:totPak/k,:),1);
        berMMSE = mean(berVecMMSE(1:totPak/k,:),1);
        semilogy(EbNo_Vec, berZF,pointStyleZF)
        hold on
        semilogy(EbNo_Vec, berMMSE,pointStyleMMSE)

    end
    %6zf, 6mmse, 12zf,12mmse, 24zf, 24mmse
    legend('6zf', '6mmse', '12zf','12mmse', '24zf', '24mmse');
    xlabel('EbNo');
    ylabel('Bit Error Rate');
    title(titl);
    axis([0 20 1e-5 1e0])
end