clc, clear, close all;
%% 802.11a Spec - Stephen Leone, Noah Santacruz
totPak = 2^8;
packetSize = 2304;
EbNo_Vec = 0:20;
traceback=12;

berVec = zeros(totPak,size(EbNo_Vec,2));
trel = poly2trellis(7,[133,171]); %[1011011,1111001] in binary
spect = distspec(trel,2);

for datarate=[6, 9, 12, 18, 24, 36, 48, 54]
    switch datarate
        case 6,
            M=2;
            R=1/2;
            pointStyle='*b-';
        case 9,
            M=2;
            R=3/4;
            pointStyle='^b-';
        case 12,
            M=4;
            R=1/2;
            pointStyle='*r-';
        case 18,            
            M=4;
            R=3/4;
            pointStyle='^r-';
        case 24,
            M=16;
            R=1/2;
            pointStyle='*g-';
        case 36,
            M=16;
            R=3/4;
            pointStyle='^g-';
        case 48,
            M=64;
            R=2/3;
            pointStyle='xk-';
        case 54,
            M=64;
            R=3/4;
            pointStyle='^k-';
        otherwise,
            disp(['Rate ' num2str(datarate) ' not supported']);
            exit(1);
    end
    switch R
    case 2/3,
        puncpat = [1,1,1,0];
    case 3/4,
        puncpat = [1,1,1,0,0,1];
    otherwise,
        puncpat = [];
    end
    disp(['Simulating Rate ' num2str(datarate) ' Mbps']);
    k = log2(M);
    snrs = EbNo_Vec + 10.*log10(k*R*48/64);
    
    for pak=1:totPak/k

        for ebno = 1:length(EbNo_Vec)
            bits = randi([0 1],k*packetSize,1);
            if ~isempty(puncpat)
                txVitBits = convenc(bits,trel,puncpat);
            else
                txVitBits = convenc(bits,trel);
            end
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


            %receive
            noisyTx = awgn(OFDMsig,snrs(ebno),'measured');
            rxOFDMsym = reshape(noisyTx,80,frameCount);
            rxFFTin = rxOFDMsym(17:80,:);
            rxFFTout = fft(rxFFTin,64,1);
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
            if ~isempty(puncpat)
                rxVitBits = vitdec(rxBits,trel,traceback,'cont','hard',puncpat); % Decode
            else
                rxVitBits = vitdec(rxBits,trel,traceback,'cont','hard'); % Decode
            end
            [~,berVec(pak,ebno)] = biterr(bits(1:end-traceback), rxVitBits(traceback+1:end));


        end
    end

    ber = mean(berVec(1:totPak/k,:),1);
    semilogy(EbNo_Vec, ber,pointStyle)
    if M == 2
        berTheory = berawgn(EbNo_Vec,'psk',2,'nondiff');
    else
        berTheory = berawgn(EbNo_Vec, 'qam', M);
    end
    hold on
    %semilogy(EbNo_Vec,berTheory,'b')
end
%6, 9, 12, 18, 24, 36, 48, 54
legend('BPSK R1/2 (6Mbps)','BPSK R3/4 (9Mbps)','4QAM R1/2 (12Mbps)',...
       '4QAM R3/4 (18Mbps)','16QAM R1/2 (24Mbps)','16QAM R3/4 (36Mbps)',...
       '64QAM R2/3 (48Mbps)','64QAM R3/4 (54Mbps)')
xlabel('EbNo');
ylabel('Bit Error Rate');
title('802.11a BER');
axis([0 20 1e-5 1e0])