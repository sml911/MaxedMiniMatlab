function [ OFDMsig ] = OFDMmod( txMod,frameCount )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
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
    OFDMsig = OFDMsig/sqrt(var(OFDMsig));

end

