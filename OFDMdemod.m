function [ rxDemod ] = OFDMdemod( rxFFT,frameCount )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    rxOFDM = [rxFFT(39:64,:);rxFFT(1:27,:)];
    rxtones = [rxOFDM(1:5,:);
               rxOFDM(7:19,:);
               rxOFDM(21:26,:);
               rxOFDM(28:33,:);
               rxOFDM(35:47,:);
               rxOFDM(49:53,:)];

    rxDemod = reshape(rxtones,frameCount*48,1);
end

