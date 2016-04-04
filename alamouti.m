%% ALAMOUTI

clc, clear, close all;

%2 tx 2 rx

numbits = 1e5;
T = 1e-5;
dopshift = 200; %Hz
EbNo = 0:2:50;
EbNoHalfPow = EbNo - 10*log10(2); %2 transmitters at 1/2 power 

M = 2;
k = log2(M);

h0 = rayleighchan(T,dopshift);
h0.StorePathGains = true;

h1 = rayleighchan(T,dopshift);
h1.StorePathGains = true;

h2 = rayleighchan(T,dopshift);
h2.StorePathGains = true;

h3 = rayleighchan(T,dopshift);
h3.StorePathGains = true;

% Create an AWGNChannel and ErrorRate calculator System object
hAWGN = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
hErrorCalc_1rx = comm.ErrorRate;
hErrorCalc_2rx = comm.ErrorRate;

hMod = comm.BPSKModulator;       % Create a PSK modulator
hDemod = comm.BPSKDemodulator(); % Create a PSK demodulator
sMouti = randi([0 1],numbits,1);
s0 = sMouti(1:numbits/2);
s1 = sMouti(numbits/2+1:numbits);



s0Mod = step(hMod, s0);       % DPSK modulate the signal
s1Mod = step(hMod, s1);       % DPSK modulate the signal



s0DelayMod = -conj(step(hMod,s1));
s1DelayMod = conj(step(hMod,s0));

r0Faded = filter(h0,s0Mod);   % Apply the channel effects
r1Faded = filter(h1,s1Mod);
r2Faded = filter(h2,s0Mod);
r3Faded = filter(h3,s1Mod);

%Make sure alamouti signals have the same path gains as the original
%symbols

r0DelayFaded = h0.PathGains .* s0DelayMod;
r1DelayFaded = h1.PathGains .* s1DelayMod;
r2DelayFaded = h2.PathGains .* s0DelayMod;
r3DelayFaded = h3.PathGains .* s1DelayMod;


r0FadedComb = r0Faded + r1Faded;
r1FadedComb = r0DelayFaded + r1DelayFaded;
r2FadedComb = r2Faded + r3Faded;
r3FadedComb = r2DelayFaded + r3DelayFaded;

for n = 1:length(EbNo)
   hAWGN.SNR = EbNoHalfPow(n);
   r0AWGN = step(hAWGN,r0FadedComb);   % Add Gaussian noise
   r1AWGN = step(hAWGN,r1FadedComb);
   r2AWGN = step(hAWGN,r2FadedComb);   % Add Gaussian noise
   r3AWGN = step(hAWGN,r3FadedComb);
   
   h0pg = (h0.PathGains);
   h1pg = (h1.PathGains);
   h2pg = (h2.PathGains);
   h3pg = (h3.PathGains);
   
   sTilde0 = conj(h0pg).*r0AWGN + h1pg.*conj(r1AWGN);
   sTilde1 = conj(h1pg).*r0AWGN - h0pg.*conj(r1AWGN);
   
   sTilde0_2rx = conj(h0pg).*r0AWGN + h1pg.*conj(r1AWGN) + conj(h2pg).*r2AWGN + h3pg.*conj(r3AWGN);
   sTilde1_2rx = conj(h1pg).*r0AWGN - h0pg.*conj(r1AWGN) + conj(h3pg).*r2AWGN - h2pg.*conj(r3AWGN);
   
   
   rMouti0 = step(hDemod,sTilde0);
   rMouti1 = step(hDemod,sTilde1);
   
   rMouti0_2rx = step(hDemod,sTilde0_2rx);
   rMouti1_2rx = step(hDemod,sTilde1_2rx);
   
   %rMouti = [rMouti0;rMouti1];
   
   %r0AWGN = r0AWGN./h0.PathGains;
   
   
   
   %r0 = step(hDemod, r0AWGN);  % Demodulate
   %r1AWGN = r1AWGN./h1.PathGains;
   %r1 = step(hDemod, r1AWGN);
   
   

   
   reset(hErrorCalc_1rx);
   reset(hErrorCalc_2rx);
   % Compute error rate.
   berVec_1rx(:,n) = step(hErrorCalc_1rx,sMouti,[rMouti0;rMouti1]);
   berVec_2rx(:,n) = step(hErrorCalc_2rx,sMouti,[rMouti0_2rx;rMouti1_2rx]);
   %berVec1(:,n) = step(hErrorCalc1,s0,r1);
end

% Compute theoretical performance results, for comparison.
BERtheory = berfading(EbNoHalfPow,'psk',M,2);
figure;
% Plot BER results.
semilogy(EbNo,berVec_1rx(1,:),'g*',EbNo,berVec_2rx(1,:),'c*');
hold on;
semilogy(EbNo,BERtheory,'c-',EbNo,BERtheory,'g-');

legend('new scheme (2 Tx, 1 Rx)','new scheme (2 Tx, 2 Rx)');
xlabel('EbNo (dB)'); ylabel('BER');
title('Alamouti - Binary PSK over Rayleigh Fading Channel');
