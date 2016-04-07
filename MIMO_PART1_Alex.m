function [ mu_ZF,mu_MMSE,mu_PRECODING,mu_BASELINE ] = MIMO_PART1_Alex( H, SNR,str_type )

%http://www.mathworks.com/help/comm/examples/spatial-multiplexing.html
% channel side info 
%generate QPSK data
N = 2;                  % Number of transmit antennas
M = 2;                  % Number of receive antennas
 
modOrd = 2;             % constellation size = 2^modOrd
EbNoVec = SNR - 10*log10(modOrd);        % Eb/No in dB
%a local random stream to be used by random number generators for
% repeatability.
hStr = RandStream('mt19937ar');

% Create PSK modulator and demodulator System objects
hMod   = comm.PSKModulator(...
            'ModulationOrder',  2^modOrd, ...
            'PhaseOffset',      0, ...
            'BitInput',         true);
hDemod = comm.PSKDemodulator( ...
            'ModulationOrder',  2^modOrd, ...
            'PhaseOffset',      0, ...
            'BitOutput',        true);



Ntrials = 1e2 ; 
% Pre-allocate variables to store BER results for speed
BER_ZF =  zeros(length(EbNoVec), Ntrials) ;
BER_PRECODING =  zeros(length(EbNoVec), Ntrials) ;
BER_MMSE = zeros(length(EbNoVec), Ntrials) ;
BER_BASELINE = zeros(length(EbNoVec), Ntrials) ;



%%
%SNR = EbNoVec + 10*log10(modOrd);

tic;
% Loop over selected EbNo points
for idx = 1:length(EbNoVec)
    
    for p = 1:1:Ntrials    
        % Calculate SNR from EbNo for each independent transmission link
        snrIndB = SNR(idx);
        snrLinear = 10^(0.1*snrIndB);


        % Create random bit vector to modulate
        msg = randi(hStr, [0 1], [N*modOrd, 1]);

        % Modulate data
        txSig = step(hMod, msg);


        % Add noise to faded data   y = H* x + n
        n = 10^(-snrIndB/20)*(randn(M, 1) + j*randn(M, 1))/sqrt(2) ; 
        y  = H*txSig +  n ;

        rx_BASELINE = step(hDemod, y)  ;
        BER_BASELINE(idx,p) = biterr(rx_BASELINE,msg)/length(msg);

        %PRECODING (SVD)
        [U, SIGMA, V] = svd(H) ;
        x_tilde = V* txSig;


        y_tilde = SIGMA* x_tilde +  U'*n;

        %Demodulate data
        rxSig_PRECODING = step(hDemod, U*y_tilde);
        %BER
        BER_PRECODING(idx,p)  = biterr(rxSig_PRECODING,msg)/length(msg);




        %ZERO FORCING    
        W_ZF = inv(H'*H)*H' ; 
        %Demodulate data
        rxSig_ZF = step(hDemod, W_ZF*y);
        BER_ZF(idx,p)  = biterr(rxSig_ZF,msg)/length(msg);


        %MMSE 
        N_0 = 1; 
        W_MMSE = inv(H'*H + N_0)*H' ; 
        rxSig_MMSE = step(hDemod, W_MMSE*y);   
        BER_MMSE(idx,p)  = biterr(rxSig_MMSE,msg)/length(msg) ;  

     
    
    end    
end    

toc;
    %%
       
% Set up a figure for visualizing BER results
figure;
grid on;
xlim([SNR(1)-0.01, SNR(end)]);
ylim([ 0 1 ])
xlabel('SNR (dB)');
ylabel('AVERAGE BER');
title(['2x2 Uncoded QPSK System',' ' , str_type]);  
hold on ;
plot(SNR(:), mean(BER_ZF,2), 'r*') 
plot(SNR(:), mean(BER_MMSE,2), 'b*') 
plot(SNR(:), mean(BER_PRECODING,2), 'mo') 
plot(SNR(:), mean(BER_BASELINE,2), 'k-','Markersize',12);
legend('ZF', 'MMSE', 'PRECODING','NO SCHEME');
hold off;
%%



mu_ZF = mean(BER_ZF,2) ;
mu_MMSE = mean(BER_MMSE,2);
mu_PRECODING = mean(BER_PRECODING,2);
mu_BASELINE = mean(BER_BASELINE,2);
      


end

