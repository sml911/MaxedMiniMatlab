function [ berzf,bermm,bersvd ] = mimo_fn( snr,n )

M = 2; 
rng(20)
msgD = randi([0 M-1],3,n);
msgQAM = qammod(msgD,M);

chan11 = rayleighchan(1e-5,0);
chan12 = rayleighchan(1e-5,0);
chan13 = rayleighchan(1e-5,0);
chan23 = rayleighchan(1e-5,0);
chan22 = rayleighchan(1e-5,0);
chan33 = rayleighchan(1e-5,0);
chan21 = rayleighchan(1e-5,0);
chan31 = rayleighchan(1e-5,0);
chan32 = rayleighchan(1e-5,0);

H = [chan11.PathGains, chan12.PathGains, chan13.PathGains;
     chan21.PathGains, chan22.PathGains, chan23.PathGains;
     chan22.PathGains, chan23.PathGains, chan33.PathGains];
 
[U,S,V] = svd(H);
 
N0 = sum(abs(msgQAM(:)).^2)./length(msgQAM(:))./10^(snr/10)

 for ii = 1:n
    y(:,ii) =  H*msgQAM(:,ii);
    ysvd(:,ii) =  H*V*msgQAM(:,ii);
 end
 for ii = 1:3
     z(ii,:) = y(ii,:) +randn(size(y(ii,:))).*sqrt(N0);
     zsvd(ii,:) = ysvd(ii,:) +randn(size(ysvd(ii,:))).*sqrt(N0);
 end
 for ii = 1:n
    zf(:,ii) =  inv(H'*H)*H'*z(:,ii);
    mmse(:,ii) = inv(H'*H + N0)*H'*z(:,ii);
    svdp(:,ii) = inv(S)*U'*zsvd(:,ii);
 end
 
 msgunczf = qamdemod(zf,M);
 msguncmm = qamdemod(mmse,M);
 msguncsvd = qamdemod(svdp,M);

 berzf = biterr(msgunczf,msgD)./3./n;
 bermm = biterr(msguncmm,msgD)./3./n;
 bersvd = biterr(msguncsvd,msgD)./3./n;
 end

