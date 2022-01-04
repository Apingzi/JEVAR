function [ tauOut,xiOut,thetaOut,yRefactorOut] =  est_single(ys1,ys2,ys,lTilde,zcLen,r,thetaRange,xiRange,M,d,lambda)
% estimate the parameters of a path in single-case

fftLen1 = 1024;                  % temporal domain 
fftLen2 = 128;                   % spatial domain  
%% coarse estimation using 2D-fft
[tauEst(1) ,xiEst(1),thetaEst(1)] = coarse_est_2Dfft(ys1,ys2,thetaRange,lTilde,r,fftLen1,fftLen2);

%% refine the estimates using the Newton's method
ysVec = reshape(ys.',[length(thetaRange)*length(xiRange),1]);

[ tauEst(2),  xiEst(2),  thetaEst(2) ] = fine_est_Newton(r,lTilde,zcLen,ysVec,xiEst(1),thetaEst(1),tauEst(1),thetaRange,xiRange,d,lambda,M);

%% the estimate of the JEVAR parameters
tauOut = tauEst(2);
xiOut = xiEst(2);
thetaOut = thetaEst(2);

%% signal reconstruction
atheta = exp(thetaRange*sin(thetaOut));
acfo = exp(xiOut*xiRange);
zc1 = exp(1j*pi*r/lTilde*tauOut*(tauOut+zcLen))*exp(-2j*pi*r/lTilde*tauOut*(0:zcLen-1));
b = [zc1 conj(zc1)].';
yRefactorOut = kron((acfo.*b),atheta);
end


