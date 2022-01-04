function [ tauOut,xiOut,thetaOut,yRefactorOut] =  est_projection(ysTransform,lTilde,zcLen,r,thetaRange,xiRange,M,yRefactor,d,lambda)
% estimate the parameters of a path in multipath case

fftLen1 = 1024;                                             % temporal domain 
fftLen2 = 128;                                               % spatial domain

%% 
yRefactorhYs = yRefactor'*ysTransform;
% LambdahYs = ysTransform - yRefactor*inv(yRefactor'*yRefactor)*yRefactorhYs;
inter = yRefactor/(yRefactor'*yRefactor);
LambdahYs = ysTransform - inter*yRefactorhYs;
ysAfterLam = reshape(LambdahYs,[M,zcLen*2]).';
ys1AfterLam = ysAfterLam(1:zcLen,:);
ys2AfterLam = ysAfterLam(zcLen+1:end,:);

%% coarse estimation using 2D-fft
% [thetaE1 ,freq,tauEst] = twoDFourier1009New(ys1AfterLam,ys2AfterLam,thetaRange,lTilde,zcLen,r,roCoeff,T,fftLen1,fftLen2,proInt);
[tauEst(1) ,xiEst(1),thetaEst(1)] = coarse_est_2Dfft(ys1AfterLam,ys2AfterLam,thetaRange,lTilde,r,fftLen1,fftLen2);


%% refine the estimates using the Newton's method

[tauEst(2), xiEst(2), thetaEst(2)] = fine_est_Newton_projection(yRefactor,LambdahYs,inter,r,lTilde,zcLen,xiEst(1),thetaEst(1),tauEst(1) ,thetaRange,xiRange,d,lambda,M);

%% the estimate of the JEVAR parameters
tauOut = tauEst(2);
xiOut = xiEst(2);
thetaOut = thetaEst(2);
%% signal reconstruction
atheta = exp(thetaRange*sin(thetaOut));
dxi = exp(xiOut*xiRange);
zc1 = exp(1j*pi*r/lTilde*tauOut*(tauOut+zcLen))*exp(-2j*pi*r/lTilde*tauOut*(0:zcLen-1));
b = [zc1 conj(zc1)].';
yRefactorOut = kron((dxi.*b),atheta);
end


