function [tau, xi, theta] = coarse_est_2Dfft(ys1,ys2,thetaRange,lTilde,r,fftLen1,fftLen2)
% coarse estimation by 2D-fft

xiMtx = fft(conj(ys1),fftLen1,1);
thetaxiMatx(:,:) = abs(fft(xiMtx, fftLen2, 2));

%%
[~,indx] = max(thetaxiMatx(:));

[xi1Indx, thetaIndx] = ind2sub([fftLen1,fftLen2],indx);
thetaPhase = (thetaIndx-1)/fftLen2;
if thetaPhase > 1/2
    thetaPhase = thetaPhase - 1;
end
theta = asin(2*thetaPhase);
atheta = exp(thetaRange*sin(theta));
xiVec =  atheta.'*ys2';

tmp = abs(fft(xiVec,fftLen1));
[~,xi2Indx] = max(tmp);

%%
fx1 = (xi1Indx-1)/fftLen1;  %k/N
fx1 = -fx1;
if fx1 < -1/2
    fx1 = fx1 + 1;
end
xi1 =  fx1;
fx2 = (xi2Indx-1)/fftLen1;  %k/N
fx2 = -fx2;
if fx2 < -1/2
    fx2 = fx2 + 1;
end
xi2 =  fx2;

xi = (xi1+xi2)/2;
tau = (xi1-xi2)/(-2*r/lTilde);
%%
% zcLen = length(ys1);
% thetaRange = -1j*2*pi*1/2*(0:6-1).';
% aCompen = exp(thetaRange*sin(theta));
% figure
% Contourmap_1211( 1024,ys1*conj(aCompen),ys2*conj(aCompen),r,lTilde,zcLen)
% hold on, plot(0,xi1,'ro','MarkerFaceColor','r','markersize',8)
% hold on, plot(0,xi2,'ro','MarkerFaceColor','r','markersize',8)
% hold on, plot(tauEst,xi,'ro','MarkerFaceColor','k','markersize',8)

end
