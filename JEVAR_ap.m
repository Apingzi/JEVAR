function [tauOut,xiOut,thetaOut] = JEVAR_ap(miu,IterMax,ys1,ys2,ys,lTilde,zcLen,r,thetaRange,xiRange,M,d,lambda)
% the JEVAR scheme
% inputs - miu: number of multipath,             IterMax: maximum iterations of the alternating projection
%          ys1: the former ZC sequence,         ys1: the latter ZC sequence
%          ys: the pair of conjugate sequences, lTilde：cf. the JEVAR paper
%          zcLen: length of a ZC sequence       r: index of a ZC sequence
%          thetaRange: indics range of the DOA vector  
%          xiRange: indics range of the frequency offset vector
%          M: number of the antennas at the receiver 
%          d: inter-antenna spacing             lambda: carrier wavelength
% output - tauOutPut:   delays of the all the multipath
%          xiOutPut:    frequency offsets of the all the multipath
%          thetaOutPut: DOAs of the all the multipath
% programmer - Zhiyu Yang
% copyright - CSRL@Fudan,2021/05/10

ysTransform = reshape(ys.',[M*zcLen*2,1]);
tauEst = zeros(miu,IterMax); xiEst = zeros(miu,IterMax); thetaEst = zeros(miu,IterMax);
%% The first round iteration of the AP
% apMtx = []; % part of the construction of orthogonal projection matrix
[ tauEst(1,1),xiEst(1,1),thetaEst(1,1),yRefactor1] =  est_single(ys1,ys2,ys,lTilde,zcLen,r,thetaRange,xiRange,M,d,lambda);
apMtx = yRefactor1;
cc = 1;
% s = apMtx;
% hHat =  (s'*s)\s'*ysTransform;
if miu > 1
    for uu = 2:miu
        [ tauEst(uu,1),xiEst(uu,1),thetaEst(uu,1),yRefactor2] =  est_projection(ysTransform,lTilde,zcLen,r,thetaRange,xiRange,M,apMtx,d,lambda);
        apMtx = [apMtx yRefactor2];
        
    end
%     s = apMtx;
%     betaHat =  (s'*s)\s'*ysTransform;
%     cosfunc(1) = norm(ysTransform - s*hHat);      % cosfunction of the objective function
%% The following iterations of the AP
    for cc = 2:IterMax
        for uu = 1:miu
             tmpMtx = apMtx;
             tmpMtx(:,uu )= [];
            [ tauEst(uu,cc),xiEst(uu,cc),thetaEst(uu,cc),yRefactor2] =  est_projection(ysTransform,lTilde,zcLen,r,thetaRange,xiRange,M,tmpMtx,d,lambda);
            apMtx(:,uu) = yRefactor2;
            
        end
%         s = apMtx;
%         betaHat =  (s'*s)\s'*ysTransform;
%         cosfunc(cc) = norm(ysTransform - s*hHat);
            congncCond = abs(tauEst(:,cc)-tauEst(:,cc-1))<1e-4 & abs(xiEst(:,cc)-xiEst(:,cc-1))<1e-7 & rad2deg(abs(thetaEst(:,cc)-thetaEst(:,cc-1)))<5e-4;
        if  sum(congncCond) == miu
            break;
        end
    end
end
%%
tau = tauEst(:,cc);
xi = xiEst(:,cc);
theta = thetaEst(:,cc);
% [hList, indx] = sort(hHat,'descend');  % the path with the maximum amplitude is the first path
[~, indx] = sort(tau);               % the path with the minimum delay is the first path
tauOut = tau(indx);
xiOut = xi(indx);
thetaOut = theta(indx);
% figure;plot(cosfunc);
end

