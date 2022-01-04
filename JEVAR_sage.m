function [tauOut,xiOut,thetaOut] = JEVAR_sage(miu,IterMax,ys,lTilde,zcLen,r,thetaRange,xiRange,M,d,lambda)
% the SAGE method

Yset = zeros(2*zcLen,M,miu);
tauEst = zeros(miu,IterMax); xiEst = zeros(miu,IterMax); thetaEst = zeros(miu,IterMax);

for cc = 1:IterMax
    
    for uu = 1:miu
        Ysum = 0;
        for kk = 1:miu
            Ysum = Ysum + Yset(:,:,kk);
        end
        Ysum = Ysum - Yset(:,:,uu);
        Ysig = ys - Ysum;
        ys1 = Ysig(1:zcLen,:);
        ys2 = Ysig(zcLen+1:end,:);
        
        [ tauEst(uu,cc),xiEst(uu,cc),thetaEst(uu,cc),Ysigs] =  est_single_sage(ys1,ys2,Ysig,lTilde,zcLen,r,thetaRange,xiRange,M,d,lambda);
        
        Yset(:,:,uu) = Ysigs;
        
    end
    
    
    if miu>1
        if cc>1
            congncCond = abs(tauEst(:,cc)-tauEst(:,cc-1))<1e-4 & abs(xiEst(:,cc)-xiEst(:,cc-1))<1e-7 & rad2deg(abs(thetaEst(:,cc)-thetaEst(:,cc-1)))<5e-4;
            if  sum(congncCond) == miu
                break
            end
        end
    else
        break;
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

end

