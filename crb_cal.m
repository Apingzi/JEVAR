function [crbTau,crbCfo,crbtheta,crbbetar,crbbetai,crbbeta] = crb_cal(lTilde,r,miu,M,tau,xi,theta,beta,atheta,snrDbSet,zcLen,giLen,d,lambda,nRange)
%% calculate CRBs for each parameters
% param = [hr,hi,tauR,cfo,theta].';
crbTau = zeros(miu,length(snrDbSet));
crbCfo = zeros(miu,length(snrDbSet));
crbtheta = zeros(miu,length(snrDbSet));
paramEstNum = 5;
% y = zeros(M*L,paramEstNum*K);
for kk = 1:miu
%     zcwithTau =  exp(1j*pi*r/lTilde*(nRange-zcLen/2-tauR(kk)).^2);
%     zcPilot = [zcwithTau, conj(zcwithTau)];
    zc1 = exp(1j*pi*r/lTilde*(nRange-zcLen/2-tau(kk)).^2);
    zcPilot = [zc1, conj(zc1)].';
    sd1tau1 = (1j*pi*r/lTilde*(2*tau(kk)+zcLen-2*nRange)).*zc1;
    sd1tau = [sd1tau1 conj(sd1tau1)].';
    dxi = exp(1j*2*pi*xi(kk)*(0:length(nRange)*2-1)).';
    y1 = beta(kk)*kron(sd1tau.*dxi,atheta(:,kk));
    athetaPartial = -1j*2*pi*d*cos(theta(kk))/lambda*(0:M-1).'.*atheta(:,kk);
    y2 = beta(kk)*kron(zcPilot.*dxi,athetaPartial);
    dxid1 = (1j*2*pi*(0:length(nRange)*2-1)).'.*dxi;
    y3 = beta(kk)*kron(zcPilot.*dxid1,atheta(:,kk));
    y4 = kron(zcPilot.*dxi,atheta(:,kk));
    y5 = 1j*y4;
    y(:,(kk-1)*paramEstNum+(1:paramEstNum)) = [y1,y2,y3,y4,y5];
end


rxZc = zeros(zcLen*2*M,paramEstNum*miu);
% rxZc = zeros(zcLen*upSampRate*M,paramEstNum*K);
rx = y((1:1:(zcLen+giLen)*2*M),:);
rxZc(1:zcLen*M,:) = rx((giLen/2)*M+(1:zcLen*M),:);
rxZc(zcLen*M+1:zcLen*2*M,:) = rx((3/2*giLen+zcLen)*M+(1:zcLen*M),:);

%% fisher 
Fisher = real(rxZc'*rxZc)*2;
% tmp = diag(Fisher);
% Fisher = zeros(K,K);
% Fisher = diag(tmp);
FisherInv = inv(Fisher);
for kk = 1:miu
    crbTau(kk,:) = FisherInv((kk-1)*paramEstNum+1,(kk-1)*paramEstNum+1).*db2pow(-snrDbSet);     % CRB of time delays
    crbtheta(kk,:) = FisherInv((kk-1)*paramEstNum+2,(kk-1)*paramEstNum+2).*db2pow(-snrDbSet);   % CRB of DOAs
    crbCfo(kk,:) = FisherInv((kk-1)*paramEstNum+3,(kk-1)*paramEstNum+3).*db2pow(-snrDbSet);     % CRB of frequency offsets
    crbbetar(kk,:) = FisherInv((kk-1)*paramEstNum+4,(kk-1)*paramEstNum+4).*db2pow(-snrDbSet);      % CRB of the real part of channel gains
    crbbetai(kk,:) = FisherInv((kk-1)*paramEstNum+5,(kk-1)*paramEstNum+5).*db2pow(-snrDbSet);      % CRB of the imaginary part of channel gains
    crbbeta(kk,:) = crbbetar(kk,:).^2+crbbetai(kk,:).^2;
end

end

