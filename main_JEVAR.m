close all;
clc;clear;
%% 
zcLen = 200;                                        % the effective length for the first/seconde half pilot
lTilde = 500;
r = 1;
giLen = 40;                                         % the length of prefix+suffix
nRange = -giLen/2:zcLen+giLen/2-1;
zcRoot =  exp(1j*pi*r/lTilde*(nRange-zcLen/2).^2);
temRange = 0:zcLen-1;
zcSeq = exp(1j*pi*r/lTilde*(temRange-zcLen/2).^2);  % the first half pilot  

seqSend = [zcRoot, conj(zcRoot)];                   % transmitted pilot
roCoeff = 0.3;                                      % roll-off factor    
upSampRate = 1;
oneSideWinLen = 20;                                 % half length of the raised-cosine filter
M = 40;                                             % number of antennas
IterMax = 30;                                       % maximum number of itrerations 
numMonte = 300;                                      % number of MonteCarlo trials
BW = 20e6; c = 3e8; fc = 2.4e9; lambda = 3e8/2.4e9; d = lambda/2;               %BW is the bandwidth; c is the light speed; fs is the carrier frequency;lambda is the carrier wavelength; d is the inter-antenna spacing
thetaRange = -1j*2*pi*d/lambda*(0:M-1).';
xiRange = 1j*2*pi*[(0:zcLen-1) (zcLen+giLen:2*zcLen+giLen-1)].';
% xiRange = 1j*2*pi*[(giLen:giLen+zcLen-1) (giLen*3+zcLen:giLen*3+2*zcLen-1)].';
miu = 5;                                           % number of multipaths
proInt =  0;
%%
snrDbSet = -5:5:30;                                 %SNR in dB
% snrDbSet = 30;
snrLen = length(snrDbSet);
tauErr = zeros(numMonte,snrLen); xiErr = zeros(numMonte,snrLen); thetaErr = zeros(numMonte,snrLen);
%% 
supSamp = kron(seqSend,[1, zeros(1,upSampRate-1)]).';
sigLen = 2*oneSideWinLen*upSampRate+1+length(supSamp)-1;
[rx,tauList,betaList,thetaList,xiList] = gen_multipathSig(miu,supSamp,upSampRate,roCoeff,oneSideWinLen,sigLen,M,d,lambda);            % generate multipath signal
for mm = 1:numMonte
     
    for ss = 1:snrLen
        snrDb = snrDbSet(ss);
        y = rx + 1*(randn(M,sigLen) + 1j*randn(M,sigLen))./sqrt(2)/db2mag(snrDb);
        sampleSeq = y(:,oneSideWinLen*upSampRate+(1:upSampRate:(zcLen+giLen)*2*upSampRate));
        
        y1 = sampleSeq(:,giLen/2+(1:zcLen)-proInt);    
        y2 = sampleSeq(:,giLen*3/2+zcLen+(1:zcLen)-proInt);
        ys1 = (y1.*conj(zcSeq)).';
        ys2 = (y2.*zcSeq).';
        ys = [ys1.' ys2.'].';
        [tauEstList,xiEstList,thetaEstList] = JEVAR_ap(miu,IterMax,ys1,ys2,ys,lTilde,zcLen,r,thetaRange,xiRange,M,d,lambda);

        tauErr(mm,ss) =  tauEstList(1) - tauList(1)-proInt;
        xiErr(mm,ss) =  xiEstList(1) - xiList(1);
        thetaErr(mm,ss) = thetaEstList(1) - thetaList(1);
        
    end
    
    if mod(mm,10)==0
        disp(mm)
    end
end
%%

tauMse = mean(tauErr(1:mm,:).^2,1);
xiMse = mean(xiErr(1:mm,:).^2,1);
thetaMse=  mean(thetaErr(1:mm,:).^2,1);

% CRB & figures
athetaList = exp(thetaRange*sin(thetaList));
[crbTau,crbXi,crbTheta,crbHr,crbHi,crbH] = crb_cal(lTilde,r,miu,M,tauList,xiList,thetaList,betaList,athetaList,snrDbSet,zcLen,giLen,d,lambda,nRange);
% [crbTauR,crbCfo,crbTheta,crbHr,crbHi,crbH] = crb_compute_DOA_multipath1(miu,sigLen1,M,tauList,xiList,thetaList,hList,athetaList,supSamp,rcfCoeffList,roCoeff,snrDbSet,upSampRate,zcLen,giLen,oneSideWinLen,d,lambda,proInt);

figure
semilogy(snrDbSet,sqrt(crbTheta(1,:))/pi*180,'LineWidth',1,'MarkerSize',8);hold on
semilogy(snrDbSet,sqrt(thetaMse)/pi*180,'o','LineWidth',1,'MarkerSize',8)
% plot(SNR,multiUserEPC(1,:),'g-+','LineWidth',2,'MarkerSize',10);
xlabel('SNR(dB)');
ylabel('RMSE of AoA Estimation (\circ)')
legend('CRB','the JEVAR scheme')
title('AOA Estimation of the First Path');
grid on;

figure
semilogy(snrDbSet,sqrt(crbXi(1,:))*BW*c/fc,'LineWidth',1,'MarkerSize',8);hold on
semilogy(snrDbSet,sqrt(xiMse)*BW*c/fc,'o','LineWidth',1,'MarkerSize',8)
xlabel('SNR(dB)');
ylabel('RMSE of Velocity Estimation (m/s)')
legend('CRB','the JEVAR scheme')
title('Velocity Estimation of the First Path');
grid on;

figure
semilogy(snrDbSet,sqrt(crbTau(1,:))*c/BW,'LineWidth',1,'MarkerSize',8);hold on
semilogy(snrDbSet,sqrt(tauMse)*c/BW,'o','LineWidth',1,'MarkerSize',8);
xlabel('SNR(dB)');
ylabel('RMSE of Range Estimation (m)')
legend('CRB','the JEVAR scheme')
title('Range Estimation of the First Path');
grid on;
