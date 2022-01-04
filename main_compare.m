% Code for joint estimation of frequency offsets, DOAs and time delays of multipath signal. 
% The detailed algorithm discreption can be seen in the paper --- Z. Yang, R. Wang, Y. Jiang and J. Li, "Joint Estimation of Velocity, Angle-of-Arrival and Range (JEVAR) Using a Conjugate Pair of Zadoff-Chu Sequences," in IEEE Transactions on Signal Processing, vol. 69, pp. 6009-6022, 2021, doi: 10.1109/TSP.2021.3122907.
close all;
clc;clear;
% profile on;
%% parameters
zcLen = 150;                                        % the effective length for the first/seconde half pilot
lTilde = 400;
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
M = 45;                                             % number of antennas
miu = 30;                                           % number of multipaths
IterMax_ap = 200;                                    % maximum number of itrerations (AP)
IterMax_sage = 200;                                  % maximum number of itrerations (SAGE)
numMonte = 200;                                      % number of MonteCarlo trials
proInt =  0;
BW = 20e6; c = 3e8; fc = 2.4e9; lambda = 3e8/2.4e9; d = lambda/2;               %BW is the bandwidth; c is the light speed; fs is the carrier frequency;lambda is the carrier wavelength; d is the inter-antenna spacing
thetaRange = -1j*2*pi*d/lambda*(0:M-1).';
xiRange = 1j*2*pi*[(0:zcLen-1) (zcLen+giLen:2*zcLen+giLen-1)].';
% xiRange = 1j*2*pi*[(giLen:giLen+zcLen-1) (giLen*3+zcLen:giLen*3+2*zcLen-1)].';
snrDbSet = 0:5:25;                                 %SNR in dB
% snrDbSet = 15;
snrLen = length(snrDbSet);
tauErr_ap = zeros(numMonte,snrLen); xiErr_ap = zeros(numMonte,snrLen); thetaErr_ap = zeros(numMonte,snrLen);
tauErr_sage = zeros(numMonte,snrLen); xiErr_sage = zeros(numMonte,snrLen); thetaErr_sage = zeros(numMonte,snrLen);
 
%% Monte-Carlo trials for the comparion between the alternating projection (AP) method and the SAGE method 
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
        [tauEstList_ap,xiEstList_ap,thetaEstList_ap] = JEVAR_ap(miu,IterMax_ap,ys1,ys2,ys,lTilde,zcLen,r,thetaRange,xiRange,M,d,lambda);
        [tauEstList_sage,xiEstList_sage,thetaEstList_sage] = JEVAR_sage(miu,IterMax_sage,ys,lTilde,zcLen,r,thetaRange,xiRange,M,d,lambda);

        tauErr_ap(mm,ss) =  tauEstList_ap(1) - tauList(1)-proInt;
        xiErr_ap(mm,ss) =  xiEstList_ap(1) - xiList(1);
        thetaErr_ap(mm,ss) = thetaEstList_ap(1) - thetaList(1);
        
        tauErr_sage(mm,ss) =  tauEstList_sage(1) - tauList(1)-proInt;
        xiErr_sage(mm,ss) =  xiEstList_sage(1) - xiList(1);
        thetaErr_sage(mm,ss) = thetaEstList_sage(1) - thetaList(1);
        
    end
    
%     if mod(mm,10)==0
        disp(mm)
%     end
end
%%

tauMse_ap = mean(tauErr_ap(1:mm,:).^2,1);
xiMse_ap = mean(xiErr_ap(1:mm,:).^2,1);
thetaMse_ap =  mean(thetaErr_ap(1:mm,:).^2,1);

tauMse_sage = mean(tauErr_sage(1:mm,:).^2,1);
xiMse_sage = mean(xiErr_sage(1:mm,:).^2,1);
thetaMse_sage =  mean(thetaErr_sage(1:mm,:).^2,1);

% CRB & figures
athetaList = exp(thetaRange*sin(thetaList));
[crbTau,crbXi,crbTheta,crbHr,crbHi,crbH] = crb_cal(lTilde,r,miu,M,tauList,xiList,thetaList,betaList,athetaList,snrDbSet,zcLen,giLen,d,lambda,nRange);
% [crbTauR,crbCfo,crbTheta,crbHr,crbHi,crbH] = crb_compute_DOA_multipath1(miu,sigLen1,M,tauList,xiList,thetaList,hList,athetaList,supSamp,rcfCoeffList,roCoeff,snrDbSet,upSampRate,zcLen,giLen,oneSideWinLen,d,lambda,proInt);


figure;
% subplot(1,3,1)
semilogy(snrDbSet,sqrt(crbTheta(1,:).')/pi*180,'k');hold on
semilogy(snrDbSet,sqrt(thetaMse_ap)/pi*180,'*'); hold on;
semilogy(snrDbSet,sqrt(thetaMse_sage)/pi*180,'o'); hold on;
xlabel('SNR (dB)');
ylabel('RMSE (\circ)')
legend('CRB','AP','SAGE');
title('Angle Estimation of the First Path')
grid on;

figure;
% subplot(1,3,2)
semilogy(snrDbSet,sqrt(crbXi(1,:).'),'k'); hold on
semilogy(snrDbSet,sqrt(xiMse_ap),'*'); hold on;
semilogy(snrDbSet,sqrt(xiMse_sage),'o'); hold on;
xlabel('SNR (dB)');
ylabel('RMSE (1/Ts)')
legend('CRB','AP','SAGE');
title('Frequency Offset Estimation of the First Path')
grid on;

figure;
% subplot(1,3,3)
semilogy(snrDbSet,sqrt(crbTau(1,:).'),'k'); hold on
semilogy(snrDbSet,sqrt(tauMse_ap),'*'); hold on
semilogy(snrDbSet,sqrt(tauMse_sage),'o'); hold on;
xlabel('SNR (dB)');
ylabel('RMSE (Ts)')
legend('CRB','AP','SAGE');
title('Time Delay Estimation of the First Path')
grid on;

% toc;
% profile off;
% profile viewer;