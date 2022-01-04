function [rx,tauListOut,betaListOut,thetaListOut,xiListOut] = gen_multipathSig(miu,supSamp,upSampRate,roCoeff,oneSideWinLen,rxSigLen,M,d,lambda)
% generate multipath signal
% inputs - miu: number of multipaths,             
%          supSamp: the pilot sequence,         
%          upSampRate: Upsampling rate,
%          roCoeff: roll-off factor          
%          oneSideWinLen: half length of the raised-cosine filter 
%          rxSigLen: length of the filter convolved with the pilot
%          M: number of the antennas at the receiver 
%          d: inter-antenna spacing             lambda: carrier wavelength
% output - rx: multipath signal
%          tauListOut:  delays of the all the multipaths
%          betaListOut: complex gains of the all the multipaths
%          thetaListOut: DOAs of the all the multipaths
%          xiOutPut:    frequency offsets of the all the multipaths
%          thetaOutPut: DOAs of the all the multipath
% programmer - Zhiyu Yang
% copyright - CSRL@Fudan,2021/05/10

%% JEVAR parameters
tauList = linspace(0.5,13.5,miu);       % time delays
thetaList = deg2rad( kron(linspace(-75,-5,15),[-1 1]) );        % DOAs
% thetaList = deg2rad(linspace(-80,80,miu)  );        % DOAs

betaList = exp(1j*2*pi*rand(1,miu));                          % complex channel gains
% betaList = [0.384162239326054 + 0.923265603104541i,0.439554671349985 + 0.898215837588275i,0.362106337926951 - 0.932136792554147i,-0.845767315702287 + 0.533551916583333i,-0.168582943620349 + 0.985687471321564i,-0.984635793286984 + 0.174620601814654i,-0.636080808114558 + 0.771622450132402i,0.987420962938360 + 0.158113383209275i,-0.882727058270771 - 0.469886093214762i,0.892615690307231 + 0.450818399599380i,0.949105967438074 + 0.314956921774134i,-0.527689454403093 - 0.849437366561988i,0.957152238572149 - 0.289585207143467i,-0.459158239323490 + 0.888354496392827i,0.898487816457061 - 0.438998455211658i,-0.867324133693183 + 0.497743756478542i,-0.210517487319667 - 0.977590091772934i,-0.410456461983741 - 0.911880196525722i,-0.932724286602451 + 0.360590356473866i,-0.745156079400045 + 0.666890108888379i,0.624291075504978 + 0.781191815782038i,0.569795161333861 + 0.821786757085145i,-0.174627309048289 - 0.984634603766572i,0.869378023583309 - 0.494147601542677i,-0.806587243524396 + 0.591115063742853i,0.639299915623908 - 0.768957487695688i,0.267860561297041 + 0.963457689627123i,0.286636893265007 - 0.958039295341994i,0.838172273840722 - 0.545405573279806i,-0.161074533726204 - 0.986942244807104i];
xiList = (unidrnd(9,1,miu)-5).*10.^(-unidrnd(6,1,miu));       % frequency offsets
% xiList = [0,-2.00000000000000e-06,-0.000300000000000000,0.0100000000000000,0,2.00000000000000e-06,0.00400000000000000,0.200000000000000,2.00000000000000e-05,-0.0200000000000000,0,1.00000000000000e-05,-0.0300000000000000,-0.00100000000000000,-4.00000000000000e-05,2.00000000000000e-06,-0.100000000000000,0.0100000000000000,4.00000000000000e-06,-0.0300000000000000,4.00000000000000e-06,-4.00000000000000e-05,-4.00000000000000e-05,0.400000000000000,0.0400000000000000,-0.00100000000000000,-4.00000000000000e-05,0.00200000000000000,0,3.00000000000000e-05];
%%
rxSigwithTheta = zeros(M,rxSigLen,miu);
rxSigwithTauXi = zeros(rxSigLen,miu);
for ii = 1:miu
    atheta = exp(-1j*2*pi*d*sin(thetaList(ii))/lambda*(0:M-1)).';
    rcfCoeff = raised_cosine(upSampRate,roCoeff,tauList(ii),oneSideWinLen);
    rxSigwithTauXi(:,ii) = betaList(ii)*conv(rcfCoeff,supSamp).*exp(1j*2*pi*xiList(ii)/upSampRate*(0:rxSigLen-1)).';
    rxSigwithTheta(:,:,ii) = atheta*rxSigwithTauXi(:,ii).';
end

rx = sum(rxSigwithTheta,3);
tauListOut = tauList(1:miu);
betaListOut = betaList(1:miu);
thetaListOut = thetaList(1:miu);
xiListOut = xiList(1:miu);
end
