function [ tauEst, xiEst, thetaEst ] = fine_est_Newton(r,lTilde,zcLen,ysVec,xi,theta,tau,thetaRange,xiRange,d,lambda,M)
% use the Newton¡¯s iteration method to  search for a fine estimate (with backtracking line search)
alpha = 0.3;
beta = 0.7;
minStepSize = 1e-8;
maxIterNewton = 20;
maxIterBLS = 100;

temRange = 0:zcLen-1;

%%  the Newton's iterations
Value = [tau; xi; theta];
for nn = 1:maxIterNewton
    
    zc1 = exp(1j*pi*r/lTilde*tau*(tau+zcLen))*exp(-2j*pi*r/lTilde*tau*temRange);
    b = [zc1 conj(zc1)].';
    dxi = exp(xi*xiRange);
    xib = dxi.*b;
    atheta = exp(thetaRange*sin(theta));
    s = kron(xib,atheta);
    yhS = ysVec'*s;
    B = norm(yhS)^2;
    
    %% 1.The first/second derivative of delay
    sd1tau1 = (1j*pi*r/lTilde*(2*tau+zcLen-2*temRange)).*zc1;
    sd1tau = [sd1tau1 conj(sd1tau1)].';
    sd2tau1 = (1j*pi*r/lTilde*2 + (1j*pi*r/lTilde*(2*tau+zcLen-2*temRange)).^2 ).*zc1;
    sd2tau = [sd2tau1 conj(sd2tau1)].';
    
    
    s1tau = kron(dxi.*sd1tau,atheta);
    s2tau = kron(dxi.*sd2tau,atheta);
    s1tauhYs = s1tau'*ysVec;
    s2tauhYs = s2tau'*ysVec;
    B1dtau = 2*real(s1tauhYs*yhS);

    d1tau = B1dtau;
    B2dtau = 2*norm(ysVec'*s1tau)^2 + 2*real(s2tauhYs*yhS);

    d2tau = B2dtau;
    %% 2.The first/second derivative of frequency offset
    
    dxid1 = xiRange.*dxi;
    dxid2 = xiRange .* dxid1;
    s1xi = kron((dxid1.*b),atheta);
    s2xi = kron((dxid2.*b),atheta);
    s1xihYs = s1xi'*ysVec;
    s2xihYs = s2xi'*ysVec;

    B1dxi = 2*real(s1xihYs*yhS);

    d1xi = B1dxi;
    B2dxi =  2*norm(ysVec'*s1xi)^2 + 2*real(s2xihYs*yhS);

    d2xi = B2dxi;
    %% 3.The first/second derivative of DOA
    
    athetad1 = thetaRange*cos(theta).*atheta;
    athetad2 = ((1j*2*pi*d*sin(theta)/lambda*(0:M-1)).'-(4*pi^2*d^2*cos(theta)^2/lambda^2*(0:M-1).^2).').*atheta;
    s1theta = kron(xib,athetad1);
    s2theta = kron(xib,athetad2);
    s1thetahYs = s1theta'*ysVec;
    s2thetahYs = s2theta'*ysVec;
    
    B1dtheta = 2*real(s1thetahYs*yhS);

    d1theta = B1dtheta;
    B2dtheta =  2*norm(ysVec'*s1theta)^2 + 2*real(s2thetahYs*yhS);

    d2theta = B2dtheta;
    %% 4.The second partial derivative of delay and frequency offset
    s2dtauxi =  kron(dxid1.*sd1tau,atheta);
    s2dtauxihYs = s2dtauxi'*ysVec;
    B2dtauxi = 2*real(s2dtauxihYs*yhS + s1tauhYs*s1xihYs');

    d2tauxi = B2dtauxi;
    %% 5.The second partial derivative of delay and theta
    
    s2dtautheta =  kron(dxi.*sd1tau,athetad1);
    s2dtauthetahYs = s2dtautheta'*ysVec;
    B2dtautheta = 2*real(s2dtauthetahYs*yhS + s1tauhYs*s1thetahYs');

    d2tautheta = B2dtautheta;
    %% 6.The second partial derivative of frequency offset and theta
    
    s2dxitheta =  kron((dxid1.*b),athetad1);
    s2dxithetahYs = s2dxitheta'*ysVec;
    B2dxitheta = 2*real(s2dxithetahYs*yhS + s1xihYs*s1thetahYs');

    d2xitheta = B2dxitheta;
    %% calculate the Newton's direction and step size
    jacVec = [d1tau; d1xi; d1theta];
    hesMatrix = [d2tau, d2tauxi, d2tautheta; d2tauxi, d2xi, d2xitheta; d2tautheta, d2xitheta, d2theta  ];
    [V,E]=eig(hesMatrix);           % to solve the cases when the Hessian matrix is not negative definite
    e = diag(E);
    V1 = V(:,e<0);
    e1 = e(e<0);
    ss = V1*((V1'*jacVec)./e1);
    %     ss = hesMatrix\jacMatrix;
    step= 1; g = zeros(1,maxIterBLS);
    for kk = 1:maxIterBLS                 % use backtracking line search to find a proper stepsize -- step
        
        x = Value - step*ss;
        athetadamp = exp(thetaRange*sin(x(3)));
        dxidamp = exp(x(2)*xiRange);
  
        zc1damp = exp(1j*pi*r/lTilde*x(1)*(x(1)+zcLen))*exp(-2j*pi*r/lTilde*x(1)*temRange);
        bdamp = [zc1damp conj(zc1damp)].';
        sdamp = kron(dxidamp.*bdamp,athetadamp);        
        yhSdamp = ysVec'*sdamp;
        Bdamp = norm(yhSdamp)^2;
        g(kk) = Bdamp;
        if g(kk) < B - alpha*step*jacVec.'*ss
            step = beta*step;
        else
            break
        end
        if step < 1e-2
            break
        end
    end
    
    Value = Value - step*ss;
    tau = Value(1);
    xi = Value(2);
    theta = Value(3);
    if norm(step*ss) < minStepSize
        break;
    end
    
end
tauEst = Value(1);
xiEst = Value(2);
thetaEst = Value(3);
% rad2deg(thetaEst)
% figure;plot(tmp)
end
