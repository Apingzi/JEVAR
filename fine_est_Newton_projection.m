function [ tauEst, xiEst, thetaEst ] = fine_est_Newton_projection(yRefactor,LambdahYs,inter,r,lTilde,zcLen,xi,theta,tau,thetaRange,xiRange,d,lambda,M)
% use the Newton¡¯s iteration method to search for a fine estimate (with backtracking line search)
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
    yRefactorhS = yRefactor'*s;
    LambdaS = s - inter*yRefactorhS;
    yhLambdaS = LambdahYs'*s;
    B = norm(LambdahYs'*s)^2;
    A = norm(LambdaS)^2;
%     tmp(nn) = B/A;
    %% 1.The first/second derivative of delay
    sd1tau1 = (1j*pi*r/lTilde*(2*tau+zcLen-2*temRange)).*zc1;
    sd1tau = [sd1tau1 conj(sd1tau1)].';
    sd2tau1 = (1j*pi*r/lTilde*2 + (1j*pi*r/lTilde*(2*tau+zcLen-2*temRange)).^2 ).*zc1;
    sd2tau = [sd2tau1 conj(sd2tau1)].';
    
    
    s1tau = kron(dxi.*sd1tau,atheta);
    s2tau = kron(dxi.*sd2tau,atheta);
    
    yRefactorhS1 = yRefactor'*s1tau;
    LambdaS1tau = s1tau - inter*yRefactorhS1;
    s1tauhLambdahYs = s1tau'*LambdahYs;
    s2tauhLambdahYs = s2tau'*LambdahYs;
    B1dtau = 2*real(s1tauhLambdahYs*yhLambdaS);
    
    A1dtau = 2*real(s1tau'*LambdaS);
    d1numehalf1 = B1dtau*A-B*A1dtau;
    d1tau = d1numehalf1/A^2;
    B2dtau = 2*abs(LambdahYs'*s1tau)^2 + 2*real(s2tauhLambdahYs*yhLambdaS);
    A2dtau = 2*norm(LambdaS1tau)^2 + 2*real(s2tau'*LambdaS);
    d2tau = ((B2dtau*A-A2dtau*B)*A^2 - d1numehalf1*2*(A'*A1dtau))/A^4;
    
    %% 2.The first/second derivative of frequency offset
    dxid1 = xiRange.*dxi;
    dxid2 = xiRange .* dxid1;
    s1xi = kron((dxid1.*b),atheta);
    s2xi = kron((dxid2.*b),atheta);
    
    yRefactorhS1 = yRefactor'*s1xi;
    LambdaS1xi = s1xi - inter*yRefactorhS1;
    s1xihLambdahYs = s1xi'*LambdahYs;
    s2xihLambdahYs = s2xi'*LambdahYs;
    B1dxi = 2*real(s1xihLambdahYs*yhLambdaS);
    A1dxi = 2*real(s1xi'*LambdaS);
    d1numehalf2 = B1dxi*A-B*A1dxi;
    d1xi = d1numehalf2/A^2;
    B2dxi = 2*abs(LambdahYs'*s1xi)^2 + 2*real(s2xihLambdahYs*yhLambdaS);
    A2dxi = 2*norm(LambdaS1xi)^2+2*real(s2xi'*LambdaS);
    d2xi = ((B2dxi*A-A2dxi*B)*A^2 - d1numehalf2*2*(A'*A1dxi))/A^4;
    %% 3.The first/second derivative of DOA
    athetad1 = thetaRange*cos(theta).*atheta;
    athetad2 = ((1j*2*pi*d*sin(theta)/lambda*(0:M-1)).'-(4*pi^2*d^2*cos(theta)^2/lambda^2*(0:M-1).^2).').*atheta;
    s1theta = kron(xib,athetad1);
    s2theta = kron(xib,athetad2);
    
    yRefactorhS1 = yRefactor'*s1theta;
    LambdaS1theta = s1theta - inter*yRefactorhS1;
    s1thetahLambdahYs = s1theta'*LambdahYs;
    s2thetahLambdahYs = s2theta'*LambdahYs;
    B1dtheta = 2*real(s1thetahLambdahYs*yhLambdaS);
    A1dtheta = 2*real(s1theta'*LambdaS);
    d1numehalf3 = B1dtheta*A-B*A1dtheta;
    d1theta = d1numehalf3/A^2;
    B2dtheta = 2*abs(LambdahYs'*s1theta)^2 + 2*real(s2thetahLambdahYs*yhLambdaS);
    A2dtheta = 2*norm(LambdaS1theta)^2+2*real(s2theta'*LambdaS);
    d2theta = ((B2dtheta*A-A2dtheta*B)*A^2 - d1numehalf3*2*(A'*A1dtheta))/A^4;
    %% 4.The second partial derivative of delay and frequency offset
    s2dtauxi = kron(dxid1.*sd1tau,atheta);
    s2dtauxihLambdahYs = s2dtauxi'*LambdahYs;
    B2dtauxi = 2*real(s2dtauxihLambdahYs*yhLambdaS + s1tauhLambdahYs*s1xihLambdahYs');
    A2dtauxi = 2*real(s2dtauxi'*LambdaS + s1tau'*LambdaS1xi);
    d2tauxi = ((B2dtauxi*A + B1dtau*A1dxi - A2dtauxi*B - A1dtau*B1dxi)*A^2 - d1numehalf1*2*real(A1dxi'*A))/A^4;
    
    %% 5.The second partial derivative of delay and theta
    
    s2dtautheta =   kron(dxi.*sd1tau,athetad1);
    s2dtauthetahLambdahYs = s2dtautheta'*LambdahYs;
    B2dtautheta = 2*real(s2dtauthetahLambdahYs*yhLambdaS + s1tauhLambdahYs*s1thetahLambdahYs');
    A2dtautheta = 2*real(s2dtautheta'*LambdaS + s1tau'*LambdaS1theta);
    d2tautheta = ((B2dtautheta*A + B1dtau*A1dtheta - A2dtautheta*B - A1dtau*B1dtheta)*A^2 - d1numehalf1*2*real(A1dtheta'*A))/A^4;
    %% 6.The second partial derivative of frequency offset and theta
    
    s2dxitheta =  kron((dxid1.*b),athetad1);
    s2dxithetahLambdahYs = s2dxitheta'*LambdahYs;
    B2dxitheta = 2*real(s2dxithetahLambdahYs*yhLambdaS + s1xihLambdahYs*s1thetahLambdahYs');
    A2dxitheta = 2*real(s2dxitheta'*LambdaS + s1xi'*LambdaS1theta);
    d2xitheta = ((B2dxitheta*A + B1dxi*A1dtheta - A2dxitheta*B - A1dxi*B1dtheta)*A^2 - d1numehalf2*2*real(A1dtheta'*A))/A^4;
    %% calculate the Newton's direction and step size
    jacMatrix = [d1tau; d1xi; d1theta];
    hesMatrix = [d2tau, d2tauxi, d2tautheta; d2tauxi, d2xi, d2xitheta; d2tautheta, d2xitheta, d2theta  ];
    [V,E]=eig(hesMatrix);                   % to solve the cases when the Hessian matrix is not negative definite
    e = diag(E);
    V1 = V(:,e<0);
    e1 = e(e<0);
    ss = V1*((V1'*jacMatrix)./e1);
    %     ss = hesMatrix\jacMatrix;
    step= 1; g = zeros(1,maxIterBLS);
    for kk = 1:maxIterBLS                   % use backtracking line search to find a proper stepsize -- step
        
        x = Value - step*ss;
        athetadamp = exp(thetaRange*sin(x(3)));
        dxidamp = exp(x(2)*xiRange);
        zc1damp = exp(1j*pi*r/lTilde*x(1)*(x(1)+zcLen))*exp(-2j*pi*r/lTilde*x(1)*temRange);
        bdamp = [zc1damp conj(zc1damp)].';
        sdamp = kron(dxidamp.*bdamp,athetadamp);       
        Bdamp = norm(LambdahYs'*sdamp)^2;
        yRefactorSdamp  = yRefactor'*sdamp;
        Adamp = norm(sdamp - inter * yRefactorSdamp)^2;
        g(kk) = Bdamp/Adamp;
        if g(kk) < B/A - alpha*step*jacMatrix.'*ss
            step = beta*step;
        else
            break
        end
        if step < 1e-2
            break
        end
    end
%     tmp(nn) = g(kk);
    Value = Value - step*ss;
    
    if norm(step*ss) < minStepSize
        break;
    end
    
    tau = Value(1);
    xi = Value(2);
    theta = Value(3);

end
tauEst = Value(1);
xiEst = Value(2);
thetaEst = Value(3);

end
