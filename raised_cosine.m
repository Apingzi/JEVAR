function [ response ] = raised_cosine( upSampRate,rollOffFactor,taud,winLen)
% Wikipedia£ºraised-cosine filter £¬Ts=1
%os_factor=oversampling factor

if nargin <=3 
    winLen = 2;
end
a = rollOffFactor;%/beta
t = (-winLen-taud):1/upSampRate:(winLen-taud); %×ó¼ÓÓÒ¼õ,ÑÓÊ±;¼õtau

p = zeros(length(t),1);
for i=1:1:length(t)
    if t(i)==1/(2*a) || t(i)==-1/(2*a)
        p(i) = (pi/4)*sinc(1/(2*a));
    else
        p(i) = sinc(t(i)).*cos(pi*a*(t(i)))./(1-(2*a*(t(i))).^2);
    end    
end
response = p;%/norm(p); %Normalization to unit energy
% figure;plot(t,response)
end
