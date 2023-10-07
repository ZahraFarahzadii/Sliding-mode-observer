

clc
clear
close all

%%  model
y = @(t)(sin(2*pi*t));

%%  inputs
h = 0.1;
time = h:h:1;

%%  ZOH
tzoh(1:2,1) = [0;0.1];
for k = 1:length(time)
    tzoh(2*k+1:2*k+2,1) = [k*h;(k+1)*h];
    yzoh(2*k+1:2*k+2,1) = [y(k*h);y(k*h)];
end

%%  FOH
tfoh(1:2,1) = [0;0.1];
for k = 1:length(time)
    tfoh(2*k+1:2*k+2,1) = [k*h;(k+1)*h];
    yfoh(2*k+1:2*k+2,1) = y(k*h) + [(y(k*h) - y((k-1)*h))*(tfoh(2*k+1)-k*h)/h;(y(k*h) - y((k-1)*h))*(tfoh(2*k+2)-k*h)/h];
end

%%  IFOH
for k = 1:length(time)
    tifoh(2*k-1:2*k,1) = [k*h;(k+1)*h];
    yifoh(2*k-1:2*k,1) = y((k-1)*h) + [(y(k*h) - y((k-1)*h))*(tifoh(2*k-1)-k*h)/h;(y(k*h) - y((k-1)*h))*(tifoh(2*k)-k*h)/h];
end

%%  ISOH
j = 1;
for k = 1:2:length(time)+1
    tisoh(3*j-2:3*j,1) = [k*h;(k+1)*h;(k+2)*h];
    yisoh(3*j-2:3*j,1) = [lagrangePolynomial([(k-2)*h,(k-1)*h,k*h],[y((k-2)*h),y((k-1)*h),y(k*h)],(k-2)*h)
                          lagrangePolynomial([(k-2)*h,(k-1)*h,k*h],[y((k-2)*h),y((k-1)*h),y(k*h)],(k-1)*h)
                          lagrangePolynomial([(k-2)*h,(k-1)*h,k*h],[y((k-2)*h),y((k-1)*h),y(k*h)],(k)*h)];
    j = j + 1;
end

%%  results
figure, hold on
plot(0:0.01:time(end),y(0:0.01:time(end)),'g-')
plot(tzoh,yzoh,'b--')
plot(tfoh,yfoh,'r-.')
plot(tifoh-h,yifoh,'m-')
plot(tisoh-2*h,yisoh,'k--')
xlabel('t (s)'), ylabel('Amplitude')
xlim([0,1])