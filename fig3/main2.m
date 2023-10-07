

clc
clear
close all

global A B F C Gn Gl

%%  inputs
A = [-0.1738 0 36.9771 -6.2589
    0 0 0 1
    1 0 0 0
    -0.0091 0 1.9333 -1.9872];
B = [-1.0095;0;0;-0.3205];
F=[0;0;0;0-0.3205];
C = [0 -1 0 0
    0 0 1 0
    0 0 0 -1];
Gn = [0 9.83 0
    -1 0 0
    0 1 0
    0 0 -1];
Gl = [0 192.49 6.29
    -15 0 -1
    0 25.82 0
    0 1.84 1-5.01];

%%  call ode45 method
[t,X] = ode45(@odeEquations,0:0.01:10,zeros(8,1));

%%  model
y = @(t)(0.6*sin(5*t));

%%  inputs
h = 60/1000;
time = h:h:10;

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

%%  ISOH
ISOH = @(t,y,ts)(y(1)*(t*(t-h))/2 + y(2)*((t-2*h)*(t)) + y(3)*(t*(t+h))/2);
j = 1;
for k = 1:2:length(time)+1
    tisoh(3*j-2:3*j,1) = [k*h;(k+1)*h;(k+2)*h];
    yisoh(3*j-2:3*j,1) = [lagrangePolynomial([(k-2)*h,(k-1)*h,k*h],[y((k-2)*h),y((k-1)*h),y(k*h)],(k)*h)
                          lagrangePolynomial([(k-2)*h,(k-1)*h,k*h],[y((k-2)*h),y((k-1)*h),y(k*h)],(k+1)*h)
                          lagrangePolynomial([(k-2)*h,(k-1)*h,k*h],[y((k-2)*h),y((k-1)*h),y(k*h)],(k+2)*h)];
    j = j + 1;
end

%%  results
figure, hold on
plot(0:0.01:time(end),y(0:0.01:time(end)),'b--')
plot(tzoh,yzoh,'r-')
legend('f(t)','fhat(t)'), title('zoh')
xlabel('t (s)'), ylabel('Amplitude')
xlim([0,10])

figure, hold on
plot(0:0.01:time(end),y(0:0.01:time(end)),'b--')
plot([0;tfoh+h],[0;yfoh],'r-')
legend('f(t)','fhat(t)'), title('foh')
xlabel('t (s)'), ylabel('Amplitude')
xlim([0,10])

figure, hold on
plot(0:0.01:time(end),y(0:0.01:time(end)),'b--')
plot([0;tisoh+h],[0;yisoh],'r-')
legend('f(t)','fhat(t)'), title('isoh')
xlabel('t (s)'), ylabel('Amplitude')
xlim([0,10])