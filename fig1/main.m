

clc
clear
close all

%%  run simulink model
open('fig1_sim')
sim('fig1_sim')

%%  results
figure, hold on
plot(tout,f_t_1,'b-',tout,f_t_11,'c+-')
plot(tout,f_t_2,'k-',tout,f_t_22,'m-*')
plot(tout,f_t_3,'r-',tout,f_t_33,'g-^')
legend('f(t - 1)','fhat(t - 1.1)','f(t - 2)','fhat(t - 2.2)','f(t - 3)','fhat(t - 3.3)')
xlim([5,20])
xlabel('t (s)'), grid on