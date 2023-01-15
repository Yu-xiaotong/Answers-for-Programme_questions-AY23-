clear;
clc;
close all;

%%The rate constants
k1 = 100;      
k2 = 300;
k3 = 150;

%Parameter setting
h = 1e-5;               %step size
t = 0:h:1;              %the vector for the argument t

N = length(t);
a = ones(1,N);          %the concentration of E
b = ones(1,N);          %the concentration of S
c = ones(1,N);          %the concentration of ES
b(1,1) = 10;
c(1,1) = 0;

r = zeros(N,4);         %the matrix to store the rate of change 

%%The fourth-order Runge-Kutta iteration
for i=2:N
    t_n=t(i-1);
    a_n=a(i-1);
    b_n=b(i-1);
    c_n=c(i-1);

    ka1=-k1*a_n*b_n+k2*c_n+k3*c_n;
    kb1=-k1*a_n*b_n+k2*c_n;
    kc1=k1*a_n*b_n-k2*c_n-k3*c_n;
    kd1=k3*c_n;

    ka2=(-k1*(a_n+ka1*h/2)*(b_n+kb1*h/2))+(k2*(c_n+kc1*h/2))+(k3*(c_n+kc1*h/2));
    kb2=(-k1*(a_n+ka1*h/2)*(b_n+kb1*h/2))+(k2*(c_n+kc1*h/2));
    kc2=(k1*(a_n+ka1*h/2)*(b_n+kb1*h/2))+(-k2*(c_n+kc1*h/2))+(-k3*(c_n+kc1*h/2));
    kd2=k3*(c_n+kc1*h/2);

    ka3=(-k1*(a_n+ka2*h/2)*(b_n+kb2*h/2))+(k2*(c_n+kc2*h/2))+(k3*(c_n+kc2*h/2));
    kb3=(-k1*(a_n+ka2*h/2)*(b_n+kb2*h/2))+(k2*(c_n+kc2*h/2));
    kc3=(k1*(a_n+ka2*h/2)*(b_n+kb2*h/2))+(-k2*(c_n+kc2*h/2))+(-k3*(c_n+kc2*h/2));
    kd3=k3*(c_n+kc2*h/2);
    
    ka4=(-k1*(a_n+ka3*h)*(b_n+kb3*h))+(k2*(c_n+kc3*h))+(k3*(c_n+kc3*h));
    kb4=(-k1*(a_n+ka3*h)*(b_n+kb3*h))+(k2*(c_n+kc3*h));
    kc4=(k1*(a_n+ka3*h)*(b_n+kb3*h))+(-k2*(c_n+kc3*h))+(-k3*(c_n+kc3*h));
    kd4=k3*(c_n+kc3*h);
    
    a(i)=a_n+h/6*(ka1+2*ka2+2*ka3+ka4);
    b(i)=b_n+h/6*(kb1+2*kb2+2*kb3+kb4);
    c(i)=c_n+h/6*(kc1+2*kc2+2*kc3+kc4);  
    r(i-1,:)=[ka1,kb1,kc1,kd1];
end


%%Plot
% figure
% hold on;
% subplot(311);
% plot(t,a,'r');
% xlabel('t');
% ylabel('[E]');
% subplot(312);
% plot(t,b,'g');
% xlabel('t');
% ylabel('[S]');
% subplot(313);
% plot(t,c,'b');
% xlabel('t');
% ylabel('[ES]');
% hold off;

figure
hold on;
subplot(411);
plot(t(1:65000),r(1:65000,1),'r');
xlabel('t');
ylabel('r1');
title('The plot of the rate of change of E');
subplot(412);
plot(t(1:65000),r(1:65000,2),'g');
xlabel('t');
ylabel('r2');
title('The plot of the rate of change of S');
subplot(413);
plot(t(1:65000),r(1:65000,3),'b');
xlabel('t');
ylabel('r3');
title('The plot of the rate of change of ES');
subplot(414);
plot(t(1:65000),r(1:65000,4),'y');
xlabel('t');
ylabel('r4');
title('The plot of the rate of change of P');
hold off;


Vm = max(r(:,4));           %Find the greatest Vm
index = find(r(:,4)==Vm);
CS = b(index);

% Plot V as a function of the concentration of S
figure
plot(b,r(:,4),'r');
xlabel('The concentration of S');
ylabel('V');
text(CS,Vm,'(8.9764,99.9124)');
title('The plot of V as a function of [S]');






