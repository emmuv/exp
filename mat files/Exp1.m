% tf=w^2/(s^2+2*zeta*w*s+w^2)
clc; clear all;
% tf=w^2/(s^2+2*zeta*w*s+w^2)
w=1;
t=(0:0.1:20); num=(w^2)
z1=input('enter zeta value')
den=[1 2*z1*w w^2]
y=step(num,den,t);
plot(t,y);
hold
stepinfo(tf(num,den)) 