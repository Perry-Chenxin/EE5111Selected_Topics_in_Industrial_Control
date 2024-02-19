clc
clear all

B = [0;0;-1];
A = [0 1 0; 0 0 1; 0 -5 -25];
Q = [100 0 0; 0 100 0; 0 0 100];
R = 10;

K = - lqr(A,B,Q,R) % ki kp kd

