clear all;
close all;
clc;

syms a;
P=[a, 0,1;
    0,a,2;
    1,2,a];

lambda=eig(P);
