%{
Writer: Akshay S Tharval
1st draft: Sept 7, 2015
Last modified: Sept 9, 2015

Subject: Assignment Q1
%}

clear all
close all
clc

x= 1:2:20 % vector from 1 to 20 having interval of 2
a= randi(50,10,1); % Creating a vector of any 10 randoms from 1 to 50
b= a.^2; % Square of every element in vector a
c= sqrt(a); % Square root of every element in vector a
z= [a,b,c] % Keeping vectors a,b, and c as the first, second and third column of matrix z
title ('Question 1');
xlabel('Second column of Matrix');
ylabel('Vector x');
plot (x,b,'g'); % Plotting vector x and second column of matrix z