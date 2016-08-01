%{
Writer: Akshay S Tharval
1st draft: Sept 11, 2015
Last modified: Sept 11, 2015

Subject: Assignment Q3
%}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots

global Ca Cb Cc 
syms T
Cao=3;
Cbo=0;
Cco=0;


% finalVector = [Ca Cb Cc]

outputVector = [3 0 0]
T=1:0.5:20
L=length(T)

% Case 1

for a = 1:L
    Mat = [1+T(a)/20 0 0; -T(a)/20 1+T(a)/10 0; 0 -T(a)/10 1];
%    j = inv(Mat);
    finalVector = inv(Mat)*outputVector'
    Cc(a) = finalVector (3,1);
    
end

plot(T,Cc)

% Case 2

