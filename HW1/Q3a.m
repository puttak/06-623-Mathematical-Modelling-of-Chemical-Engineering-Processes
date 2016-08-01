%{
Writer: Akshay S Tharval
1st draft: Sept 11, 2015
Last modified: Sept 11, 2015

Subject: Assignment Q3a
%}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots

syms T %Defining a variable T

%Setting the initial concentration of A, B and C respectively
Cao=3;
Cbo=0;
Cco=0;

outputVector = [Cao Cbo Cco]; %Defining output vector

T=1:0.5:20; %Setting the range to vary T
L=length(T); %Length of array T


for a = 1:L
    
    Mat = [1+T(a)/20 0 0; -T(a)/20 1+T(a)/10 0; 0 -T(a)/10 1]; % Defining the matrix of equations
    
    finalVector = inv(Mat)*outputVector'; %Solving the system to find the final concentrations
    
    Cc(a) = finalVector (3,1); % Storing the final concentration of C 
    
end

plot(T,Cc,'-r') %Plotting T vs C


