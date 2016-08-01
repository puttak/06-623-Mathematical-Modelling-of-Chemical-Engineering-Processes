%{
Writer: Akshay S Tharval
1st draft: Sept 11, 2015
Last modified: Sept 11, 2015

Subject: Assignment Q3b
%}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots

global k2
T= 0.5; % The value of T is given

% Assuming Cao=3 mol/L
Cao=3; % Initial concentration of A
Cbo=0; % Initial concentration of B
Cco=0; % Initial concentration of C
k1=0.5; % Value of k1 is given

outputVector = [Cao Cbo Cco];
k2=1:0.5:20; % Range of k2 
L= length(k2); % Length of array of k2

for a = 1:L
   
    Mat = [1+k1*T 0 0; -k1*T 1+T*k2(a) 0; 0 -T*k2(a) 1]; %Defining the matrix 
    
    finalVector = inv(Mat)*outputVector'; % Solving the system for the value of final concentrations of A,B and C
    
    Cc(a) = finalVector (3,1); % Storing the value of Cc
    
end

plot (k2/k1,Cc,'-r') %Plotting k2/k1 vs Cc