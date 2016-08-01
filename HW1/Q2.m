%{
Writer: Akshay S Tharval
1st draft: Sept 9, 2015
Last modified: Sept 11, 2015

Subject: Assignment Q2
%}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots

syms a ; % Creating variable 'a'

A= [-3 -2 1;2 a 1;4 1 -2]; %converting the equations in matrix form along with variable 'a'
b= [5 -2 7]; %Output vector

% for the system to be ill-conditioned, the value of matrix m has to be 0

detA= det(A); %calculating the determinant of A

% The determinant value of A will be an equation in terms of a

s=solve(detA,a); %solving the determinant value (which is in term of a)

% s is the value at which the system is ill-conditioned
'The value at which the system is ill-conditioned is:'
s %Printing the value of 'a' where the system is ill-conditioned

k= 1:0.1:2*s; %Make an array of values from 1 to twice the value of alpha for plotting
L= length(k); %Finding the length of array L

for c=1:L
    M= [-3 -2 1;2 k(c) 1;4 1 -2];
    Condi (c) = cond(M); %Storing the value of condition of matrix M in 'Condi'
   
end
    
% Specifying the line color, marker-size,marker color)
plot(k,Condi,'--r') % Plotting condition value vs value of alpha 

xlabel('Condition of system') %Label of x-axis

ylabel('Value of range of Alpha') % Label of y-axis

grid on % To display axis grid lines


