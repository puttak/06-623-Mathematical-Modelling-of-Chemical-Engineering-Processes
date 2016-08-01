%{
Writer: Akshay S Tharval
1st draft: Sept 27, 2015
Last modified: Sept 27, 2015

Subject: Assignment2 Q3
 %}

clear all % Clear stored variables
clc % Clear the screen
close all % Close all previously created plots

% Any matrix A

A = [3 6 9; 1 5 3; 7 5 2]
% Starting vector
x0 = [1 0 0]
x0 = x0';

% Randomly assigning the value of error (Tolerance)
err = 5;

while err > (10^-6)
    % Algorithm for Power Method
    % (i+1)th value of eigen vector
    x1 = A*x0;
    lambda = norm(x1,Inf);
    x1 = x1 / lambda;
    % Finding the value of error
    err = norm(x1) - norm(x0);
    x0 = x1;
end

'Largest Eigen Vector is'
x1 = A*x0;
lambda = norm(x1,Inf);
x1 = x1 / lambda

'Corresponding Largest Eigen Value is'
z = (x1'*(A*x1))/(x1'*x1)
