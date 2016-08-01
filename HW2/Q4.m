%{
Writer: Akshay S Tharval
1st draft: Sept 27, 2015
Last modified: Sept 27, 2015

Subject: Assignment2 Q4
 %}

clear all % Clear stored variables
clc % Clear the screen
close all % Close all previously created plots

disp('Analytical method')
disp('Roots using analytical method')
fx = [1 -5 7 -3];
s = roots(fx)

% Newton Raphson Method

disp('Using the standard update function')

% Initial guess (Given)
x = 0;

%Initial value of iteration
iter = 1;

for iter = (1:5)
    % Original functions
    f = x^3 - 5*x^2 + 7*x -3;
    f1 = 3*x^2 - 10*x + 7;
    f2 = 6*x - 10;
    % Standard Newton Raphson
    u = -(f/f1);
    % Updating the value of x
    x = x + u;
    % Iteration number
    disp(iter)
    % Value of x for the corresponding iteration
    disp('The value of x on above iteration is ')
    disp(x)
    iter = iter + 1;
end

disp('Using the modified update function')

% Initial guess (Given)
y = 0;

%Initial value of iteration
iter = 1;

for iter = (1:5)
    % Original function in terms of y
    f = y^3 - 5*y^2 + 7*y -3;
    f1 = 3*y^2 - 10*y + 7;
    f2 = 6*y - 10;
    % Modified update function of Newton Raphson is given by
    u = -f*f1/(f1^2-f*f2);
    %Updating the value of y
    y = y + u;
    %Iteration number
    disp(iter)
    %Value of y for the corresponding iteration
    disp('The value of y on above iteration is ')
    disp(y)
    iter = iter + 1;
end

%Comparing the two update functions

disp('Using the standard update function')

% Initial guess (Given)
x = 4;

%Initial value of iteration
iter = 0;

err = 5;

while err > (10^-6)
    % Original functions
    f = x^3 - 5*x^2 + 7*x -3;
    f1 = 3*x^2 - 10*x + 7;
    f2 = 6*x - 10;
    % Standard Newton Raphson
    u = -(f/f1);
    % Updating the value of x
    x1 = x + u;
    err = abs(x-x1);
    x = x1;
    iter = iter + 1;
end

disp('Number of iterations needed')
disp(iter)
disp('Value of x after above mentioned iterations is')
disp(x)


disp('Using the modified update function')

% Initial guess (Given)
y = 4;

%Initial value of iteration
iter1 = 0;

err1 = 5;

while err1 > (10^-6)
    % Original function in terms of y
    f = y^3 - 5*y^2 + 7*y -3;
    f1 = 3*y^2 - 10*y + 7;
    f2 = 6*y - 10;
    % Modified update function of Newton Raphson is given by
    u = -f*f1/(f1^2-f*f2);
    %Updating the value of y
    y1 = y + u;
    err1 = abs(y - y1);
    y = y1;
    iter1 = iter1 + 1;
end

disp('Number of iterations needed')
disp(iter1)
disp('Value of x after above mentioned iterations is')
disp(x)
