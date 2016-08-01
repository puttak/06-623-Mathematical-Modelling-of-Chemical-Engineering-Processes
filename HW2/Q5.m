%{
Writer: Akshay S Tharval
1st draft: Sept 27, 2015
Last modified: Sept 27, 2015

Subject: Assignment2 Q5
 %}

clear all % Clear stored variables
clc % Clear the screen
close all % Close all previously created plots

% In this function to reduce the derivative complexity, we will assume x =
% sqrt(2*g*h)
% For using fzero, we need to use the following format
func = @ (x) x*tanh(0.3*x)-5;
z = fzero(func,5)

%Finding the value of height from the substituted x= sqrt(2*g*h)
ht = (z^2)/2/9.81

%Defining the function and its derivative
func = inline('x*tanh(0.3*x)-5','x');
funci = inline('(1-tan(0.3*x)*tan(0.3*x)*0.3 + tanh(0.3*x))','x');

% Newton's method:
err= 5;
x0 = 5;
func = inline('x*tanh(0.3*x)-5','x');
funci = inline('(1-tan(0.3*x)*tan(0.3*x)*0.3 + tanh(0.3*x))','x');

while err>1*10^(-6)
    x1 = x0 - func(x0)/funci(x0);
    err = abs(x1-x0);
    x0 = x1;
end

disp('The height after Newtons Method')
height = (x0^2)/2/9.81

