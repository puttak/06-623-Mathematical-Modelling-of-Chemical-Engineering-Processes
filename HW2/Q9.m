%{
Writer: Akshay S Tharval
1st draft: Sept 27, 2015
Last modified: Sept 27, 2015

Subject: Assignment2 Q9
 %}

function Ass2Prob9 ()

clear all % Clear stored variables
clc % Clear the screen
close all % Close all previously created plots

% Calls the function defined after the main function
fun = @eval;
% Initial guess of all q values
q0 = ones(1,7);
% Solve the non-linear equations simultaneously
q = fsolve(fun, q0)
end


function f = eval(q)

% Given values of li
l1 = 100;
l2 = 100;
l3 = 200;
l4 = 75;
l5 = 100;
l6 = 75;
l7 = 50;

% Non-linear function defined from the system
f(1) = q(1) - q(6) - q(2);
f(2) = q(2) - q(4) - q(3);
f(3) = q(3) + q(4) - q(5);
f(4) = q(5) + q(6) - q(7);
f(5) = l3 * q(3)^2 - l4*q(4)^2;
f(6) = l2 * q(2)^2 + l4 * q(4)^2 + l5 * q(5)^2 - l6 * q(6)^2;
f(7) = l1 * q(1)^2 + l6 * q(6)^2 + l7 * q(7)^2 - (5.2*10^5)*(((3.14^2)*(0.2)^5)/8/0.02/998);

end

