%{
Writer: Akshay S Tharval
1st draft: Sept 26, 2015
Last modified: Sept 26, 2015

Subject: Assignment2 Q1
 %}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots

% tau = V/nu
%tau = input('Please enter the value of tau ') 
tau = 10; %Since V = 10L and nu = 1L/sec

% Rate constants
%k = input('Please enter the values of rate constants along with square brackets ')
k = [0.1 0.2 0.8 0.1];
k = k;

%Initial Concentrations
%ci = input('Please enter the initial concentration of A, B, C and D along with the square brackets ')
ci = [5 0 0 1];
ci = ci';

%After modelling we get the matrix as
A = [1+k(1,1)*tau 0 0 0; -k(2,1)*tau 1 k(4,1)*tau k(2,1)*tau-k(3,1)*tau; 0 -k(4,1)*tau 1 0; 0 k(3,1)*tau-k(2,1)*tau 0 1];

% LU Decomposition
[L,U] = lu(A);

%Displays matrix L
disp('After LU decomposition, we get L matrix as')
L

%Displays matrix U
disp('After LU decomposition, we get U matrix as')
U

'LU Decomposition'
%We know Y=(L^-1)*b
'Vector formed due to first multiplication of a triangular matrix by the feed vector '
y = inv(L)*ci
% X= inv(U)*Y

'The Steady state reactor concentrations are '
css = inv(U)*y