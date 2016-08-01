%{
Writer: Akshay S Tharval
1st draft: Sept 26, 2015
Last modified: Sept 26, 2015

Subject: Assignment2 Q2
 %}

clear all % Clear stored variables
clc % Clear the screen
close all % Close all previously created plots

% The given matrix
A = [3 2 2 1; 2 3 1 2; -1 1 2 0; 2 4 3 5]

%Part 1
% Function eig gives matrix D of eigenvalues and matrix V as columns are
% corresponding right eigenvectors
[V,D] = eig(A);

% Gives the dimensions of matrix A
[m,n] = size(A);

% Saving each norm in variable n
for i = 1:m
    n(i,1) = norm(V(:,i));
end

% Prints the norm of each eigen value after normalizing
'Norm of each eigenvalue is'
n

% Part 2
% Matrix A as symbolic variable
' Matrix A using symbolic function sym'
S = sym(A)

'Exact Solution for the eigenspace using fumction sym'
[v,d] = eig(A)




