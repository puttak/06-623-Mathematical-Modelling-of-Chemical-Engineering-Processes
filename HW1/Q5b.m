%{
Writer: Akshay S Tharval
1st draft: Sept 11, 2015
Last modified: Sept 11, 2015

Subject: Assignment Q5b
%}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots


n = 15 % Given number of nodes

deltaPL= 101325/3/25400; % The value of DeltaP/L in terms of Pascal/Micron

deltaY= 150/(n+1); % Calculating the value of Delta Y

u= 8.9*10^-4; % Viscosity of water at 25 Celcius in Pascal.Seconds

G = deltaY^2*deltaPL/u; %calculating G, which is basically a contains all the constants

% Considering the boundary condition

%Creating the sparse matrix using sparse function
D = sparse(1:n,1:n,2*ones(1,n),n,n); 
E = sparse(2:n,1:n-1,-1*ones(1,n-1),n,n);
matV = E+D+E'

GMat(1:n)=G; %Creating a vector of G

Vel=inv(matV)*GMat' %Solving the system of equations

Vel = [0;Vel;0] %Because we consider the boundary condition, we need to put a zero row on the top and bottom of the output vector

N1= 0:n+1; %Creating an array of for plotting

plot(Vel,N1,'-g') % Creating velocity profile by plotting Velocity at each node


