%{
Writer: Akshay S Tharval
1st draft: Sept 11, 2015
Last modified: Sept 11, 2015

Subject: Assignment Q3b
%}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots

%Case 1
N = 4
deltaPL= 101325/3/25400; % The value of DeltaP/L in terms of Pascal/Micron
deltaY= 150/(N+1); % Calculating the value of Delta Y
u= 8.9*10^-4; % Viscosity of water at 25 Celcius in Pascal.Seconds
G = deltaY^2*deltaPL/u; %calculating G, which is basically a contains all the constants

%Since N=4, there will be 4 equations
% Considering the boundary condition

matV= [2 -1 0 0;-1 2 -1 0;0 -1 2 -1; 0 0 -1 2];
GMat=[G; G; G; G;];
Vel=inv(matV)*GMat
N1= 1:4;
plot(Vel,N1,'-g')
