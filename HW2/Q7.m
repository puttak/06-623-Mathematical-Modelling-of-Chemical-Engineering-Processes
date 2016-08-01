%{
Writer: Akshay S Tharval
1st draft: Sept 27, 2015
Last modified: Sept 27, 2015

Subject: Assignment2 Q7
 %}

clear all % Clear stored variables
clc % Clear the screen
close all % Close all previously created plots

%Defining the values of Antoine's Constants for two compounds
A1 = -3848.09;
A2 = -4328.12;
B1 = 17.5318;
B2 = 17.913;

% Boiling points of coumpound Alpha and Beta
Tbalpha = (A1/(log(760)-B1))-273.15 %79.93 degree celcius
Tbbeta = (A2/(log(760)-B2))-273.15 % 110.56 degree celcius

%Range of x
x = 0:0.1:1

%Initiate T
T = ones(1,10)

 % Calculation1 for P_bubble
for i = (1:11)
    p1 = @(T) (x(i).*(exp(A1./(T+273.15) + B1)))+((1-x(i)).*(exp(A2./(T+273.15)+B2)))-760;
    T(i) = fzero(p1,60);
end

%Calculation2 for P_dew
for i = (1:11)
    p2 = @(y) ((y./(exp(A1./(T(i)+273.15)+B1))) + ((1-y)./exp(((A2)./(T(i)+273.15))+ B2)) - (1/760));
    y(i) = fzero(p2,0);
    
end

%Plotting the values of x (from 0 to 1) and y (found from the above
%calculation2) against the different values of T from calculation1
plot(x,T,'r',y,T,'b')
hold on;
 
% Incase x1 = x2
x1 = 0.5;
p1 = @(T) x1*(exp(A1/(T+273.15) + B1))+((1-x1).*(exp(A2/(T+273.15)+B2)))-760;
T1 = fzero(p1,0)

% For the value of y if x1 = x2
y1 = x1*(exp(A1./(T1+273.15)+B1))/760
plot(x1,T1,'o',y1,T1,'co')
grid on
legend('x','y')

