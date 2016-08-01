%{
Writer: Akshay S Tharval
1st draft: Sept 11, 2015
Last modified: Sept 11, 2015

Subject: Assignment Q6
%}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots

global Ea lnko
R= 8.314; %Given value of R
k= [0.0026 0.0118 0.046 0.0873 0.18]; %From the given tbale of values of 'K'
T = [430 450 470 480 490]; %From the given values of T from table
L = length(T); %Length of the array T (same as k)

%Calculating the values of (RT)^-1 for every given value
for i=1:L
A(i) = [(T(i)*R)^-1];
end

%Forming the system of matrix (A)
A1= [1 A(1); 1 A(2);1 A(3); 1 A(4); 1 A(5)];

%Calculating the value of lnK for every given value of K and making a
%vector to store the values
for i=1:L
b(i) = [log(k(i))];
end

%Calculating by pseudo-inverse method
x= pinv(A1)*b'

%Calculating by normal solving method
y=A1'/b


KoCal1 = exp(x(1,1)) %by pseudo-inverse
KoCal2 = exp(y(1,1)) %by normal method

%for plotting 
kCal = zeros(1,5);
kCal1 = zeros(1,5);

% Storing the value of calculated 'K' from pseudo-inverse method
for i=1:5
    kCal(i) = KoCal1*exp(-x(2,1)/R/T(i));
    i=i+1
   
end

%Storing the value of calculated 'K' from normal method
for i=1:5
    kCal1(i) = KoCal2*exp(-y(1,1)/R/T(i));
    i=i+1
   
end


figure %initialize a figure window
plot(k,kCal,'--gs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5]) %Plotting Calculated value of 'K' vs the given values of 'K'
xlabel ('Given value of K') %Label on x-axis
ylabel ('Calculated value of K') %Label on y-axis
title ('Using pseudo-inverse') %Title of the figure

figure %initialize a figure window
plot(k,kCal1,'--gs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])%Plotting Calculated value of 'K' vs the given values of 'K'
xlabel ('Given value of K') %Label on x-axis
ylabel ('Calculated value of K')%Label on y-axis
title ('Using normal method to solve a system of equations') %Title of the figure




% After comparison, it is found that the error is very high and thats
% caused due to the use of pseudo-inverse function.
%This function selects the value ehich has the least square and in this
%case that caused the calculated value of 'K' to rise substantially.
%In comparison,normal method of solving yields a better result.