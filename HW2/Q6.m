%{
Writer: Akshay S Tharval
1st draft: Sept 27, 2015
Last modified: Sept 27, 2015

Subject: Assignment2 Q6
 %}

clear all % Clear stored variables
clc % Clear the screen
close all % Close all previously created plots

%Assuming the initial values of x and y (here x1 and x2)
x1 = 0;
x2 = 0;
x1x2 = [x1;x2];
% Defining a Jacobian matrix
Jn = zeros(2,2);
% Initial value of Iterations:
iter = 1;
% Setting a while loop for max no of iterations:
err = 10;

while err>(10^(-6))
    f1 = x1x2(2)-(x1x2(1)-1)^2;
    f2 = (x1x2(2)+4)^2 - tan(x1x2(1));
    F = [f1;f2];
    % Defining the Jacobian by manual calculations
    Jn(1,1) = -2*(x1x2(1)-1);
    Jn(1,2) = 1.0;
    Jn(2,1) = -(sec(x1x2(1))^2);
    Jn(2,2) = 2*(x1x2(2)+4);
    % Difference value of x and y
    dx = inv(Jn)*F;
    % Calculating the error
    err = F'*F;
    % Next values of x and y
    x1x2 = x1x2 - dx;
    %Displaying the iteration and the value of x and y
    disp(iter)
    disp('The value of x and y for the above mentioned iteration is')
    disp(x1x2)
    %Incrementing the value of iteration number
    iter = iter + 1;
end
disp('No of iterations required:')
disp(iter-1)
disp('Final value of x and y which satisfies the equations above:')
disp(x1x2)

% For plotting
x1 = linspace(0,10);
z = sqrt(tan(x1))-4;
x2 = (x1-1).^2;
figure
plot(x1,x2,'b', x1, z,'r');
xlabel('x')
ylabel('y')

