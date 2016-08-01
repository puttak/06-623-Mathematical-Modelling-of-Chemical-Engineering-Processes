%{
Writer: Akshay S Tharval
1st draft: Sept 27, 2015
Last modified: Sept 27, 2015

Subject: Assignment2 Q8
 %}

function Prob8()
clear all % Clear stored variables
clc % Clear the screen
close all % Close all previously created plots

% Given values of K1 and K2 in case a
K11 = 0.1071;
K21 = 0.01493;

% Given initial concentrations of components
NoCO21 = 0.33333333333;
NoO21 = 0.33333333333;
NoN21 = 0.33333333333;
NoCO1 = 0;
NoNO1 = 0;

% Given values of K1 and K2 in case a
K12 = 0.1071;
K22 = 0.01493;

% Given initial concentrations of components
NoCO22 = 2;
NoO22 = 0.33333333333;
NoN22 = 0.33333333333;
NoCO2 = 0;
NoNO2 = 0;
E = [0.5; 0.5];

%Calling functioneval1 for case a
fun1 = @(E) eval1 (E, K11, K21, NoCO21, NoO21, NoN21, NoCO1, NoNO1);
fun1(E);
%Calling functioneval1 for case b
fun2 = @(E) eval1(E, K12, K22, NoCO22, NoO22, NoN22, NoCO2, NoNO2);
fun2(E);
options = optimoptions('fsolve','Display','iter');
[Ea] = fsolve(fun1, E, options);
[Eb] = fsolve(fun2, E, options);


disp('For case a the extent of reaction is (E1 and E2)')
Ea
disp('For case b the extent of reaction is (E1 and E2)')
Eb
end

function func = eval1(E, K1, K2, NoCO2, NoO2, NoN2, NoCO, NoNO)

%Defining the functions which will be returned when called
func(1) = ((((NoCO + 2*E(1))^2)*(NoO2 - E(2) + E(1)))/(((NoCO2 - 2*E(1))^2)*(NoCO2 + NoCO + NoO2 + NoN2 + NoNO + E(1))))-K1;
func(2) = (((NoNO + 2*E(2))^2)/((NoO2 + E(1) - E(2))*(NoN2 - E(2))))-K2;

end

