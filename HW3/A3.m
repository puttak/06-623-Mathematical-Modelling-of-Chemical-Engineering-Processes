
%{
Writer: Akshay S Tharval
1st draft: Oct 15, 2015
Last modified: Oct 15, 2015

Subject: Assignment3 Q1
 %}

function main1()

clear all
clc
close all

disp('====================================================')
disp('Q1a')

%Case a
E0 = 0.1;
Km1 = 0.3;

% Calling the ode45 funtion to solve the system
[T1, Y1] = ode45(@eval1, [0:0.5:150], [1,0]); 

figure();
% Plotting S
plot(T1, Y1(:,1));
hold all
% Plotting P
plot (T1, Y1(:,2));
hold all
ES1 = (E0 * Y1(:,1))/(Km1 + Y1(:,1));
%Plotting ES
plot(T1, ES1(:,1));
hold all
E1 = E0 - ES1;
% Plotting E
plot(T1, E1(:,1));
hold off

% Graph details
xlabel(' Time ');
ylabel(' Concentration in mM');
title('Concentration in mM v/s Time according to case A');
legend('S','P','ES','E');

disp('Please refer graph')

disp('==================================================')
disp('Q1b')

%Case b
E0 = 0.1;       % mM
Km2 = 0.005;    % mM

% Calling the ode45 funtion to solve the system
[T2, Y2] = ode45(@eval2, [0:0.5:150], [1,0]); 

figure();
% Plotting S
plot(T2, Y2(:,1));
hold all
% Plotting P
plot (T2, Y2(:,2));
hold all
ES2 = (E0 * Y2(:,1))/(Km2 + Y2(:,1));       % mM
%Plotting ES
plot(T2, ES2(:,1));
hold all
E2 = E0 - ES2;      % mM
% Plotting E
plot(T2, E2(:,1));
hold off

% Graph details
xlabel(' Time ');
ylabel(' Concentration in mM');
title('Concentration in mM v/s Time according to case B');
legend('S','P','ES','E');
disp('Please refer graph')

% Comparing P(t) for different parameters
figure();
plot(T1,Y1(:,2),T2,Y2(:,2));

xlabel(' Time ');
ylabel(' Concentration of product (P) (mM)');
title('Comparing P(t) for different parameters');
legend('P: part a','P: part b');

end



% Case a
function S = eval1(t, S)
%The given values of constants
kcat = 0.5;     % s^-1
Km = 0.3;       % mM
S0 = 1;         % mM
E0 = S0 * 0.1;  % mM
p = kcat * E0 * ((S(1))/(Km + S(1)));
s = - p;
S = [s; p];
end

% Case b
function R = eval2(t, R)
% The given values of constants
kcat = 0.08;    % s^-1
Km = 0.005;     % mM
S0 = 1;         % mM
E0 = S0 * 0.1;  % mM
p = (kcat * E0 * ((R(1))/(Km + R(1))));
s = - p;
R = [s; p];
end

%{
Here we use equations as
[ES] = (E0 * S)/(Km + S)
[E] = E0 - ES
[P] = (kcat * E0 * (S/(Km + S)))
[S] = - [P]

[S] is -[P] because the amount of P produced will be the amount of S
reacted.
 %}

%%
%{
Writer: Akshay S Tharval
1st draft: Oct 15, 2015
Last modified: Oct 15, 2015

Subject: Assignment3 Q2
 %}

function main2()

clear all
clc
close all

% Calling ode45 function using tau range 0 to 1
tau = linspace(0, 1);
[T, Y] = ode45(@eval,tau, 0.004);

% Graph details
figure();
plot(T, Y,'b.');
xlabel('tau (t/t*)');
ylabel('theta (T/Tf)');
ylim([0,1]); % Limit of y
title('Theta v/s Tau');
grid on
disp('Q2.')
disp('Please refer to graph')
end

% Defininf the function
function dTdtau = eval(tau, theta)
% Given values of constants
sigma = 5.676 * 10^-8;  % W/m^2*K^4
rho = 8933;             % kg/m^3
d = 0.002;              % m
TF = 1200;              % K
tstar = 100;
dTdtau = (tstar *2 * sigma * (1 - theta^4)* TF^4)/(TF *rho * d * (355.2 + (2 * 0.1004 * theta * TF)));
end

%%
%{
Writer: Akshay S Tharval
1st draft: Oct 16, 2015
Last modified: Oct 16, 2015

Subject: Assignment3 Q3
 %}

function main3()

% We have value of x2 as either 1 or -1
clear all
clc
close all

M = [1, 0 ; 0,0];
options = odeset('Mass',M);

tspan = [0,100];

% x2 can be 1 or -1
[T1, X1] = ode15s(@eval1,tspan,[0, 1],options); 
[T2, X2] = ode15s(@eval1,tspan,[0, -1],options);
% Given condition for x2 (x2 = 0.8) 
[T3, X3] = ode15s(@eval1, tspan,[0, 0.8],options);

% Plotting x2 vs x1 with x2 = 1 and x2 = -1
figure();
plot( X1(:,1) , X1(:,2))
hold on
plot(X2(:,1), X2(:,2));
xlabel('x1');
ylabel('x2');
title('x2 v/s x1');
legend('Plot with x2 = 1','Plot with x2 = -1','location','bestoutside');
disp('Q3.')
disp(' Case a has x2 value as 1 and -1, please refer graph1 for representation')
grid on

% With new value of x2
figure();
plot(X3(:,1) , X3(:,2));
xlabel('x1');
ylabel('x2');
title(' Problem 3: x2 v/s x1');
disp('Please refer graph 2 for represeantation when x2 = 0.8')
grid on
end

% Defining the function eval for evaluating system
function f = eval1(t, x)
theta = pi/3;       % Radians
f1 = (-x(1)+cos(theta))+ (x(2) - sin(theta)); 
f2 = (x(1)^2 + x(2)^2 - 1);                               
f = [f1;f2]; 
end

%%
%{
Writer: Akshay S Tharval
1st draft: Oct 16, 2015
Last modified: Oct 16, 2015

Subject: Assignment3 Q4
 %}

function main4()
clear all
clc
close all

% Given initial values of concentrations of every component
C0 = [10 0 0 0];
%Calling ode45 to solve the system
[T, C] = ode45(@eval, [0 25], C0);

disp('Refer Graph 1 and its subplots for individual concentraion profile')
% Graphical details
figure();
subplot(2,2,1) % Indicates that this is subplot 1
plot(T, C(:,1),'r');
ylabel('Concentration of A in Mol');
xlabel('Time in s');
title('Graph 1a: Conc of A vs Time');
legend('Concentration of A')

subplot(2,2,2) % Indicates that this is subplot 2
plot(T, C(:,2),'y');
ylabel('Concentration of B in Mol');
xlabel('Time in s');
title('Graph 1b: Conc of B vs Time');
legend('Concentration of B')

subplot(2,2,3) % Indicates that this is subplot 3
plot(T, C(:,3),'b');
ylabel('Concentration of C in Mol');
xlabel('Time in s');
title('Graph 1c: Conc of C vs Time');
legend('Concentration of C')

subplot(2,2,4) % Indicates that this is subplot 4
plot(T, C(:,4),'g');
ylabel('Concentration of D in Mol');
xlabel('Time in s');
title('Graph 1d: Conc of D vs Time');
legend('Concentration of D')

disp('Refer graph 2 for combined concentration profile of components')
% New graph which includes the concentration profile of every component
figure();
plot(T,C);
ylabel('Concentration of A, B, C and D in Mol');
xlabel('Time in s');
title('Graph 2: Concentration of all components vs Time');
legend('Concentration of A','Concentration of B','Concentration of C','Concentration of D')
end

function dC = eval(t, C)
k = [1; 1; 0.5; 1]; % Given values of rate constants in min^-1
dCa = -k(1) * C(1);                  
dCb = k(1) * C(1) - k(2) * C(2) + k(3) * C(3) - k(4) * C(2);
dCc = k(2) * C(2) - k(3) * C(3);
dCd = k(4)* C(2);
dC = [dCa, dCb, dCc, dCd]';
end

%%
%{
Writer: Akshay S Tharval
1st draft: Oct 16, 2015
Last modified: Oct 16, 2015

Subject: Assignment3 Q5
 %}

function main5()
clear all
clc
close all

% Initial conditions
C0 = [10 0 0 0];
% Calling ode45 for solving the system
[T, C] = ode45(@eval1, [0 30], C0);

% Initial conditions
C1 = [10 0 0 0];
% Calling ode45 for solving the system
[T1, C1] = ode45(@eval2, [0 30], C1);

disp('Q5')
disp('Please refer graph1 and its subplot for concentraions in Batch reactor')
figure();

% For subplot 1
subplot(2,2,1)
plot(T, C(:,1),'r');
xlabel('Concentration of A in Mol');
ylabel('Time in s');
title('Graph 1a: Concentration of A vs Time');
legend('Concentration of A')

% For subplot 2
subplot(2,2,2)
plot(T, C(:,2),'y');
xlabel('Concentration of B in Mol');
ylabel('Time in s');
title('Graph 1b: Concentration of B vs Time');
legend('Concentration of B')

% For subplot 3
subplot(2,2,3)
plot(T, C(:,3),'b');
xlabel('Concentration of C in Mol');
ylabel('Time in s');
title('Graph 1c: Concentration of C vs Time');
legend('Concentration of C')

% For subplot 4
subplot(2,2,4)
plot(T, C(:,4),'g');
xlabel('Concentration of D in Mol');
ylabel('Time in s');
title('Graph 1d: Concentration of D vs Time');
legend('Concentration of D')

disp('Please refer Graph 2 for comparison of concentration profile in CSTR and Batch Reactor')
% Comparison between CSTR and Batch reactor for concentration A
figure();

% Plot of Batch Reactor
plot(T,C(:,1));
hold all
% Plotting of CSTR Reactor
plot(T1,C1(:,1));
legend('CSTR-Concentration of A','BATCH-Concentration of A','location','bestoutside');
hold off
xlabel('Time in min')
ylabel('Concentrations of A, B, C, D in Mol');
title('Graph 2: Concentrations of A, B, C, D in Mol v/s Time in min');
end

% Function to evaluate system in case of CSTR
function dC = eval1(t, C)
k = [1; 1; 0.5; 1];     % GIiven value of rate constants in min^-1            
dCa = 10 -C(1) - k(1) * C(1);
dCb = 0 -C(2) + k(1) * C(1) - k(2) * C(2) + k(3) * C(3) - k(4) * C(2);
dCc = 0 -C(3) + k(2) * C(2) - k(3) * C(3);
dCd = 0 -C(4)+ k(4)* C(2);
dC = [dCa; dCb; dCc; dCd];
end

% Function to evaluate system in case of Batch
function dC = eval2(t, C)
k = [1; 1; 0.5; 1];     % GIiven value of rate constants in min^-1
dCa = -k(1) * C(1);                  
dCb = k(1) * C(1) - k(2) * C(2) + k(3) * C(3) - k(4) * C(2);
dCc = k(2) * C(2) - k(3) * C(3);
dCd = k(4)* C(2);
dC = [dCa; dCb; dCc; dCd];
end

%%
%{
Writer: Akshay S Tharval
1st draft: Oct 16, 2015
Last modified: Oct 16, 2015

Subject: Assignment3 Q6
 %}

function main6()
clear all
clc
close all

% Initial value
C = [150;150];
% Calling ode45 to solve the system
[T,Y] = ode45(@eval, [0 100], C);

% Plotting Temperature of liquid vs time
plot(T, Y(:,1),'b');
% Specifying the limits of x and y axis
axis([0 5 0 100]);
hold on
% Plotting Temperature of container vs time
plot(T,Y(:,2),'y');
axis([0 5 0 100]);
hold off
xlabel('Time in sec');
ylabel('Temperature of Liquid and Container in F');
legend('Liquid','Container')
title('Temperature of Container/Liquid v/s Time');
grid on
end


function DL = eval(t, C)
% Given initial values of constants
rhoL = 62;                       % lbm/ft^3
rhoC = 139;                      % lbm/ft^3
cpL = 1;                         % Btu/lbm*F
cpC = 0.2;                       % Btu/lbm*F
volL = 0.03;                     % ft^3
volC = 0.003;                    % ft^3
AL = 0.4;                        % ft^2
AC = 0.5;                        % ft^2
h = 8.8;                         % Btu/hr*ft^3*F

% Given equation of dL/dt
dL = AL* h *(C(2) - C(1))/(rhoL * cpL * volL);

% Given equation of dC/dt
dC = (AC* h *(32 - C(2))/(rhoC * cpC * volC))+(AC*(C(1) - C(2))/(rhoC * cpC * volC));
% Returning the value when function is called
DL = [dL; dC];
end

%%
%{
Writer: Akshay S Tharval
1st draft: Oct 18, 2015
Last modified: Oct 18, 2015

Subject: Assignment3 Q7
 %}

function main7()

clc
clear all
close all

% Initial values of y at given values of x
y0 = [-1 0 -4];
% x ranges from 0 to 1
x0 = 0:0.01:1;
% Calling ode45 for solving the system
[T,Y] = ode45(@eval,x0,y0);

figure();
plot(T,Y(:,1));
xlabel('X');
ylabel('Y');
title('Plot of x vs y');
grid on

% Initial value of error
error = 5;
% Same range of X
x1 = 0:0.01:1;
a = 0.1;
y1 = [a 0 -4];
% Iterating while the error is less than 0.005
    while error>(5*10^-3)
    [T1,Y1] = ode45(@eval, x1, y1);
    error = abs(Y1(end,1)+1);
    % Updating the initial value of a
    a = a+0.01;
    y1 = [a 0 0-4];
    end
    
disp('The end value of y at x = 0 is ')
disp(Y1(end,1))
% Calling the function ode45 for solving the system
[T, Y2] = ode45(@eval, x1, [a 0 -4]);

figure();
plot(T, Y2);
xlabel('X');
ylabel('Y');
title('X vs Y');
grid on
end

function D = eval(x,y)
dy = 2*x^2 + 2*x + 2*y(1) - y(3);
D = [y(2); y(3); dy];
end

%%
%{
Writer: Akshay S Tharval
1st draft: Oct 18, 2015
Last modified: Oct 18, 2015

Subject: Assignment3 Q8
 %}

function main8()
clc
clear all
close all

% Range of 1 to 2 with 20 intervals in between
xspan = linspace(1,2,20);

% Initialize the vlaue of c with first value (Given)
c0 = [6.308447;10];


% Using the function ode45
[X,Y1] = ode45(@eval,xspan,c0);

% Initialize the value of c with second value (Given)
c1 = [6.308447;15];

% Using the function ode45
[X,Y2] = ode45(@eval,xspan,c1);


% Next value of c
c2 = c0 + ((c1 - c0)*(55.430456 - Y1(end,1))/(Y2(end,1)-Y1(end,1)));

% Evaluate the value of the function at new value of c
[X,Y] = ode45(@eval,xspan,c2);

% Graphical representation
figure();
plot(xspan,Y(:,1));
xlabel('x');
ylabel('y(x)');
title('y(x) v/s x');
grid on

% Displaying the solution at each iteration
disp('solution is');
disp(Y(:,1));
end

% Defining the function for evaluating the value of the given system
function S = eval(x,U)
dy = U(2);
dt =(3 * exp(2*x))- (2*sin(x))+ U(2) - U(1);
S = [dy ; dt];
end

%%
%{
Writer: Akshay S Tharval
1st draft: Oct 18, 2015
Last modified: Oct 18, 2015

Subject: Assignment3 Q9
 %}

%{
Q9a.
-u" + (pi ^2 *u) = 2 * (pi ^ 2) * sin(pi * x)

- Non - Homogenous
    -u" + (pi ^2 *u) = 2 * (pi ^ 2) * sin(pi * x
    Taking u' as y1 and therefore u'' will be y1'
    -y1' + (pi ^2 *u) = 2 * (pi ^ 2) * sin(pi * x)

- Homogenous
    u" + (pi ^2 *u) = 0
    Taking u' as y2 and therefore u'' will be y2'
    -y2' + (pi ^2 *u) = 0
    u(0) = u(1) = 0

Thus we have converted a boundary condition problem into two initial value
problem
 %}


% Q9b
function main9()
clc
clear all
close all

disp('Q9a is in the commments of the program')
disp('=====================================================')
disp('Q9b')
% Non Homogenous
% Given value of x according to its range
xspan = [0:0.25:1];
% Initial values
u0 = [0;0];
u1 = [0;4];

% Calling ode45 to solve the system
[X1,Y1] = ode45(@eval1,xspan,u0);
[X2,Y2] = ode45(@eval1,xspan,u1);

% Calculating the new value of u
u2 = u0 + ((u1 - u0)*(0 - Y1(end,1))/(Y2(end,1)-Y1(end,1)));
% Calling ode45 to solve the system
[X3,Y3] = ode45(@eval1,xspan,u2);

% Displaying the value of solution at each iteration
disp('Value of u1(x) is');
disp(Y3(:,1))

% For Homogenous
% Here we use the same values of X and u
u11 = [0;1];
u12 = [0;4];
xspan = 0:0.25:1;

% Calling ode45 to solve the system
[X4,Y4] = ode45(@eval2,xspan,u11);
[X5,Y5] = ode45(@eval2,xspan,u12);

% Calculating the new value of u
u13 = u11 + ((u12 - u11)*(0 - Y4(end,1))/(Y5(end,1)-Y4(end,1)));
% Calling ode45 to find new value of u
[X6,Y6] = ode45(@eval2,xspan,u13);
% Displaying the solution at each iteration
disp('Value of u2(x) is')
disp(Y6(:,1))

% To check for Dirichlet condition
if (Y3(end,1) == 0 && Y6(end,1)== 0) 
    disp('The prediction of Dirichlet condition is correct @ x =1');
else
    disp('The prediction of Dirichlet condition is incorrect @ x =1');
end

disp('=======================================================')
disp('Q9c')
%Q9c
u = 0;
% For calculaiton of c
c = (u - Y3(end,1))/Y6(end,1);
disp('Value of c is')
disp(c)

% Calculation for w(x)
wx = Y3(:,1) + (c * Y6(:,1));
disp('Value of w(x) as asked is');
disp(wx)

disp('=====================================================')
disp('Q9d')
%Q9d
ua = sin((pi * xspan));
% For calculation of error
error = ua - wx';
disp('Exact error at each point is');
disp(error')
end

% Defining the funtion for evaluation of Non-Homogenous system
function M = eval1(x,u)
du1 = u(2);
dt1 = ((pi)^2*u(1)) - (2*((pi)^2)*sin(pi*x));
M = [du1 dt1]';
end

% Defining the function for evaluation of Homogenous system
function N = eval2(x,u)
du2 = u(2);
dt2 =(((pi)^2)*u(1));
N = [du2 dt2]';
end