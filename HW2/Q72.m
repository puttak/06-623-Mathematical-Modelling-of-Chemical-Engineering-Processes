clc
clear all

A1 = -3848.09
A2 = -4328.12
B1 = 17.5318
B2 = 17.913
 
Tba = (A1/(log(760)-B1))-273.15
Tbb = (A2/(log(760)-B2))-273.15
 
x = 0:0.1:1
T = ones(1,10)

 
for i = (1:11)
    p1 = @(T) (x(i).*(exp(A1./(T+273.15) + B1)))+((1-x(i)).*(exp(A2./(T+273.15)+B2)))-760;
    T(i) = fzero(p1,60);
end
disp(x)
disp(T)
 
for i = (1:11)
    p2 = @(p) ((p./(exp(A1./(T(i)+273.15)+B1))) + ((1-p)./exp(((A2)./(T(i)+273.15))+ B2)) - (1/760));
    p(i) = fzero(p2,0);
    
end
disp(p)
 
 
plot(x,T,'r',p,T,'b')
hold on;
 
% For equimolar concentrations,
x1 = 0.5;
p1 = @(T) x1*(exp(A1/(T+273.15) + B1))+((1-x1).*(exp(A2/(T+273.15)+B2)))-760;
T1 = fzero(p1,0)
 
plot(x1,T1,'o')
y1 = x1*(exp(A1./(T1+273.15)+B1))/760
 
plot(y1,T1,'co')


