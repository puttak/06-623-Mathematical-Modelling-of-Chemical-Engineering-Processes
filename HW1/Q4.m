%{
Writer: Akshay S Tharval
1st draft: Sept 12, 2015
Last modified: Sept 13, 2015

Subject: Assignment Q4
%}

clear all %clear stored variables
clc %clear the screen
close all %close all previously created plots



%A= input('Please enter system Matrix A (equations) along with square brackets: ')
%b= input('Please enter system Matrix b (output) along with square brackets: ')

A = [3 -1 3 1; 6 0 9 -2; -12 0 -10 5; 72 -8 48 -19]
b = [6; 13; 17; 93]

%Testing if the value of determinent of A is zero
if det(A) == 0
    'The entered system is ill-conditioned, please enter a valid system of matrix'
end

%To test the consistency of the matrix
s1 = size (A);
s2 = size (b);
n = s1(1,1); %saving the value of size of A in n (very important)
if s1(1,1) ~= s2(1,1)
    'The system is not consistent, please re-enter the matrix'
end

%Generating and printing the Augmented Matrix
'Augmented matrix'
C=[A b]

index = 0; % The position of the max number in a column, initially zero
count=0; % Count of number of times pivoting was done

for c=1:n-1 %columns going from 1 to n
    
    for h=c:n-1 %Rows going from 'c' to 'n-1' (because once we are done with the 1st row and 1st column, we should start with 2nd row and column and go upto second last row and column
        %finding the maximum number and its position in the column
        maxNum = max(abs(C(h:n,c)));
        index = find(abs(C(1:n,c)) == maxNum);

        % Pivoting the row of maximum value with the first row
        C([h index],:)= C([index h],:);
        count=count+1;
        
        %Now the pivoting is done and the first row is having its first value as highest in that column we can do the elimination
        
        %Performing elimination
        for i= h:(n-1)
            C(i+1,:) = C(i+1,:)- (C(h,:).*(C(i+1,c)/C(h,c)));
            
        end
        c=c+1; %Incrementing in order to go to the next column
        
    end
end

'The matrix after Gauss Elimination is'
C %Printing the augmented matrix after elimination

%The matrix is converted into required form and now we can perform back
%substitution in order to find the value of variables

x= zeros (n,1); %Creating an array for outputs

%Back Substitution
D = C(:,1:n); %Separating the augmented matrix, first nxn elements in a matrix 
F = C(:,n+1); % Seaparating the last column in other matrix
x(n) = F(n)/D(n,n); %for the solution of last variable (This is done because the calculation required for the last variable is different form other variables

%for calculating the value of other variables
for o= n-1:-1:1
    x(o)= (F(o)-D (o,o+1:n)*x(o+1:n))/D(o,o);
end

%Printing the number of times we performed pivoting
'Number of times we performed pivoting'
count

%calculating the error between calculated 'b' and actual 'b'
error= A*x-b

%Printing the output vector
'Output vector is'
x

%Printing the error
'Error in calculation is'
error

