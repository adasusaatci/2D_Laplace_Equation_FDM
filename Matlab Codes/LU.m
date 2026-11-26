%% Cleaning
close all;
clear all;
clc;

%% Parameter declaration
r = 81; 
theta = 41;
MN = r-1;
NM = theta-1;
r_in = 1;
r_out = 2;
theta_min = 0;
theta_max = 2*pi;
dr=(r_out - r_in)/(r-1); %Grid
dt=(theta_max - theta_min)/(theta-1); % Grid
R=r_in:dr:r_out;
T=theta_min:dt:theta_max; 
[r,theta]=meshgrid(r,theta);
NC = r*theta;
A=sparse(NC);   
RHS=sparse(NC,1); 

%% Building A Matrix
% Loop to set up matrix A and RHS
for i=2:MN   
    for j=2:NM    
        k=i+(j-1)*r;    
         % Coefficients for interior nodes
        A(k,k)=(-2/dr^2) - (2/((R(i)^2)*dt^2));       
        A(k,k-1)=(1/dr^2) - (1/(R(i)*2*dr));          
        A(k,k+1)=(1/dr^2) + (1/(R(i)*2*dr));                    
        A(k,k-r)=1/((R(i)^2)*dt^2);            
        A(k,k+r)=1/((R(i)^2)*dt^2); 
        % Right-hand side values
        RHS(k,1)=1;                  
    end
end
for i = 1:r
    for j = 1:theta
        k = i + (j - 1) * r;
        % Setting boundary conditions
        if i == 1 || i == r || j == 1 || j == theta
            A(k, k) = 1;
        end
    end
end
%% LU Factorization
[L,U,P] = lu(A);
x_0 = U \ ( L \ (P*RHS) );

for i=1:r
    for j=1:theta
        k=i+(j-1)*r;
        Uij(i,j)=x_0(k);
    end
end