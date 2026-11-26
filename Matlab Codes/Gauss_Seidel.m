%% Cleaning
close all;
clear all;
clc;

%% Parameter decleration
r_in = 1;
r_out = 2;
theta_min = 0;
theta_max = 2*pi;
r = 81;
theta = 41;
MN = r-1;
NM = theta-1;
dr=(r_out - r_in)/(r-1); 
dt=(theta_max - theta_min)/(theta-1); 
R=r_in:dr:r_out;
T=theta_min:dt:theta_max; 
[r,theta]=meshgrid(r,theta);
node_count = r*theta;
A=sparse(node_count);   
RHS=sparse(node_count,1); 

%% Building A Matrix
% Loop to set up matrix A and RHS
for i=2:MN    
    for j=2:NM   
        k=i+(j-1)*r;    
        A(k,k)=(-2/dr^2) - (2/((R(i)^2)*dt^2));       %  Uij
        A(k,k-1)=(1/dr^2) - (1/(R(i)*2*dr));          %  U_(i-1,j)
        A(k,k+1)=(1/dr^2) + (1/(R(i)*2*dr));          %  U_(i+1,j)           
        A(k,k-r)=1/((R(i)^2)*dt^2);             % U_(i,j-1)
        A(k,k+r)=1/((R(i)^2)*dt^2);             % U_(i,j+1)     
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

% Drichlet Boundary Conditions
for j=1:theta
        k=r+(j-1)*r;   
        A(k,k)=1;
        RHS(k,1)=(2-1/2)*sin(T(j));
end
%% Building Sub - Matrices
diagonal = diag(diag(A));
lower = tril(-A,-1);
upper = triu(-A,1);

%% Gauss - Siedel Algorithm

err = 1e-6;
itermax = 400000;
it_count= 1;
x_0 = sparse(node_count,1);


while it_count <= itermax
	
	x_1 = x_0;
	x_0 = ((diagonal - lower) \ upper) * x_1 + ((diagonal - lower) \ RHS);
	jacobi_error = norm( x_0 - x_1 );
	if jacobi_error< err
		fprintf('Converged at %d iterations', it_count)
		break
	end
	it_count = it_count + 1;
    convergeplot(1,it_count)=jacobi_error;
end

for i=1:r
    for j=1:theta
        m=i+(j-1)*r;
        Uij(i,j)=x_0(m);
    end
end