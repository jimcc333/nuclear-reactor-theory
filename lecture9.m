%% Lecture 9 Script %%
% Author: Cem Bagdatlioglu
% Date: 2017-09-15
% Determines the neutron flux distribution in a 1D geometry using two
% methods and compares the outputs

%% Matrix Solution
mesh = 31;                      % Number of mech points
dx = 3/(mesh-1);
x = linspace(0,3,mesh);         % populate x with points spanning [0,3]
S = ones(1,mesh);               % uniform volume source
A = zeros(mesh,mesh);
A(1,1) = 1;  S(1) = 0;          % boundary condition f0=0
A(mesh,mesh) = 1;  S(mesh) = 0; % boundary condition fN=0
% populate the inner diagonal of A
for i = 2:mesh-1
    A(i,i-1) = -1/(dx)^2;       % coefficient of fi-1
    A(i,i+1) = -1/(dx)^2;       % coefficient of fi+1
    A(i,i)   = 2/(dx)^2 + 2;    % coefficient of fi
end
phi=linsolve(A,S');             % MATLAB calculates f=A-1S.
plot(x, phi, '+')               % Plot numerical solution with plus signs
hold on;

%% Exact (Theoretical) Solution
exact_phi = 1 - cosh(x-1.5)/cosh(1.5);  % Find the symmetric and exact solution
plot(x,exact_phi)                       % Plot the exact solution
legend('numerical solution', 'exact solution')

