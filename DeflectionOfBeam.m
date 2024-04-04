% Parameters
L = 10000; % Length of the bridge (mm)
W = 4000; % Load per unit length (N/mm)
E = 200e3; % Young's modulus of A36 steel (MPa)
I = 1.8225e11; % Moment of inertia (mm^4)
FOS = 1.5; % Factor of safety
N = 5; % Number of nodes
h = L/(N-1); % Distance between two consecutive nodes
K = (h^2)/(E*I); % Constant factor

x = linspace(0,L,N); % Discretize the beam

% Bending Moment Values at each nodes
M = zeros(1, N);
for i = 1:N
    M(i) = ((W*L*(i-1)*h)/2) - ((W*((i-1)*h)^2)/2);
end

% Coefficient Matrix
A = zeros(N, N);

for j = 1:N
    A(j,j) = -2;
end

for j = 1:N-1
    A(j,j+1) = 1;
    A(j+1,j) = 1;
end

% Boundary conditions
A(1,:) = 0;
A(1,1) = 1;
A(N,:) = 0;
A(N,N) = 1;

% Constant right-hand side matrix
B = zeros(N, 1);
for i = 1:N
    B(i) = K * M(i);
end
B(1) = 0; % Deflection at the left support is zero
B(N) = 0; % Deflection at the right support is zero

% Deflection matrix
D = A \ B;

% Plot deflection curve
plot(x, D, '-o');
xlabel('Length (mm)');
ylabel('Deflection (mm)');
title('Deflection of Simply Supported Beam');
grid on;
