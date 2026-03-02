% Example input skeleton
% ---- Geometry ----
L = 1.0;              % choose L = 1 for simplicity
EA = 1000;            % choose EA = 1000 (any consistent value works)

A = 1.0;              % pick A = 1
E = EA/A;             % so EA stays 1000

nodes = [ 0   0;      % Node 1
          L   0;      % Node 2
          0   L ];    % Node 3

% ---- Connectivity ----
elems = [1 2;         % bottom member
         1 3;         % vertical member
         3 2];        % diagonal member

% ---- Loads ----
loads = zeros(6,1);

F3 = 10;              % horizontal load at node 2
F4 = 5;               % vertical load at node 2

loads(3) = F3;        % DOF 3 = node 2 x
loads(4) = F4;        % DOF 4 = node 2 y

% ---- Boundary Conditions ----
% Node 1 fixed (DOF 1,2)
% Node 3 fixed (DOF 5,6)
fixedDofs = [1 2 5 6];

% ---- Solve ----
res = truss2d_solve(nodes, elems, E, A, loads, fixedDofs);

% ---- Display Results ----
disp('Displacements (u):')
disp(res.u)

disp('Reactions (R):')
disp(res.R)

disp('Member Axial Forces (N):')
disp(res.N)