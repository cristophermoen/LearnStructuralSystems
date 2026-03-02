function results = truss2d_solve(nodes, elems, E, A, loads, fixedDofs)
%TRUSS2D_SOLVE 2D truss solver (direct stiffness method).
%
%   results = truss2d_solve(nodes, elems, E, A, loads, fixedDofs)
%
% INPUTS
%   nodes      : (nNode x 2) [x y] coordinates
%   elems      : (nElem x 2) [n1 n2] connectivity (1-based node indices)
%   E          : (nElem x 1) Young's modulus per element OR scalar
%   A          : (nElem x 1) Area per element OR scalar
%   loads      : (2*nNode x 1) global nodal load vector [Fx1; Fy1; Fx2; Fy2; ...]
%   fixedDofs  : vector of constrained DOF indices (1..2*nNode),
%                e.g., node 1 fixed in x,y -> [1 2]
%
% OUTPUT (struct)
%   results.u        : (2*nNode x 1) displacement vector
%   results.R        : (2*nNode x 1) reaction vector (nonzero at fixed DOFs)
%   results.K        : (2*nNode x 2*nNode) assembled global stiffness
%   results.N        : (nElem x 1) member axial forces (tension +)
%   results.delta    : (nElem x 1) member axial extensions (local axis)
%   results.geom     : struct with L,c,s per element
%
% NOTES
%   - Units must be consistent.
%   - Will error if structure is a mechanism / insufficient supports.

    % ---- Basic checks ----
    nNode = size(nodes, 1);
    nElem = size(elems, 1);
    nDof  = 2 * nNode;

    if numel(loads) ~= nDof
        error('loads must have length 2*nNode.');
    end

    if isscalar(E), E = repmat(E, nElem, 1); end
    if isscalar(A), A = repmat(A, nElem, 1); end
    if numel(E) ~= nElem || numel(A) ~= nElem
        error('E and A must be scalar or length nElem.');
    end

    fixedDofs = unique(fixedDofs(:));
    if any(fixedDofs < 1) || any(fixedDofs > nDof)
        error('fixedDofs must be between 1 and 2*nNode.');
    end

    % ---- Assemble global stiffness ----
    K = zeros(nDof, nDof);

    Ls = zeros(nElem,1);
    cs = zeros(nElem,1);
    ss = zeros(nElem,1);

    for e = 1:nElem
        n1 = elems(e,1);
        n2 = elems(e,2);

        x1 = nodes(n1,1); y1 = nodes(n1,2);
        x2 = nodes(n2,1); y2 = nodes(n2,2);

        dx = x2 - x1;
        dy = y2 - y1;
        L  = hypot(dx, dy);
        if L <= 0
            error('Element %d has zero length.', e);
        end
        c = dx / L;
        s = dy / L;

        Ls(e) = L; cs(e) = c; ss(e) = s;

        ke = (E(e)*A(e)/L) * ...
            [ c*c,  c*s, -c*c, -c*s;
              c*s,  s*s, -c*s, -s*s;
             -c*c, -c*s,  c*c,  c*s;
             -c*s, -s*s,  c*s,  s*s ];

        dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2]; % [uix uiy ujx ujy]

        % scatter-add
        K(dofs, dofs) = K(dofs, dofs) + ke;
    end

    % ---- Apply boundary conditions (partition) ----
    allDofs = (1:nDof).';
    freeDofs = setdiff(allDofs, fixedDofs);

    Kff = K(freeDofs, freeDofs);
    Ff  = loads(freeDofs);

    % singularity check (rank deficiency => mechanism/insufficient BCs)
    if rank(Kff) < size(Kff,1)
        error('Kff is singular. Likely insufficient supports or a mechanism.');
    end

    % ---- Solve ----
    u = zeros(nDof,1);
    u(freeDofs) = Kff \ Ff;   % fixed DOFs assumed u=0

    % ---- Reactions ----
    R = K*u - loads;

    % ---- Member axial forces (tension +) ----
    delta = zeros(nElem,1);
    N = zeros(nElem,1);

    for e = 1:nElem
        n1 = elems(e,1);
        n2 = elems(e,2);

        c = cs(e); s = ss(e); L = Ls(e);

        dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
        ue = u(dofs); % [uix uiy ujx ujy]

        % axial extension along member local axis:
        % delta = [-c -s c s] * ue
        d = (-c)*ue(1) + (-s)*ue(2) + (c)*ue(3) + (s)*ue(4);
        delta(e) = d;

        % axial force N = EA/L * delta
        N(e) = (E(e)*A(e)/L) * d;
    end

    % ---- Pack results ----
    results.u = u;
    results.R = R;
    results.K = K;
    results.delta = delta;
    results.N = N;
    results.geom.L = Ls;
    results.geom.c = cs;
    results.geom.s = ss;
end
