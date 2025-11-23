function [U, x, t, stats] = nls(alphaA, ee, vv, Tfinal, InitialCondition, varargin)
% nls  CNLS 1D (Galerkín P1) + punto medio implícito + Newton
%   Sistema acoplado (condiciones de Neumann homogéneas):
%     i ψ1_t + (1/2) ψ1_xx + (|ψ1|^2 + e|ψ2|^2) ψ1 = 0
%     i ψ2_t + (1/2) ψ2_xx + (e|ψ1|^2 + |ψ2|^2) ψ2 = 0
%
% Entradas:
%   alphaA, ee, vv, Tfinal, InitialCondition(x)-> [psi1; psi2]
%
% Opcionales (Name-Value):
%   'Nx',8001, 'k',1e-2, 'xL',-20, 'xR',60, 'tolNewton',1e-10, 'maxNewtonIt',30,
%   'useSparse',true, 'saveEvery',1, 'showNewton',false, 'ExactSolution',[]
%
% Salidas:
%   U: (2*(N+1) x Ns), x: nodos, t: tiempos guardados
%   stats: .err_table, .cons_table (si nls_invariants disponible), .I1_analytic

% --------- Parseo
p = inputParser; p.KeepUnmatched = true;
addParameter(p,'Nx',8001);
addParameter(p,'k',1e-2);
addParameter(p,'xL',-20);
addParameter(p,'xR',60);
addParameter(p,'tolNewton',1e-10);
addParameter(p,'maxNewtonIt',30);
addParameter(p,'useSparse',true);
addParameter(p,'saveEvery',1);
addParameter(p,'showNewton',false);
addParameter(p,'ExactSolution',[]);
parse(p,varargin{:});
Nx=p.Results.Nx; k=p.Results.k; xL=p.Results.xL; xR=p.Results.xR;
tolNewton=p.Results.tolNewton; maxNewtonIt=p.Results.maxNewtonIt;
useSparse=p.Results.useSparse; saveEvery=p.Results.saveEvery;
showNewton=p.Results.showNewton; psi_exact=p.Results.ExactSolution;

% --------- Mallado y matrices FEM P1
N  = Nx;
h  = (xR - xL)/N;
x  = linspace(xL, xR, N+1).';
tgrid = 0:k:Tfinal; Nt = numel(tgrid)-1;

dM = (h/6) * [2; 4*ones(N-1,1); 2];
oM = (h/6) * ones(N,1);
M  = spdiags([ [oM;0], dM, [0;oM] ], [-1,0,1], N+1, N+1);

% K = matriz de rigidez: (φ'_i, φ'_j)  -> corresponde al término + (1/2) ψ_xx en la forma débil
dK = (1/h) * [1; 2*ones(N-1,1); 1];
oK = -(1/h) * ones(N,1);
K  = spdiags([ [oK;0], dK, [0;oK] ], [-1,0,1], N+1, N+1);

if ~useSparse, M = full(M); K = full(K); end

% --------- Bloques reales (u1,u2,u3,u4) = (Re ψ1, Im ψ1, Re ψ2, Im ψ2)
I  = speye(N+1);
A2 = [sparse(N+1,N+1), -I,  sparse(N+1,N+1), sparse(N+1,N+1);   % multiplica por i
      I,                sparse(N+1,N+1),     sparse(N+1,N+1), sparse(N+1,N+1);
      sparse(N+1,N+1),  sparse(N+1,N+1),     sparse(N+1,N+1), -I;
      sparse(N+1,N+1),  sparse(N+1,N+1),     I,               sparse(N+1,N+1)];

MB = blkdiag(M, M, M, M);
KB = blkdiag(K, K, K, K);

% --------- Estado inicial
U0 = InitialCondition(x);                 % [psi1; psi2], 2*(N+1) x 1
psi1_0 = U0(1:N+1);
psi2_0 = U0(N+2:end);
unRI = [real(psi1_0); imag(psi1_0); real(psi2_0); imag(psi2_0)];

takeSnap = @(jj) (mod(jj, saveEvery) == 0);
Ns   = sum(arrayfun(@(jj) takeSnap(jj), 0:Nt));
U    = complex(zeros(2*(N+1), Ns));
tsave = zeros(1, Ns);
sidx = 1;
U(:,sidx) = U0; tsave(sidx) = tgrid(1); sidx = sidx + 1;

norminf = @(v) norm(v, inf);

% --------- Avance temporal
for n = 0:Nt-1
    wRI = unRI;  % inicial para Newton

    for it = 1:maxNewtonIt
        avg = 0.5*(wRI + unRI);

        u1 = avg(1:N+1);              % Re(ψ1)
        u2 = avg(N+2:2*(N+1));        % Im(ψ1)
        u3 = avg(2*(N+1)+1:3*(N+1));  % Re(ψ2)
        u4 = avg(3*(N+1)+1:end);      % Im(ψ2)

        z1  = u1.^2 + u2.^2;          % |ψ1|^2
        z2  = u3.^2 + u4.^2;          % |ψ2|^2
        z1t = z1 + ee*z2;             % coef ψ1
        z2t = ee*z1 + z2;             % coef ψ2

        % ---- Residuo: i M (w-u)/k + (1/2) K * avg + M * f(avg) = 0
        lin_mid = -0.5 * (KB * avg);

        nonlin  = MB * [ z1t .* u1;
                         z1t .* u2;
                         z2t .* u3;
                         z2t .* u4 ];

        F = A2 * (MB*((wRI - unRI)/k)) + lin_mid + nonlin;

        if showNewton
            fprintf('n=%d it=%d ||F||_inf=%.3e\n', n, it, norminf(F));
        end
        if norminf(F) <= tolNewton, break; end

        % ---- Jacobiano de f(avg)
        dz1_u1 = 2*u1;        dz1_u2 = 2*u2;        dz1_u3 = 2*ee*u3;    dz1_u4 = 2*ee*u4;
        dz2_u1 = 2*ee*u1;     dz2_u2 = 2*ee*u2;     dz2_u3 = 2*u3;       dz2_u4 = 2*u4;

        % comp1: z1t*u1
        D11 = spdiags(dz1_u1.*u1 + z1t, 0, N+1, N+1);
        D12 = spdiags(dz1_u2.*u1,      0, N+1, N+1);
        D13 = spdiags(dz1_u3.*u1,      0, N+1, N+1);
        D14 = spdiags(dz1_u4.*u1,      0, N+1, N+1);
        % comp2: z1t*u2
        D21 = spdiags(dz1_u1.*u2,      0, N+1, N+1);
        D22 = spdiags(dz1_u2.*u2 + z1t,0, N+1, N+1);
        D23 = spdiags(dz1_u3.*u2,      0, N+1, N+1);
        D24 = spdiags(dz1_u4.*u2,      0, N+1, N+1);
        % comp3: z2t*u3
        D31 = spdiags(dz2_u1.*u3,      0, N+1, N+1);
        D32 = spdiags(dz2_u2.*u3,      0, N+1, N+1);
        D33 = spdiags(dz2_u3.*u3 + z2t,0, N+1, N+1);
        D34 = spdiags(dz2_u4.*u3,      0, N+1, N+1);
        % comp4: z2t*u4
        D41 = spdiags(dz2_u1.*u4,      0, N+1, N+1);
        D42 = spdiags(dz2_u2.*u4,      0, N+1, N+1);
        D43 = spdiags(dz2_u3.*u4,      0, N+1, N+1);
        D44 = spdiags(dz2_u4.*u4 + z2t,0, N+1, N+1);

        DB = [D11 D12 D13 D14;
              D21 D22 D23 D24;
              D31 D32 D33 D34;
              D41 D42 D43 D44];

        % Jacobiano total de F
        J = A2*(MB/k) - 0.5*KB + MB*DB;

        % Paso de Newton
        dw  = J \ (-F);
        wRI = wRI + dw;
    end

    unRI = wRI;

    if takeSnap(n+1)
        wC = [complex(unRI(1:N+1), unRI(N+2:2*(N+1)));                   % ψ1
              complex(unRI(2*(N+1)+1:3*(N+1)), unRI(3*(N+1)+1:end))];    % ψ2
        U(:,sidx)   = wC;
        tsave(sidx) = tgrid(n+2);
        sidx = sidx + 1;
    end
end

t = tsave;

% --------- Error L∞ si hay solución exacta (magnitud)
if ~isempty(psi_exact)
    errLinf  = zeros(numel(t),1);
    for j = 1:numel(t)
        psi_ex_j = psi_exact(x, t(j));  % [psi1_exact; psi2_exact]
        err1 = abs( abs(psi_ex_j(1:N+1)) - abs(U(1:N+1,j)) );
        err2 = abs( abs(psi_ex_j(N+2:end)) - abs(U(N+2:end,j)) );
        errLinf(j) = max([err1; err2]);
    end
    stats.err_table = table(t(:), errLinf(:), 'VariableNames', {'t','Linf_err'});
else
    stats.err_table = table(t(:), nan(numel(t),1), 'VariableNames', {'t','Linf_err'});
end

% --------- Cantidades conservadas (si la función existe)
try
    cons = nls_invariants(x, U, t, 1, 'e', ee);
    stats.cons_table = cons.table;
catch
    stats.cons_table = table();
end

% --------- Masa analítica típica para solitón simétrico (opcional)
stats.I1_analytic = 2*sqrt(2*alphaA) / (1 + ee);
end