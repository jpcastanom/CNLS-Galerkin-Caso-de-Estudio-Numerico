function cons = nls_invariants(x, U, t, varargin)
% nls_invariants  Invariantes para CNLS según el paper de Ismail (2008).
% Calcula I1, I2, I3, I4 usando trapecio nodal:
%   I1 = ∫ |ψ₁|² dx          (masa componente 1)
%   I2 = ∫ |ψ₂|² dx          (masa componente 2)
%   I3 = ∫ i·Σⱼ(ψ̄ⱼ·∂ψⱼ/∂x - ∂ψ̄ⱼ/∂x·ψⱼ) dx   (momento)
%   I4 = ∫ {(1/2)·Σⱼ|∂ψⱼ/∂x|² - (1/2)·Σⱼ|ψⱼ|⁴ - e·|ψ₁|²·|ψ₂|²} dx  (energía)
%
% Uso:
%  cons = nls_invariants(x, U, t);
%  cons = nls_invariants(x, U, t, 'PaperStep',10, 'e',2/3, 'alpha',[1.0]);
%  cons = nls_invariants(x, U, t, 'e',2/3, 'alpha',[1.0, 0.6, 0.3]);  % 3 solitones

% -------- Parseo de parámetros opcionales
PaperStep  = 10;         % intervalo para tabla impresa
PaperTimes = [];
e_cpl      = 0;          % coeficiente de acoplamiento
alpha_vec  = [];         % vector de alphas [alpha1, alpha2, ...] para múltiples solitones
PrintTable = true;       % imprimir tabla en consola

if ~isempty(varargin)
    for kk = 1:2:numel(varargin)
        name = lower(varargin{kk});
        val  = varargin{kk+1};
        switch name
            case 'paperstep',   PaperStep  = val;
            case 'papertimes',  PaperTimes = val;
            case 'e',           e_cpl      = val;
            case 'ecpl',        e_cpl      = val;
            case 'alpha',       alpha_vec  = val(:);  % forzar vector columna
            case 'printtable',  PrintTable = val;
        end
    end
end

% -------- Datos de malla
N  = numel(x) - 1;
ns = size(U, 2);

% -------- Reservas para almacenar resultados
I1 = zeros(ns, 1);      % masa ψ₁
I2 = zeros(ns, 1);      % masa ψ₂
I3 = zeros(ns, 1);      % momento
I4 = zeros(ns, 1);      % energía

% -------- Cálculo para cada snapshot
for j = 1:ns
    Uj = U(:, j);
    
    % Separar componentes
    if numel(Uj) == 2*(N+1)
        psi1 = Uj(1:N+1);
        psi2 = Uj(N+2:end);
        two_comp = true;
    elseif numel(Uj) == (N+1)
        psi1 = Uj;
        psi2 = [];
        two_comp = false;
    else
        error('nls_invariants: tamaño de U incompatible con x.');
    end
    
    % ---- Derivadas espaciales (diferencias centradas en malla uniforme)
    h = x(2) - x(1);  % espaciado uniforme
    dpsi1 = zeros(N+1, 1);
    dpsi1(1)   = (psi1(2) - psi1(1)) / h;
    dpsi1(N+1) = (psi1(N+1) - psi1(N)) / h;
    for m = 2:N
        dpsi1(m) = (psi1(m+1) - psi1(m-1)) / (2*h);
    end
    
    if two_comp
        dpsi2 = zeros(N+1, 1);
        dpsi2(1)   = (psi2(2) - psi2(1)) / h;
        dpsi2(N+1) = (psi2(N+1) - psi2(N)) / h;
        for m = 2:N
            dpsi2(m) = (psi2(m+1) - psi2(m-1)) / (2*h);
        end
    else
        dpsi2 = zeros(N+1, 1);
    end
    
    % ============================================================
    % CANTIDADES CONSERVADAS 
    % ============================================================
    
    % I1: Masa de ψ₁ (ecuación 7)
    I1(j) = trapz(x, abs(psi1).^2);
    
    % I2: Masa de ψ₂ (ecuación 7)
    if two_comp
        I2(j) = trapz(x, abs(psi2).^2);
    else
        I2(j) = 0;
    end
    
    % I3: Momento (ecuación 8)
    % I3 = ∫ (i)·Σⱼ(ψ̄ⱼ·∂ψⱼ/∂x - ∂ψ̄ⱼ/∂x·ψⱼ) dx
    term1 = 1i*(conj(psi1).*dpsi1 - conj(dpsi1).*psi1);
    
    if two_comp
        term2 = 1i*(conj(psi2).*dpsi2 - conj(dpsi2).*psi2);
        integrand_I3 = term1 + term2;
    else
        integrand_I3 = term1;
    end
    
    I3(j) = real(trapz(x, integrand_I3));
    
    % I4: Energía (ecuación 9)
    % I4 = ∫ {(1/2)·Σⱼ|∂ψⱼ/∂x|² - (1/2)·Σⱼ|ψⱼ|⁴ - e·|ψ₁|²·|ψ₂|²} dx
    if two_comp
        Ekin = 0.5 * (abs(dpsi1).^2 + abs(dpsi2).^2);
        Epot = 0.5 * (abs(psi1).^4 + abs(psi2).^4) + ...
               e_cpl * (abs(psi1).^2 .* abs(psi2).^2);
    else
        Ekin = 0.5 * abs(dpsi1).^2;
        Epot = 0.5 * abs(psi1).^4;
    end
    
    integrand_I4 = Ekin - Epot;
    I4(j) = trapz(x, integrand_I4);
end

t = t(:);

% I3_1   = real(trapz(x, term1));
% I3_2   = real(trapz(x, term2));
% I3_tot = real(trapz(x, term1 + term2));
% 
% fprintf('I3_1   = %.10f\n', I3_1);
% fprintf('I3_2   = %.10f\n', I3_2);
% fprintf('I3_tot = %.10f\n', I3_tot);

% -------- Errores respecto al valor inicial
dI1 = I1 - I1(1);
dI2 = I2 - I2(1);
dI3 = I3 - I3(1);
dI4 = I4 - I4(1);

% Errores relativos
rel_err_I1 = abs(dI1) ./ max(abs(I1(1)), eps);
rel_err_I2 = abs(dI2) ./ max(abs(I2(1)), eps);
rel_err_I3 = abs(dI3) ./ max(abs(I3(1)), eps);
rel_err_I4 = abs(dI4) ./ max(abs(I4(1)), eps);

% Máximos errores relativos
max_rel = struct(...
    'I1', max(rel_err_I1), ...
    'I2', max(rel_err_I2), ...
    'I3', max(rel_err_I3), ...
    'I4', max(rel_err_I4));

% -------- Valor analítico de I1 (si se proporcionó alpha)
% Para N solitones: I1_analytic = (2/(1+e)) * Σⱼ sqrt(2*αⱼ)
% Referencia: paper ecuación (34) para 1 solitón, después de (51) para N solitones
if ~isempty(alpha_vec) && e_cpl > 0
    % Calcular suma de las contribuciones de cada solitón
    I1_analytic = (2/(1+e_cpl)) * sum(sqrt(2*alpha_vec));
    err_I1_analytic = abs(I1 - I1_analytic);
    rel_err_I1_analytic = err_I1_analytic / abs(I1_analytic);
else
    I1_analytic = NaN;
    err_I1_analytic = nan(size(I1));
    rel_err_I1_analytic = nan(size(I1));
end

% -------- Tabla completa
cons.table = table(...
    t(:), I1, I2, I3, I4, ...
    dI1, dI2, dI3, dI4, ...
    'VariableNames', {'t','I1','I2','I3','I4','dI1','dI2','dI3','dI4'});

% -------- Selección de tiempos para imprimir (paper table)
if ~isempty(PaperTimes)
    tpoints = PaperTimes(:).';
else
    t0 = t(1);
    tf = t(end);
    if abs(t0) < 10*eps, t0 = 0; end
    
    if abs(t0) < 10*eps
        tpoints = 0:PaperStep:tf;
    else
        kf = floor((tf - t0)/PaperStep);
        tpoints = t0 + (0:kf)*PaperStep;
        if abs(tf - tpoints(end)) > 1e-12
            tpoints = [tpoints, tf];
        end
    end
end

idx = arrayfun(@(tau) find(abs(t - tau) == min(abs(t - tau)), 1, 'first'), tpoints);
idx = unique(idx, 'stable');

% Crear tabla estilo paper
cons.paper_table = table(...
    t(idx), I1(idx), I2(idx), I3(idx), I4(idx), ...
    'VariableNames', {'Time','I1','I2','I3','I4'});

% Si hay valor analítico, crear tabla de errores
if ~isnan(I1_analytic)
    cons.error_table = table(...
        t(idx), I1(idx), repmat(I1_analytic, numel(idx), 1), ...
        err_I1_analytic(idx), rel_err_I1_analytic(idx), ...
        'VariableNames', {'Time','I1_numerical','I1_analytic','Error_I1','Rel_Error_I1'});
end

% -------- Guardar datos en estructura
cons.I1 = I1;
cons.I2 = I2;
cons.I3 = I3;
cons.I4 = I4;

cons.dI1 = dI1;
cons.dI2 = dI2;
cons.dI3 = dI3;
cons.dI4 = dI4;

cons.rel_err_I1 = rel_err_I1;
cons.rel_err_I2 = rel_err_I2;
cons.rel_err_I3 = rel_err_I3;
cons.rel_err_I4 = rel_err_I4;

cons.max_rel_error = max_rel;

cons.I1_analytic = I1_analytic;
cons.err_I1_analytic = err_I1_analytic;
cons.rel_err_I1_analytic = rel_err_I1_analytic;

% -------- Imprimir tablas en consola
if PrintTable
    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('  CANTIDADES CONSERVADAS\n');
    fprintf('========================================================\n');
    fprintf('  Coeficiente de acoplamiento e = %.4f\n', e_cpl);
    if ~isnan(I1_analytic)
        fprintf('  Valor analítico I1 = %.6f', I1_analytic);
        if numel(alpha_vec) > 1
            fprintf('  (suma de %d solitones)\n', numel(alpha_vec));
            fprintf('  Alphas: [');
            fprintf('%.2f ', alpha_vec);
            fprintf(']\n');
        else
            fprintf('  (solitón simple, alpha = %.2f)\n', alpha_vec);
        end
    end
    fprintf('--------------------------------------------------------\n');
    fprintf('  Time      I1          I2          I3          I4\n');
    fprintf('--------------------------------------------------------\n');
    for kk = 1:numel(idx)
        fprintf('%7.2f  %10.6f  %10.6f  %10.6f  %10.6f\n', ...
            t(idx(kk)), I1(idx(kk)), I2(idx(kk)), I3(idx(kk)), I4(idx(kk)));
    end
    fprintf('========================================================\n');
    
    % Imprimir errores máximos
    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('  ERRORES RELATIVOS MÁXIMOS\n');
    fprintf('========================================================\n');
    fprintf('  I1: %.3e\n', max_rel.I1);
    fprintf('  I2: %.3e\n', max_rel.I2);
    fprintf('  I3: %.3e\n', max_rel.I3);
    fprintf('  I4: %.3e\n', max_rel.I4);
    fprintf('========================================================\n\n');
    
    % Si hay valor analítico, imprimir tabla de errores
    if ~isnan(I1_analytic)
        fprintf('\n');
        fprintf('========================================================\n');
        fprintf('  ERROR EN I1 (respecto al valor analítico)\n');
        fprintf('========================================================\n');
        fprintf('  Time      I1_num      I1_exact    Error       Rel_Error\n');
        fprintf('--------------------------------------------------------\n');
        for kk = 1:numel(idx)
            fprintf('%7.2f  %10.6f  %10.6f  %.7e  %.6e\n', ...
                t(idx(kk)), I1(idx(kk)), I1_analytic, ...
                err_I1_analytic(idx(kk)), rel_err_I1_analytic(idx(kk)));
        end
        fprintf('========================================================\n\n');
    end
end

end