function conv_table = convergencia_espacial()
% CONVERGENCIA_ESPACIAL
% Calcula la convergencia espacial (L∞ vs h) para el sistema acoplado CNLS
% usando la función nls() y devuelve la tabla de resultados.

    % --- ACTIVAR PARALELISMO ---
    % Verificar si existe un pool paralelo, si no, crearlo
    if isempty(gcp('nocreate'))
        parpool('local', 6);  % Ajusta el número de workers según tu máquina
        fprintf('✓ Pool paralelo creado\n\n');
    else
        fprintf('✓ Pool paralelo ya activo\n\n');
    end

    % --- Parámetros físicos ---
    alphaA = 1;
    ee     = 1;
    vv     = 1;
    Tfinal = 4;

    % --- Solución exacta (solitón simétrico para sistema acoplado) ---
    % ψ1 = ψ2 = sqrt(2α/(1+e)) * sech(sqrt(2α)*(x-vt)) * exp(i*(vx - (v²/2 - α)t))
    psi_exact = @(x,t) [ ...
        sqrt(2*alphaA/(1+ee)) .* sech(sqrt(2*alphaA)*(x - vv.*t)) .* ...
            exp(1i*(vv*x - (0.5*vv^2 - alphaA)*t)); ...
        sqrt(2*alphaA/(1+ee)) .* sech(sqrt(2*alphaA)*(x - vv.*t)) .* ...
            exp(1i*(vv*x - (0.5*vv^2 - alphaA)*t)) ...
    ];

    % --- Condición inicial ---
    InitialCondition = @(x) psi_exact(x, 0);

    % --- Dominio espacial ---
    xL = -20;
    xR = 60;

    % --- Mallas a probar (según paper) ---
    h_list = [0.2, 0.1, 0.05, 0.025];

    % --- Regla para Δt ---
    c_dt  = 0.1;
    k_of_h = @(h) c_dt * h^2;

    % --- Parámetros de Newton (TOLERANCIA RELAJADA) ---
    tolNewton   = 1e-8;   % Relajado de 1e-10 para converger más rápido
    maxNewtonIt = 20;     % Reducido de 30
    useSparse   = true;
    showNewton  = false;

    % --- Almacenar resultados ---
    Linf_err = zeros(numel(h_list), 1);
    x_cells  = cell(numel(h_list), 1);  % Guardar las mallas x (por si luego las quieres)

    fprintf('Calculando convergencia espacial EN PARALELO...\n');
    fprintf('═══════════════════════════════════════\n');
    fprintf('Parámetros \n');
    fprintf('Tolerancia Newton: %.0e\n', tolNewton);
    fprintf('Ejecutando %d simulaciones en paralelo...\n\n', numel(h_list));

    % --- BUCLE PARALELO ---
    parfor r = 1:numel(h_list)
        h = h_list(r);
        k = k_of_h(h);

        % Número de intervalos espaciales
        Nx = round((xR - xL) / h);

        fprintf('[Worker] Malla h = %.4f, Δt = %.6f, Nx = %d - INICIANDO...\n', h, k, Nx);

        % --- Llamar a la función nls ---
        [U, x, t, stats] = nls(alphaA, ee, vv, Tfinal, InitialCondition, ...
            'Nx', Nx, ...
            'k', k, ...
            'xL', xL, ...
            'xR', xR, ...
            'tolNewton', tolNewton, ...
            'maxNewtonIt', maxNewtonIt, ...
            'useSparse', useSparse, ...
            'showNewton', showNewton, ...
            'saveEvery', round(Tfinal/k));  % Solo guardar estado final

        %#ok<NASGU> % stats y t no se usan explícitamente aquí

        % --- Calcular error L∞ al tiempo final ---
        % U contiene [psi1; psi2] en la última columna
        N = numel(x) - 1;
        psi1_num = U(1:N+1, end);
        psi2_num = U(N+2:end, end);

        % Solución exacta al tiempo final
        psi_ex  = psi_exact(x, Tfinal);
        psi1_ex = psi_ex(1:N+1);
        psi2_ex = psi_ex(N+2:end);

        % Error L∞ máximo entre ambas componentes
        err1 = max(abs(psi1_ex - psi1_num));
        err2 = max(abs(psi2_ex - psi2_num));
        Linf_err(r) = max(err1, err2);
        x_cells{r}  = x;  %#ok<NASGU>

        fprintf('[Worker] Malla h = %.4f - COMPLETADA. Error L∞ = %.6e\n', h, Linf_err(r));
    end

    fprintf('\n✓ Todas las simulaciones completadas\n\n');

    %% --- Calcular orden de convergencia observado ---
    hs   = h_list(:);
    errs = Linf_err(:);

    % Orden entre mallas consecutivas: log(e_i/e_{i+1}) / log(h_i/h_{i+1})
    order = [NaN; diff(log(errs)) ./ diff(log(hs))];

    %% --- Tabla de resultados ---
    conv_table = table(hs, errs, order, ...
        'VariableNames', {'h', 'L_inf_Error', 'Order'});

    fprintf('═══════════════════════════════════════\n');
    fprintf('TABLA DE CONVERGENCIA ESPACIAL\n');
    fprintf('═══════════════════════════════════════\n\n');
    disp(conv_table);

end
