function conv_table = convergencia_temporal()
% CONVERGENCIA_TEMPORAL
% Calcula la convergencia temporal (L∞ vs Δt) para el sistema acoplado CNLS
% usando la función nls() y devuelve ÚNICAMENTE la tabla de resultados:
%   columnas: Delta_t, L_inf_Error, Order

    % --- ACTIVAR PARALELISMO ---
    if isempty(gcp('nocreate'))
        parpool('local', 6);  % Ajusta el número de workers si quieres
        fprintf('✓ Pool paralelo creado con 6 workers\n\n');
    else
        fprintf('✓ Pool paralelo ya activo\n\n');
    end

    % --- Parámetros físicos ---
    alphaA = 1;
    ee     = 1;
    vv     = 1;
    Tfinal = 4;

    % --- Solución exacta (solitón simétrico para sistema acoplado) ---
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

    % --- Espaciado espacial fijo (muy fino para no contaminar error temporal) ---
    h_fixed = 0.001;

    % --- Pasos de tiempo a probar ---
    k_list = [0.2, 0.1, 0.05, 0.025];

    % --- Parámetros de Newton ---
    tolNewton   = 1e-12;
    maxNewtonIt = 30;
    useSparse   = true;
    showNewton  = false;

    % --- Almacenar resultados ---
    Linf_err = zeros(numel(k_list), 1);

    % --- Mensajes iniciales (mismo estilo que convergencia_espacial) ---
    fprintf('Calculando convergencia temporal EN PARALELO...\n');
    fprintf('═══════════════════════════════════════\n');
    fprintf('Parámetros \n');
    fprintf('Tolerancia Newton: %.0e\n', tolNewton);
    fprintf('Espaciado espacial fijo: h = %.4f\n', h_fixed);
    fprintf('Ejecutando %d simulaciones en paralelo...\n\n', numel(k_list));

    % --- BUCLE PARALELO ---
    parfor r = 1:numel(k_list)
        k = k_list(r);

        % Número de intervalos espaciales
        Nx = round((xR - xL) / h_fixed);

        fprintf('[Worker] Paso de tiempo Δt = %.4f, h = %.4f, Nx = %d - INICIANDO...\n', ...
                k, h_fixed, Nx);

        % Resolver sistema acoplado CNLS
        [U, x, t, stats] = nls(alphaA, ee, vv, Tfinal, InitialCondition, ...
            'Nx', Nx, ...
            'k', k, ...
            'xL', xL, ...
            'xR', xR, ...
            'tolNewton', tolNewton, ...
            'maxNewtonIt', maxNewtonIt, ...
            'useSparse', useSparse, ...
            'showNewton', showNewton, ...
            'saveEvery', round(Tfinal/k)); %#ok<NASGU>

        % --- Error L∞ al tiempo final ---
        N = numel(x) - 1;
        psi1_num = U(1:N+1, end);
        psi2_num = U(N+2:end, end);

        psi_ex  = psi_exact(x, Tfinal);
        psi1_ex = psi_ex(1:N+1);
        psi2_ex = psi_ex(N+2:end);

        err1 = max(abs(psi1_ex - psi1_num));
        err2 = max(abs(psi2_ex - psi2_num));
        Linf_err(r) = max(err1, err2);

        fprintf('[Worker] Paso de tiempo Δt = %.4f - COMPLETADA. Error L∞ = %.6e\n', ...
                k, Linf_err(r));
    end

    fprintf('\n✓ Todas las simulaciones completadas\n\n');

    % --- Orden de convergencia observado ---
    ks   = k_list(:);
    errs = Linf_err(:);
    order = [NaN; diff(log(errs)) ./ diff(log(ks))];

    % --- Tabla de resultados ---
    conv_table = table(ks, errs, order, ...
        'VariableNames', {'Delta_t', 'L_inf_Error', 'Order'});

    fprintf('═══════════════════════════════════════\n');
    fprintf('TABLA DE CONVERGENCIA TEMPORAL\n');
    fprintf('═══════════════════════════════════════\n\n');
    disp(conv_table);
end
