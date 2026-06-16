function [t_vec, H, iter_hist] = solver_implicit(h0, tspan, dt, p, Tol, MaxIter)
%SOLVER_IMPLICIT  Euler Implicito con Newton-Raphson embebido (sistema 2x1).
%
% DESCRIPCION DE LA LOGICA:
%   El metodo de Euler Implicito (backward Euler) define el avance de
%   un paso de tiempo de forma IMPLICITA:
%
%       h_(n+1) = h_n + dt * f(t_(n+1), h_(n+1))
%
%   A diferencia de Euler explicito, h_(n+1) aparece en ambos lados de
%   la ecuacion. Para sistemas no lineales (como f(h) = c*sqrt(h)) no
%   existe forma cerrada, por lo que en cada paso de tiempo se resuelve
%   la ecuacion algebraica no lineal residual:
%
%       G(w) = w - h_n - dt * f(w, p_(n+1)) = 0
%
%   mediante un bucle de Newton-Raphson EMBEBIDO (independiente de
%   newton_multivar.m, que resuelve un problema distinto: el equilibrio
%   global f(h)=0). La Jacobiana analitica de G respecto de w es:
%
%       J_G(w) = I - dt * J_f(w)
%
%   donde J_f es la Jacobiana analitica del campo vectorial
%   (calc_jacobian.m) e I es la identidad 2x2. La actualizacion de
%   Newton es:
%
%       delta_w = -J_G(w)^-1 * G(w)
%       w = w + delta_w
%
%   hasta que norm(delta_w, Inf) < Tol o se agoten MaxIter iteraciones.
%   Como estimacion inicial w^(0) se usa el propio h_n (estrategia
%   robusta: al ser el metodo implicito incondicionalmente estable, la
%   convergencia de Newton no depende criticamente del punto de
%   partida, incluso con pasos de tiempo grandes en escenarios rigidos).
%
%   Se PREASIGNA la matriz H completa antes de iterar en el tiempo. El
%   numero de iteraciones de Newton consumidas en cada paso se registra
%   en iter_hist para poder comparar el COSTO COMPUTACIONAL del metodo
%   implicito frente al explicito (Experimento 3: rigidez/estabilidad).
%
%   ROBUSTEZ FRENTE A LA SATURACION POR DESBORDE: la saturacion definida
%   en campo_vectorial.m introduce, exactamente en h_i = Hmax, un salto
%   (discontinuidad) genuino del campo vectorial: justo debajo de Hmax
%   el tanque sigue llenandose a su tasa natural, y exactamente en Hmax
%   la derivada se fuerza a 0 (pared rigida). Este salto es una
%   propiedad FISICA correcta del fenomeno de rebose (no un error de
%   programacion), pero provoca que la ecuacion residual G(w)=0 no
%   tenga raiz clasica exactamente en ese instante de contacto con la
%   pared, por lo que un Newton puro puede oscilar en un ciclo de
%   periodo 2 alrededor de Hmax sin reducir su paso por debajo de Tol.
%   Para resolverlo de forma robusta, el Newton embebido incorpora dos
%   salvaguardas estandar de "Newton globalizado":
%     a) BACKTRACKING (line search): se acepta un paso completo
%        (alpha=1) solo si reduce la norma del residuo ||G||_inf; en
%        caso contrario se reduce alpha a la mitad sucesivamente. Esto
%        amortigua el ciclo y acerca la iteracion al punto de contacto
%        con una precision de varios ordenes de magnitud mejor que el
%        Newton puro.
%     b) MEJOR ITERADA + ESTANCAMIENTO: se conserva en w_best el
%        iterado con menor ||G||_inf observado hasta el momento. Si
%        durante 3 iteraciones consecutivas no se logra mejorar dicho
%        residuo (estancamiento, tipico de un punto no suave), se
%        detiene el bucle y se reporta w_best en lugar de la ultima
%        iterada, garantizando un resultado preciso y monotono incluso
%        en el instante exacto de saturacion.
%   En el regimen normal (sin saturacion activa), estas salvaguardas no
%   tienen efecto: el paso completo de Newton ya reduce el residuo, por
%   lo que la convergencia cuadratica original se preserva.
%
% ENTRADAS:
%   h0      : vector columna 2x1, condicion inicial [m]
%   tspan   : vector [t0, tf] con limites de integracion [s]
%   dt      : paso de tiempo fijo [s]
%   p       : struct de parametros, o function handle p(t) -> struct
%   Tol     : tolerancia del Newton embebido, norm(delta_w, Inf) < Tol
%   MaxIter : maximo de iteraciones de Newton permitidas por paso
%
% SALIDAS:
%   t_vec     : vector fila 1 x (N+1) con los instantes de tiempo [s]
%   H         : matriz 2 x (N+1) con la trayectoria [h1(t); h2(t)] [m]
%   iter_hist : vector fila 1 x N con el numero de iteraciones de
%               Newton utilizadas en cada paso de tiempo

    N = round((tspan(2) - tspan(1)) / dt);
    t_vec = tspan(1) + (0:N) * dt;

    H = zeros(2, N + 1);      % preasignacion
    H(:, 1) = h0;
    iter_hist = zeros(1, N);

    es_variable = isa(p, 'function_handle');
    I2 = eye(2);
    MAX_ESTANCAMIENTO = 3;     % iteraciones sin mejora antes de aceptar w_best
    ALPHA_MIN = 1/1024;        % piso del backtracking

    for k = 1:N
        h_n = H(:, k);
        t_np1 = t_vec(k + 1);

        if es_variable
            p_np1 = p(t_np1);
        else
            p_np1 = p;
        end

        w = h_n;   % estimacion inicial del Newton embebido

        f0 = campo_vectorial(w, p_np1);
        G0 = w - h_n - dt * f0;
        w_best = w;
        normG_best = norm(G0, Inf);
        estancado = 0;

        for j = 1:MaxIter
            f_w = campo_vectorial(w, p_np1);
            G   = w - h_n - dt * f_w;
            normG = norm(G, Inf);

            if normG < normG_best - 1e-12
                normG_best = normG;
                w_best = w;
                estancado = 0;
            else
                estancado = estancado + 1;
            end

            Jf  = calc_jacobian(w, p_np1);
            JG  = I2 - dt * Jf;
            delta_w = -JG \ G;

            % --- backtracking: solo se acepta el paso si reduce ||G||
            alpha = 1.0;
            while alpha > ALPHA_MIN
                w_trial = w + alpha * delta_w;
                f_trial = campo_vectorial(w_trial, p_np1);
                G_trial = w_trial - h_n - dt * f_trial;
                if norm(G_trial, Inf) < normG
                    break;
                end
                alpha = alpha / 2;
            end
            w = w + alpha * delta_w;

            iter_hist(k) = j;
            paso = norm(alpha * delta_w, Inf);

            if paso < Tol || estancado >= MAX_ESTANCAMIENTO
                break;
            end
        end

        % se reporta siempre la mejor iterada observada (robusta ante kinks)
        f_w = campo_vectorial(w, p_np1);
        G   = w - h_n - dt * f_w;
        if norm(G, Inf) < normG_best
            w_best = w;
        end

        H(:, k + 1) = w_best;
    end
end
