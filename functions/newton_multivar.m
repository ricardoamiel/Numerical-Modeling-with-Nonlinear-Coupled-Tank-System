function [h_star, n_iter, historial, convergio] = newton_multivar(h0, p, Tol, MaxIter)
%NEWTON_MULTIVAR  Resuelve F(h) = 0 (estado estacionario) por Newton-Raphson.
%
% DESCRIPCION DE LA LOGICA:
%   El estado estacionario (h1*, h2*) anula simultaneamente ambas
%   derivadas, es decir, es la raiz del sistema algebraico no lineal
%   F(h) = f(h) = 0, donde f es exactamente el campo vectorial fisico
%   (campo_vectorial.m). Se aplica la iteracion clasica de
%   Newton-Raphson multivariable:
%
%       h^(k+1) = h^(k) - J(h^(k))^-1 * F(h^(k))
%
%   En cada iteracion:
%     1) Se evalua F = campo_vectorial(h, p_inf), usando una copia de
%        los parametros con Hmax = Inf para garantizar que se busca el
%        equilibrio "libre" del modelo, sin que la saturacion por
%        desborde interfiera con la busqueda de la raiz.
%     2) Se evalua la Jacobiana analitica J = calc_jacobian(h, p).
%     3) Se resuelve el sistema lineal J*delta_h = -F (en lugar de
%        invertir J explicitamente, por estabilidad numerica).
%     4) Se actualiza h = h + delta_h.
%     5) Criterio de paro: norm(delta_h, Inf) < Tol, o se alcanza
%        MaxIter sin convergencia (en cuyo caso convergio = false).
%
%   Dado que F es la misma funcion de las EDOs, la convergencia
%   cuadratica de Newton en un entorno del equilibrio permite obtener
%   h* con alta precision en pocas iteraciones, sirviendo como
%   referencia exacta contra la cual comparar la trayectoria temporal
%   integrada por RK4 / Euler Implicito.
%
% ENTRADAS:
%   h0      : vector columna 2x1, estimacion inicial [m]
%   p       : struct con A1, A2, c1, c2, Qin (Hmax es ignorado aqui)
%   Tol     : tolerancia para norm(delta_h, Inf)
%   MaxIter : numero maximo de iteraciones permitidas
%
% SALIDAS:
%   h_star    : vector columna 2x1, estado estacionario aproximado [m]
%   n_iter    : numero de iteraciones efectivamente realizadas
%   historial : matriz 2 x n_iter con la trayectoria de h^(k)
%   convergio : booleano, true si se alcanzo Tol antes de MaxIter

    p_inf = p;
    p_inf.Hmax = Inf;   % se busca el equilibrio libre, sin saturacion

    h = h0;
    historial = zeros(2, MaxIter);
    convergio = false;
    n_iter = 0;

    for k = 1:MaxIter
        F = campo_vectorial(h, p_inf);
        J = calc_jacobian(h, p);

        delta_h = -J \ F;
        h = h + delta_h;

        n_iter = k;
        historial(:, k) = h;

        if norm(delta_h, Inf) < Tol
            convergio = true;
            break;
        end
    end

    historial = historial(:, 1:n_iter);
    h_star = h;
end
