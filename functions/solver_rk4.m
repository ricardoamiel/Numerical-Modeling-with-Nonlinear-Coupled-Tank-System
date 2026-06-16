function [t_vec, H] = solver_rk4(h0, tspan, dt, p)
%SOLVER_RK4  Integrador explicito de Runge-Kutta de 4to orden (sistema 2x1).
%
% DESCRIPCION DE LA LOGICA:
%   Se discretiza el intervalo [tspan(1), tspan(2)] con paso fijo dt,
%   generando N+1 nodos temporales igualmente espaciados. La matriz de
%   resultados H (2 x N+1) se PREASIGNA con zeros() antes del bucle
%   principal para evitar el crecimiento dinamico de arreglos dentro
%   del ciclo de integracion (clave para la eficiencia computacional).
%
%   En cada paso se evalua el campo vectorial 4 veces (k1..k4) sobre el
%   estado actual y estados intermedios, combinandolos con los pesos
%   clasicos de RK4:
%
%       k1 = f(t_n,        h_n)
%       k2 = f(t_n+dt/2,   h_n + dt/2*k1)
%       k3 = f(t_n+dt/2,   h_n + dt/2*k2)
%       k4 = f(t_n+dt,     h_n + dt*k3)
%       h_(n+1) = h_n + dt/6*(k1 + 2*k2 + 2*k3 + k4)
%
%   SOPORTE DE PARAMETROS VARIABLES EN EL TIEMPO: el argumento p puede
%   ser (a) un struct constante, o (b) un function handle p(t) que
%   retorna el struct de parametros evaluado en el instante t. Esto
%   permite reutilizar el mismo solver tanto para el Experimento 1
%   (parametros constantes) como para el Experimento 2 (escalon de
%   Qin en t = 500 s), sin necesidad de subfunciones ni de bifurcar el
%   bucle de integracion: simplemente se evalua p_t = p(t) si p es un
%   handle, o p_t = p si es un struct fijo.
%
% ENTRADAS:
%   h0    : vector columna 2x1, condicion inicial [h1(0); h2(0)] [m]
%   tspan : vector [t0, tf] con los limites de integracion [s]
%   dt    : paso de tiempo fijo [s]
%   p     : struct de parametros, o function handle p(t) -> struct
%
% SALIDAS:
%   t_vec : vector fila 1 x (N+1) con los instantes de tiempo [s]
%   H     : matriz 2 x (N+1) con la trayectoria [h1(t); h2(t)] [m]

    N = round((tspan(2) - tspan(1)) / dt);
    t_vec = tspan(1) + (0:N) * dt;

    H = zeros(2, N + 1);     % preasignacion: evita realloc dentro del bucle
    H(:, 1) = h0;

    es_variable = isa(p, 'function_handle');

    for k = 1:N
        t = t_vec(k);
        h = H(:, k);

        if es_variable
            p1 = p(t);
            p2 = p(t + dt/2);
            p4 = p(t + dt);
        else
            p1 = p; p2 = p; p4 = p;
        end

        k1 = campo_vectorial(h,              p1);
        k2 = campo_vectorial(h + dt/2 * k1,  p2);
        k3 = campo_vectorial(h + dt/2 * k2,  p2);
        k4 = campo_vectorial(h + dt   * k3,  p4);

        H(:, k + 1) = h + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end
