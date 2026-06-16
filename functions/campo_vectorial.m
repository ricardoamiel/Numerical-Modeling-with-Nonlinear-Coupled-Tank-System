function [dh, Qover] = campo_vectorial(h, p)
%CAMPO_VECTORIAL  Campo vectorial f(h) del sistema de dos tanques acoplados.
%
% DESCRIPCION DE LA LOGICA:
%   El sistema fisico esta gobernado por:
%       dh1/dt = (Qin - c1*sqrt(h1)) / A1
%       dh2/dt = (c1*sqrt(h1) - c2*sqrt(h2)) / A2
%
%   Esta funcion es PURA (mismo h y p -> mismos dh, Qover; sin estado
%   oculto ni efectos secundarios) para poder ser invocada de forma
%   segura tanto dentro de los solucionadores temporales (RK4, Euler
%   Implicito) como dentro del bucle de Newton-Raphson multivariable.
%
%   Se implementan dos correcciones fisicas/numericas sobre el modelo
%   ideal:
%     1) RAICES IMAGINARIAS: si h_i < 0 (puede ocurrir transitoriamente
%        por overshoot numerico en pasos explicitos grandes), la altura
%        efectiva para el calculo de la raiz se trunca a 0. Esto evita
%        que sqrt() devuelva un numero complejo y mantiene el campo
%        vectorial bien definido en todo R^2.
%     2) SATURACION POR DESBORDE: si h_i ya alcanzo o supero Hmax, el
%        tanque no puede seguir llenandose. Se calcula primero la
%        derivada "libre" (sin restriccion); si dicha derivada es
%        positiva (el tanque seguiria llenando), se satura dh_i a 0 y
%        el caudal que hubiera causado ese crecimiento se reporta como
%        flujo excedente Qover_i [m^3/s] (caudal de rebose, no
%        volumen). Este segundo argumento de salida es solo
%        informativo: los solucionadores integran unicamente dh.
%
% ENTRADAS:
%   h : vector columna 2x1 = [h1; h2]  [m]
%   p : struct con campos escalares A1, A2, c1, c2, Qin, Hmax
%
% SALIDAS:
%   dh    : vector columna 2x1 = [dh1/dt; dh2/dt]      [m/s]
%   Qover : vector columna 2x1 = [Qover1; Qover2]      [m^3/s]

    % --- Paso 1: alturas efectivas no negativas para evitar raices complejas
    h_eff = max(h, 0);
    sqrt_h1 = sqrt(h_eff(1));
    sqrt_h2 = sqrt(h_eff(2));

    % --- Paso 2: campo vectorial "libre" (sin saturacion de desborde)
    dh1_libre = (p.Qin - p.c1 * sqrt_h1) / p.A1;
    dh2_libre = (p.c1 * sqrt_h1 - p.c2 * sqrt_h2) / p.A2;
    dh = [dh1_libre; dh2_libre];

    % --- Paso 3: saturacion por desborde (Hmax) y reporte de caudal excedente
    Qover = [0; 0];
    A = [p.A1; p.A2];
    sobre_limite = (h >= p.Hmax) & (dh > 0);
    if sobre_limite(1)
        Qover(1) = dh(1) * A(1);   % caudal que se pierde por rebose [m^3/s]
        dh(1) = 0;                 % la altura no puede seguir creciendo
    end
    if sobre_limite(2)
        Qover(2) = dh(2) * A(2);
        dh(2) = 0;
    end
end
