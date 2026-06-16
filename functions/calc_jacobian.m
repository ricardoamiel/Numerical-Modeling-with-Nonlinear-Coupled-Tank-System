function J = calc_jacobian(h, p)
%CALC_JACOBIAN  Jacobiana analitica J = df/dh del campo vectorial f(h).
%
% DESCRIPCION DE LA LOGICA:
%   Derivando analiticamente las dos ecuaciones del modelo respecto de
%   h1 y h2 se obtiene una matriz triangular inferior 2x2:
%
%       J = [ -c1/(2*A1*sqrt(h1))         0                  ]
%           [  c1/(2*A2*sqrt(h1))   -c2/(2*A2*sqrt(h2))       ]
%
%   Esta Jacobiana analitica (no aproximada por diferencias finitas) es
%   reutilizada en dos contextos distintos del proyecto:
%     a) newton_multivar.m, para hallar el estado estacionario h* tal
%        que f(h*) = 0.
%     b) solver_implicit.m, dentro del Newton-Raphson embebido en cada
%        paso de Euler Implicito, donde se necesita la matriz
%        I - dt*J para linealizar el sistema algebraico no lineal.
%
%   Para evitar division entre cero cuando h_i = 0 (singularidad fisica
%   de la raiz cuadrada, tipica en el instante inicial con tanques
%   vacios), se introduce un piso numerico h_min = 1e-8. Esto regulariza
%   la matriz sin alterar de forma apreciable el resultado en el resto
%   del dominio.
%
% ENTRADAS:
%   h : vector columna 2x1 = [h1; h2]  [m]
%   p : struct con campos escalares A1, A2, c1, c2 (Qin y Hmax no
%       afectan la Jacobiana, ya que el termino Qin es constante y la
%       saturacion se trata aparte como un evento de contorno)
%
% SALIDAS:
%   J : matriz 2x2 = df/dh evaluada en el estado h

    h_min = 1e-8;               % piso numerico para evitar division por 0
    h1 = max(h(1), h_min);
    h2 = max(h(2), h_min);

    sqrt_h1 = sqrt(h1);
    sqrt_h2 = sqrt(h2);

    J11 = -p.c1 / (2 * p.A1 * sqrt_h1);
    J12 = 0;
    J21 =  p.c1 / (2 * p.A2 * sqrt_h1);
    J22 = -p.c2 / (2 * p.A2 * sqrt_h2);

    J = [J11, J12; J21, J22];
end
