function angulo_medio = angulo_medio_axial(angulos)
% angulo_medio_axial  Calcula la orientación media axial (0–180°)
%   angulos: vector de ángulos en grados (0–180)
%   angulo_medio: ángulo medio axial en grados

    % Asegurar que es un vector columna
    angulos = angulos(:);

    % 1.Doblar los ángulos
    angulos_doblados = 2 * angulos;

    % 2.Convertir a radianes
    ang_rad = deg2rad(angulos_doblados);

    % 3.Calcular componentes x e y promedio
    x = mean(cos(ang_rad));
    y = mean(sin(ang_rad));

    % 4.Calcular ángulo medio doblado (en radianes)
    angulo_medio_doblado = atan2(y, x);

    % 5.Convertir a grados y desdoblar
    angulo_medio = rad2deg(angulo_medio_doblado) / 2;

    % 6.Ajustar al rango 0–180
    if angulo_medio < 0
        angulo_medio = angulo_medio + 180;
    end
end