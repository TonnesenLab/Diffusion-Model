function [MPfinal,Anglesfinal]=MeasuringPoints(radius,source,M,ds,op)
    
    % Define number of measuring points
    angles = linspace(0, 2*pi, 16);
%     if op==6
%         angles=linspace(0,1.5*pi,4);
%     end
    numbMP = length(angles);
    
    % Preallocate memory and define parameters
    x = zeros(numbMP, 1);
    y = zeros(numbMP, 1);
    z = zeros(numbMP, 1);
 
    % Generate points on concentric circle.
    for i = 1:numbMP
        x(i) = source(2) + round((radius /ds(2)) * cos(angles(i)));
        y(i) = source(1) + round((radius /ds(1)) * sin(angles(i)));
        z(i) = source(3);
    end

    valid_indices = x > 0 & y > 0 & x <= M(2) & y <= M(1);
    MPfinal(1:sum(valid_indices(:)), :) = [y(valid_indices), x(valid_indices), z(valid_indices)];
    Anglesfinal = angles(valid_indices);

%     figure(2)
%     scatter(MPfinal(:,1),MPfinal(:,2),'o');
%     hold on
end

