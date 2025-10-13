function [color,MOLEC] = colorcodefunction (model, molec)
    % Set legend color based on molecule and model
    colorMatrix = [1, 0.4, 0.4; 1, 0.4, 0.4; 1, 0.4, 0;...
        0, 0.8, 0.6; 0.2, 0.6, 0.6; 0, 0.4, 0.6; ...
        0.6, 0.8, 0.8; 0.2, 0.6, 0.6; 0, 0.4, 0.6; ...
        0, 0.8, 1; 0.2, 0.6, 0.6; 0, 0.4, 0.6; ...
        0.4,0,0.8; 0.4,0.4,0.8; 0.4,0.8,0.8];
    
    switch molec
        case 'TMA'
            MOLEC=1;
            moleculeColorIdx = 0;
        case 'dex70k'
            MOLEC=2;
            moleculeColorIdx = 3;
        case 'dex3k'
            MOLEC=3;
            moleculeColorIdx = 6;
        case 'BSA'
            MOLEC=3;
            moleculeColorIdx = 9;
        case 'GLUT'
            MOLEC=4;
            moleculeColorIdx = 12;
    end
     switch model
        case '3D'
            modelColorIdx = 1;
        case '2D'
            modelColorIdx = 2;
        case 'pseudo3D'
            modelColorIdx = 3;
     end
    color = colorMatrix((modelColorIdx + moleculeColorIdx), :);
end