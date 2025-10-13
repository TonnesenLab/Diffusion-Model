function []=CoeffofVariationofConcentration(op,file,model,ds,M,S0,U)
    
    % Load Data and scale from SI units to desired ones
    scaleC=convertFromSIunit(U);
    C=load(file);
    Y=C.Data.*scaleC;
    COV=zeros(16,6);
    % Get Measuring Points at each radius 0.5, 1, 2,5, 10 um
    radius=[0.2,0.5,1,2,3,5].*10^-6; % in m 
    for i=1:length(radius)
        [MP]=MeasuringPoints(radius(i),S0,M,ds,op);
        i
        MP
        size(Y)
        num_MP = size(MP,1);
        if strcmp(model,'3D')==1
            idx_end= ones(num_MP,1).*size(Y,4);
            linear_idx = sub2ind(size(Y), MP(:, 1), MP(:, 2), MP(:,3), idx_end);
        else
            Yend=squeeze(Y(MP(:,1),MP(:,2),end));
            
        end
        Yend
        Ymean = mean(Yend');
        Ymean
        COV(:,i)=((Yend-Ymean).^2)./Ymean;
    end
    
    %Calculate statistics and save data
    SUMSQ=sum(COV);
    m=min(COV,[],'all');
    m
    M=max(COV,[],'all');
    M
    normCOV=(COV-m)./(M-m);
    save('Yfinal.mat','Yend','COV','SUMSQ','normCOV');

end