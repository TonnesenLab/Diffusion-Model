function []=CoeffofVariationofConcentration(op,file,model,ds,M,S0,U)
%CoeffofVariationofConcentration - This function calculates the coefficient
%of variation of a concentration cloud at 0.2, 0.5, 1,2,3,5 um from the
%release site. 
%Inputs: op - corresponds to the identifier of the analysis that is going to be performed.
%       file - The name of the file to be analysed
%       model - 2D or 3D concentration cloud
%       ds - pixel size
%       M - image size
%       S0 - release point 
%       U - desired concentration units 
%Outputs: The resulting Coefficient of variations for the different distances are saved in the results folder as
%CoefofVariation_analysis.mat
    
    % Load Data and scale from SI units to desired ones
    scaleC=convertFromSIunit(U);
    C=load(file);
    Y=C.Data.*scaleC;
    COV=zeros(16,6);

    % Get Measuring Points at each radius 0.2, 0.5, 1,2,3,5 um
    radius=[0.2,0.5,1,2,3,5].*10^-6; % in m 
    for i=1:length(radius)
        [MP]=MeasuringPoints(radius(i),S0,M,ds,op);
        num_MP = size(MP,1);
        if strcmp(model,'3D')==1
            idx_end= ones(num_MP,1).*size(Y,4);
        else
            Yaux=squeeze(Y(:,:,end));
            for j=1:length(MP)
                Yend(j)=Yaux(MP(j,1),MP(j,2));
            end
            
        end
        
        Ymean = mean(Yend);
        COV(:,i)=((Yend-Ymean).^2)./Ymean;
    end
    
    %Calculate statistics and save data
    SUMSQ=sum(COV);
    m=min(COV,[],'all');
    M=max(COV,[],'all');
    normCOV=(COV-m)./(M-m);
    save(['C:\Users\Paula\Documents\MATLAB\DiffusionFlux Model\results\' 'CoefofVariation_analysis.mat'],'Yend','COV','SUMSQ','normCOV');

end