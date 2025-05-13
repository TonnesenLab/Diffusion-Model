function TIFFgenerator(resfolder,fname,frames,Idx,model,Clim,unitC)
% Tiffgeneration - This function generates a tiff image from a concentration cloud without background

% Input: resfolder - folder where results are stored
%        fname - file name to open
%        frames - number of frames of the concentration profile to open.
%        Idx - first and last frames to include
%        model - 3D or 2D (pseudo 3D)
%        Clim - concentration limits
%        unitC - desired concentration units. 

    Filename = fullfile(resfolder, [fname, '.mat']);
    load(Filename);
    IMname = fullfile(resfolder, [fname, '_']);
    Mz=size(Data, 3);
    Clim = Clim * convertToSIunit(unitC);

    f = round(linspace(Idx(1), Idx(2), frames));

    if strcmp(model,'3D')==1
        Idx(3)=ceil(Mz/2);
        for i=1:length(f)
            IM=Data(:,:,Idx(3),f(i));
            IM(IM<Clim(1))=0;

            if i==1
                imwrite(IM,IMname)
            else
                imwrite(IM,IMname,'WriteMode','append');
            end
        end
    else
        class(Data)
        for i=1:length(f)
            IM=Data(:,:,f(i));
            imname=[IMname,'_',num2str(i),'.tif'];
            fig=figure;
            h=axes;
            % Ajustar el tamaño de la figura a 400x400 píxeles
            set(fig, 'Units', 'pixels', 'Position', [100, 100, 400, 400]);
            set(h, 'Units', 'normalized', 'Position', [0 0 1 1]);

            image(h,IM,'CDataMapping','scaled');
            axis(h,'image');
            set(h,'ColorScale', 'log', 'XTick', [], 'YTick', []);
            h.CLim=[0.00001, 0.1];

            % Guardar la imagen en TIFF
            print(fig, imname, '-dtiff', '-r100'); % -r100 especifica 400x400 píxeles porque 400x400 es un múltiplo de 100 DPI
    
        end
    end
end
