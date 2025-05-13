function [YX,VU]=tensormap(FN,IN,M,MS,Clim1,sc,imfolder,resfolder,alpha,scaleC)
% tensormap - This function calculates the principal direcction of each diffusion patter 
%               for multiple simulations generating a tensormap

% Input: FN - Root name of files to be analyse
%        IN - SUSHI background Image name and path were the simulations
%               were run
%        M - Image size
%        MS - Binary value indicate multiple simulations (1) or single
%               simulation (0)
%        Clim1 - Concentration limits to analyse the concentration clouds
%        sc - scale value for the tensors
%        imfolder - path to SUSHI images
%        resfolde - folder were the resulting data is stored
%        alpha - transparencey of the background image
%        scaleC - scale factor to adjust concentration units to desired ones

% Output: YX - release sites of all the simulations
%         VU - arrow tips for the principal direction


    alpha=alpha/100;

    %Load Image
    im = imread(fullfile(imfolder, IN));
    dim=[M(1)-1, M(2)-1];

    
    % Open concentration files
    files = dir(fullfile(resfolder, [FN, '*.mat']));
    
    num_files = length(files);
    VU = zeros(num_files, 2);
    YX = zeros(num_files, 2);

    SMS = zeros(num_files, 1);
    col=[1,0.2,0.8]; %pink 
    


    if MS==1
        %Load source file Source(1)= X, Source(2)= Y
        Sfile=['Sources_',FN,'.xlsx'];
        Sources=xlsread(Sfile);
        for i=1:num_files

            str = split(files(i).name, '_');
            n = str2double(str{2}(1:end-4));
            F = fullfile(resfolder, files(i).name);
            load(F,'Data');
            S0 = Sources(n,[2,1]); % just for file 20230414b
  
            YX(i, :) = S0;

            AUX = squeeze(Data(:,:,end)).*scaleC;
            AUX(AUX < Clim1) = 0;

            if max(AUX(:)) == 0
                VU(i, :) = [0, 0];
                SMS(i)=0;
            else
                [Xmax, Ymax,A_perm(i),Area(i),maxD(i)] = findextremes(AUX, S0,i);
                %  convert to polar coordinates (theta -pi to pi)
                %[theta(i), rho(i)] = cart2pol(Xmax - S0(2), Ymax - S0(1));
                VU(i, 2)=(Xmax - S0(2));
                VU(i, 1)=(Ymax - S0(1));
   
            end
        end
        
        %Calculate goodness of the directionality.
        minG=min(A_perm)-0.0001;
        good_new=(A_perm-minG)./(max(A_perm)-minG);

         figure
         cdfplot(A_perm)
         title('CDF Area  under the curve: perimeter')

        VU_new=VU.*good_new';
        
        % Remove zero magnitude entries
        mag1 = sqrt(VU_new(:, 2).^2 + VU_new(:, 1).^2);
        VU_nz = VU_new(mag1 ~= 0,:);
        YX_nz = YX(mag1 ~= 0,:);
        
        % Save data
        save('TM.mat','VU_nz','good_new','Area','A_perm','Xmax','Ymax','YX_nz','Clim1','S0');
        
        % Display tensor map on top of the image
        figure;
        image([1, dim(2)], [1, dim(1)],  im, 'AlphaData', alpha);
        colormap(gray)
        axis image 
        xlim([1, dim(2)]);
        ylim([1, dim(1)]);
        hold on;
        quiver(YX_nz(:,2), YX_nz(:,1), VU_nz(:,2).* sc, VU_nz(:,1).* sc, 0, 'LineWidth', 1, 'AutoScale', 'off', 'Color', col, 'MaxHeadSize', 0.5);
        ax = gca;
        ax.Visible = 'off';
        truesize;
        hold on
        my_vertices=[152 88;152 128; 270 165; 270 129];
        h = drawpolygon(ax,'Position',my_vertices);
        h.Color= 'Yellow';
        h.FaceAlpha=0;
        hold on
        my_vertices2=[165 10;165 75; 230 75; 230 10];
        h = drawpolygon(ax,'Position',my_vertices2);
        h.Color= 'Green';
        h.FaceAlpha=0;

        %Define Rois 50x50 pixels and get principal direction in that area
%         ROIs(1,:)=[165,10,65,65]; %x, y GREEN BOX
%         ROIs(2,:)=[152,88,152,128,270,129,270,165]; polygon
%         ROIs(3,:)=[183,227,60,60]; 
        x1=165;
        xend=x1+65;
        y1=10;
        yend=y1+65;
        
        x2=152;
        xend2=x2+118;
        
        x3=183;
        xend3=x3+60;
        y3=227;
        yend3=y3+60;
        
        count1=0;
        count2=0;
        count3=0;
        for i=1:length(YX_nz)
            y=YX_nz(i,1);
            x=YX_nz(i,2);
            y2=0.429*x + 81.38;
            yend2=0.429*x + 111.37;
            if (y1<=y)&&(y<=yend)&&(x1<= x)&&(x <= xend)
                count1=count1+1;
                puntas_1(count1,:)=VU_nz(i,:); % GREEN BOX
            elseif (y2<=y)&&(y<=yend2)&&(x2<= x)&&(x <= xend2)
                count2=count2+1;
                puntas_2(count2,:)= VU_nz(i,:);% YELLOW
            elseif(y3<=y)&&(y<=yend3)&&(x3<= x)&&(x <= xend3)
                count3=count3+1;
                puntas_3(count3,:)= VU_nz(i,:);
            end
        end
        
        [theta1, ~] = cart2pol(puntas_1(:,2), puntas_1(:,1));
        theta1(theta1<0)=theta1(theta1<0)+3.1416*2;
        
        figure
        polarhistogram(theta1,12,'Normalization','probability')
        set(gca, 'ThetaDir','clockwise');
        rlim([0 0.4])

        figure
        plot(puntas_1(:,2),puntas_1(:,1),'r.','MarkerSize',5)
        set(gca, 'YDir','reverse');
        
        [theta2, ~] = cart2pol(puntas_2(:,2), puntas_2(:,1));
        theta2(theta2<0)=theta2(theta2<0)+3.1416*2;
        
        figure
        polarhistogram(theta2,12,'Normalization','probability')
        set(gca, 'ThetaDir','clockwise');
        rlim([0 0.4])
        
        figure
        plot(puntas_2(:,2),puntas_2(:,1),'r.','MarkerSize',5)
        set(gca, 'YDir','reverse');
        
        save([resfolder 'Tensromap_info.mat'],'puntas_1','puntas_2','puntas_3','theta1','theta2');

    else
        disp('Tensormap Analysis can only be runned for multiple simulations');
    end




end