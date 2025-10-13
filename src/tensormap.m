function [YX,VU]=tensormap(FN,IN,M,idx,MS,Clim1,sc,imfolder,resfolder,alpha,Xs,Ys,rsq_th,scaleC)
    %col=[0,0.9,0.3]; %green
    %colli=[0.4,1,0.6]; %green
    % col=[0,0.4,1]; %blue
    % colli=[0.4,0.9,1]; % ligth blue
    %colli=[1,0.6,0.8]; 
    
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
    savename = fullfile(resfolder, [FN, '_MAP.tif']);


    if MS==1
        %Load source file Source(1)= X, Source(2)= Y
        Sfile=['Sources_',FN,'.xlsx'];
        Sources=xlsread(Sfile);
        for i=1:num_files

            str = split(files(i).name, '_');
            n = str2double(str{2}(1:end-4));
            F = fullfile(resfolder, files(i).name);

            load(F,'Data');
            %S0 = Sources(n,[2,1]); % just for file 20230414b
            S0 = Sources(n,[1,2]);
  
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
        
        %Normalize Chi
        maxmaxD= max(maxD);
        minG=min(A_perm)-0.0001;
        maxG=max(A_perm);
        good_new=(A_perm-minG)./(max(A_perm)-minG);

        % Threshold calculation for rsq
        if rsq_th <= 0
            figure
            h = cdfplot(good_new);
            title('Cumulative distribution of good_new')
            rsqValue = h.XData;
            cf = h.YData;
            rsq_th = rsqValue(find(cf >= 0.8, 1)); 
            
        end

        %rsq_new=rsq./rsq_th;        
        %Chi_new(Chi_new>1)=1;
         figure
         cdfplot(A_perm)
         title('CDF Area  under the curve: perimeter')
        %SMS_new=SMS./max(SMS);
%         Goodfit=1-Chi_new;
%         figure
%         cdfplot(Goodfit)
        %title('Cumulative distribution of Goodness of fit')

        VU_new=VU.*good_new';

        % Convert polar coordinates to cartesian coordinates
        %[VU(:, 2), VU(:, 1)] = pol2cart(theta, mod); % [x,y]= pol2cart(theta, mod);
        
        % Remove zero magnitude entries
        mag1 = sqrt(VU_new(:, 2).^2 + VU_new(:, 1).^2);
        non_zero_indices = mag1 ~= 0;
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

        
        %Define Rois 50x50 pixels and get principal direction in that area
%         ROIs(1,:)=[165,10,65,65]; %x, y GREEN BOX
%         my_vertices2=[165 10;165 75; 230 75; 230 10];
%         ROIs(2,:)=[152,88,152,128,270,129,270,165]; polygon
%         my_vertices=[152 88;152 128; 270 165; 270 129];
%         ROIs(3,:)=[183,227,60,60];

        
        ROIs=zeros(4,8);
        %ROIS Frame 05
        ROIs(1,:)=[54, 458, 80, 80,0,0,0,0]; %BOX
        ROIs(2,:)=[101, 33, 80, 80,0,0,0,0]; %BOX
        ROIs(3,:)=[1,363,1,333,267,469,254,491];%polygon
        ROIs(4,:)=[64,165,86,137,271,279,234,311];%polygon 
        
        %ROIS Frame 04
%         ROIs(1,:)=[224, 8, 50, 50,0,0,0,0]; %BOX
%         ROIs(2,:)=[55, 230, 50, 50,0,0,0,0]; %BOX
%         ROIs(3,:)=[156,170,135,159,191,49,210,71];%polygon
%         ROIs(4,:)=[1,105,1,130,58,188,78,174];%polygon 
        
        %ROIS Frame 03
%         ROIs(1,:)=[72, 110, 70, 70,0,0,0,0]; %BOX
%         ROIs(2,:)=[350, 210, 70, 70,0,0,0,0]; %BOX
%         ROIs(3,:)=[244,77,277,80,232,229,207,219];%polygon
%         ROIs(4,:)=[6,261,108,195,132,210,15,288];%polygon     
        %ROIS Frame 02
%         ROIs(1,:)=[243, 345, 80, 80,0,0,0,0]; %BOX
%         ROIs(2,:)=[333, 19, 80, 80,0,0,0,0]; %BOX
%         ROIs(3,:)=[298,67,343,67,312,200,277,200];%polygon
%         ROIs(4,:)=[24,240,41,211,196,292,175,324];%polygon
        % ROIS Frame 01
%         ROIs(1,:)=[127, 23, 50, 50,0,0,0,0]; %x, y GREEN BOX
%         ROIs(2,:)=[197, 150, 50, 50,0,0,0,0]; %x, y GREEN BOX
%         ROIs(3,:)=[24, 247, 50, 50,0,0,0,0]; %x, y GREEN BOX
%         ROIs(4,:)=[27,10,63,10,72,94,42,94]; %polygon
%         ROIs(5,:)=[280,50,301,50,270,130,247,130];%polygon   
%         ROIs(6,:)=[101,265,139,249,186,319,150,319];%polygon
        
        % Draw the polygons in the tensormap
        for L=1:size(ROIs,1)
            if ROIs(L,5)==0
                my_vertices=[ROIs(L,1) ROIs(L,2);ROIs(L,1)  ROIs(L,2)+ROIs(L,3);ROIs(L,1)+ROIs(L,4) ROIs(L,2)+ROIs(L,3);ROIs(L,1)+ROIs(L,4) ROIs(L,2)];
                color= 'Green';
            else
                my_vertices=[ROIs(L,1) ROIs(L,2);ROIs(L,3) ROIs(L,4);ROIs(L,5) ROIs(L,6);ROIs(L,7) ROIs(L,8)];
                color= 'Yellow';
            end
        h = drawpolygon(ax,'Position',my_vertices);
        h.Color=color;
        h.FaceAlpha=0;
        hold on
        VERT{L}=my_vertices;
        end
        
        for L=1:size(ROIs,1)
          puntas=[];
          theta=[];
          rho=[];
            
          xv=VERT{L}(:,1);
          yv=VERT{L}(:,2);
          yq=YX_nz(:,1);
          xq=YX_nz(:,2);
          [in, ~] = inpolygon(xq, yq, xv, yv);
          puntas(:,:)=VU_nz(in,:);
            
        [theta, rho] = cart2pol(puntas(:,2), puntas(:,1));
        theta(theta<0)=theta(theta<0)+3.1416*2;
        figure
        polarhistogram(theta,12,'Normalization','probability')
        title(['ROI: ' num2str(L)])
        set(gca, 'ThetaDir','clockwise');
        rlim([0 0.4])
        MRho=mean(rho);
        figure
        plot(puntas(:,2),puntas(:,1),'r.','MarkerSize',5)
        set(gca, 'YDir','reverse');

        PUNTAS{L}=puntas;
        THETA{L}=theta;
        RHO{L}=rho;
            
        end
        
        save('Tensormap_info.mat','PUNTAS','THETA');

    else
        disp('Tensormap Analysis can only be runned for multiple simulations');
    end


% if MS==1
%     Xsquares=Xs;
%     Ysquares=Ys;
%     % 6. Estimate Average tensor map
%     for xi=1:Xsquares+1
%         if xi==1
%             xp(xi)=1;
%         else
%             xp(xi)=round(real_dim(2)*((xi-1)/Xsquares));
%         end
%     end
%     for yi=1:Ysquares+1
%         if yi==1
%             yp(yi)=1;
%         else
%             yp(yi)=round(real_dim(1)*((yi-1)/Ysquares));
%         end
%     end
% 
%     xstep=round(xp(2)/2);
%     ystep=round(yp(2)/2);
%     Q=cell(Xsquares*Ysquares,1);
%     Q2=cell(Xsquares*Ysquares,1);
%     IDX=zeros(Ysquares,Xsquares);
%     for k=1:Xsquares
%         for t=1:Ysquares
%             IDX(t,k)=k*t+((Xsquares-k))*(t-1);
%         end
%     end
% Xm1=zeros(1,length(X1));
% Ym1=zeros(1,length(X1));
% Xm2=zeros(1,length(X2));
% Ym2=zeros(1,length(X2));
%     for j=1:length(X1)
%          for k=1:Xsquares
%             for t=1:Ysquares 
%                 if X1(j)>=xp(k) && X1(j)<=xp(k+1)  && Y1(j)>=yp(t) && Y1(j)<=yp(t+1)
%                     count3=IDX(t,k);
%                     s=size(Q{count3},1);
%                     idx=s+1;
%                     Q{count3}(idx,:)=UV1(j,:);
%                     Xav1(count3)=xp(k)+xstep;
%                     Yav1(count3)=yp(t)+ystep;
%                     Xm1(j)=Xav1(count3);
%                     Ym1(j)=Yav1(count3);
%                     
%                 end
%             end
%         end
%     end
%     for j=1:length(X2)
%          for k=1:Xsquares
%             for t=1:Ysquares 
%                 if X2(j)>=xp(k) && X2(j)<=xp(k+1)  && Y2(j)>=yp(t) && Y2(j)<=yp(t+1)
%                     count3=IDX(t,k);
%                     s2=size(Q2{count3},1);
%                     idx2=s2+1;
%                     Q2{count3}(idx2,:)=UV2(j,:);
%                     Xav2(count3)=xp(k)+xstep;
%                     Yav2(count3)=yp(t)+ystep;
%                     Xm2(j)=Xav2(count3);
%                     Ym2(j)=Yav2(count3);
%                     
%                 end
%             end
%         end
%     end
%  save('TM.mat','Q','UV1','Q2','UV2','Xav2','Yav2','Xav1','Yav1','rsq','rsq_new','theta','mod','theta2','mod2')
%     for z=1:length(Q)
%         if isempty(Q{z})==1
%             Q{z}=[0,0];
%         end
%         if isempty(Q2{z})==1
%             Q2{z}=[0,0];
%         end
%         u=Q{z}(:,2);
%         v=Q{z}(:,1);
%         u2=Q2{z}(:,2);
%         v2=Q2{z}(:,1);
%         [Uav(z),Vav(z),Ustdp(z),Vstdp(z),Ustdn(z),Vstdn(z)]=vector_meanandstds2(u,v,u2,v2);
% 
%     end

%     % 7. Display average tensor map
%     figure;
%     quiver(Xav1,Yav1,Uav,Vav)
%     xlim([1,real_dim(2)])
%     ylim([1,real_dim(1)])
%     f2=figure;
%     quiver(Xav1,Yav1,Uav,Vav)
%     axis image 
%     hax=gca;
%     hax.YLim=[1,real_dim(1)];
%     hax.XLim=[1,real_dim(2)];
% %     set(hax, 'YDir','reverse')
%     image(hax.XLim,hax.YLim,IM,'AlphaData',alpha);
%     colormap(gray)
%     hold on
%      sc2=sc*2;
%     quiver(Xav1,Yav1,Uav.*sc2,Vav.*sc2,0,'LineWidth',1,'AutoScale','off','Color',col,'MaxHeadSize',0.2);
%     hold on
%     Uavn=Uav.*(-1);
%     Vavn=Vav.*(-1);
%     quiver(Xav1,Yav1,Uavn.*sc2,Vavn.*sc2,0,'LineWidth',1,'AutoScale','off','Color',col,'MaxHeadSize',0.2);
% %     hold on
% %     quiver(Xav1,Yav1,Ustdp.*sc2,Vstdp.*sc2,0,'LineWidth',1,'AutoScale','off','Color',colli,'MaxHeadSize',0.2);
% %     hold on
% %     quiver(Xav1,Yav1,Ustdn.*sc2,Vstdn.*sc2,0,'LineWidth',1,'AutoScale','off','Color',colli,'MaxHeadSize',0.2);
% %     hold on
% %     quiver(Xav1,Yav1,Ustdp.*-sc2,Vstdp.*-sc2,0,'LineWidth',1,'AutoScale','off','Color',colli,'MaxHeadSize',0.2);
% %     hold on
% %     quiver(Xav1,Yav1,Ustdn.*-sc2,Vstdn.*-sc2,0,'LineWidth',1,'AutoScale','off','Color',colli,'MaxHeadSize',0.2);
% %     
%     save('TM.mat','Q','UV1','Q2','UV2','Xav2','Yav2','Xav1','Yav1','rsq','rsq_new','theta','mod','theta2','mod2','Ustdp','Vstdp','Ustdn','Vstdn')
%     %%Alernative for patch
% %     %%Alernative for patch
% %     Xp(1,:)=Xav1;
% %     Xp(2,:)=Xav1+Ustdp.*sc2;
% %     Xp(3,:)=Xav1+Ustdn.*sc2;
% % %     Xp(4,:)=Xav1+Uav.*sc2;
% %     Xpr(1,:)=Xav1;
% %     Xpr(2,:)=Xav1-Ustdp.*sc2;
% %     Xpr(3,:)=Xav1-Ustdn.*sc2;
% % %     Xpr(4,:)=Xav1+Uavn.*sc2;
% % 
% %     Yp(1,:)=Yav1;
% %     Yp(2,:)=Yav1+Vstdp.*sc2;
% %     Yp(3,:)=Yav1+Vstdn.*sc2;
% % %     Yp(4,:)=Yav1+Vav.*sc2;
% %     Ypr(1,:)=Yav1;
% %     Ypr(2,:)=Yav1-Vstdp.*sc2;
% %     Ypr(3,:)=Yav1-Vstdn.*sc2;
% % %     Ypr(4,:)=Yav1+Vavn.*sc2;
% %     hold on
% %     patch(Xp,Yp,colli,'EdgeColor',colli,'FaceAlpha',.5);
% %     hold on
% %     patch(Xpr,Ypr,colli,'EdgeColor',colli,'FaceAlpha',.5);
%     
%     hax.Visible = 'off';
%     hax.XTick=[];
%     hax.YTick=[];
%     truesize(f2);
%     sname2 = [resfolder FN '_AV.tif'];
%     H2 = getframe(gca);
%     imwrite(H2.cdata, sname2)
% 
%     f3=figure;
%     quiver(Xm1',Ym1',U1,V1)
%     axis image 
%     hax=gca;
%     hax.YLim=[1,real_dim(1)];
%     hax.XLim=[1,real_dim(2)];
% %     set(hax, 'YDir','reverse')
%     image(hax.XLim,hax.YLim,IM,'AlphaData',alpha);
%      colormap(gray)
%     hold on
%     quiver(Xm1',Ym1',U1.*sc,V1.*sc,0,'LineWidth',1,'AutoScale','off','Color',col,'MaxHeadSize',0.2);
%     hold on
%     quiver(Xm2',Ym2',U2.*sc,V2.*sc,0,'LineWidth',1,'AutoScale','off','Color',col,'MaxHeadSize',0.2);
%     hax.Visible = 'off';
%     hax.XTick=[];
%     hax.YTick=[];
%     truesize(f3);
%     sname3 = [resfolder FN '_ALL.tif'];
%     H3 = getframe(gca);
%     imwrite(H3.cdata, sname3)
%end

end