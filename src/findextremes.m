% Find the furthest points of the diffusion cloud from the souce
function [Xmaxd,Ymaxd,Area_perm,Area,maxD]=findextremes(AUX,S0,i)
    Ys=S0(1,1);
    Xs=S0(1,2);

    [r,c]=find(AUX);
    
    if range(r)==0 && range(c)==0
        k=0
        Xmaxd=S0(1,2);
        Ymaxd=S0(1,1);
        Area_perm=0;
        Area=0;
        maxD=0;
        
    else  
        if range(r)==0 || range(c)==0
           k=1:length(r);
           k
        else
            k=convhull(c,r);
        end
        % Calculate distances
        distance = sqrt((r - S0(1,1)).^2 + (c - S0(1,2)).^2);
        dis_perimeter=sqrt((r(k) - S0(1,1)).^2 + (c(k) - S0(1,2)).^2);
        % Find maximum distance index
        [maxD, I] = max(distance);
        Xmaxd = c(I);
        Ymaxd = r(I);
        dis_norm=distance./max(distance);
    %     dis_norm_perm=dis_perimeter./max(dis_perimeter);
        [f,x_dis]=ecdf(dis_norm);
        [f_perm,x_disperm]=ecdf(dis_perimeter);
        %make sure the cumulative distribution of the distance goes from 0 to 1
        %and the for the area under the curve is estimated from 0-1
        limit=max(x_disperm)+0.01;
        x_disperm=[0; x_disperm;limit;18];
        f_perm=[0; f_perm;0;0];

        Area_perm=trapz(x_disperm,f_perm);
        Area=trapz(x_dis,f);
    

        %     if  i<10 % i==589|| i==230|| i==218 || i==110|| i==102 %i==100 || i==438 || i==419 ||i==326 || i==275 || i==302
        % 
        % %         figure
        % %         plot(x_dis,f,'r-')
        % %         title(['CDF  Indx ',num2str(i)])
        %         
        %         figure
        %         plot(x_disperm,f_perm,'k-')
        %         title(['CDF perimeter Indx ',num2str(i)])
        %     end

             % Calculate linear regression
            if abs(Ys - Ymaxd) > abs(Xs - Xmaxd)
                P = polyfit([Ymaxd; Ys], [Xmaxd; Xs], 1);
                yfit = (min(r):max(r))';
                xfit = polyval(P, yfit);
                %rsq = 1 - (sum((polyval(P, r) - c).^2) / sum((c - mean(c)).^2));
                %SMS=sum((polyval(P, r) - c).^2);
                %Chi2=sum(((polyval(P, r) - c).^2)./c); %Pearsons Chi square test goodnesss of the fit
            else
                P = polyfit([Xmaxd; Xs], [Ymaxd; Ys], 1);
                xfit = (min(c):max(c))';
                yfit = polyval(P, xfit);
                %rsq = (1 - sum((polyval(P, c) - r).^2) / sum((r - mean(r)).^2));
                %SMS=sum((polyval(P, c) - r).^2);
                %Chi2=sum(((polyval(P, c) - r).^2)./r); %Pearsons Chi square test goodnesss of the fit
            end

            %Display cloud and fitted line
        %     if i<10  % i==589|| i==230|| i==218 || i==110|| i==102 %i==100 || i==438 || i==419 ||i==326 || i==275 || i==302
        %         figure
        %         plot(c,r,'.')
        %         hold on 
        %         plot(Xs,Ys,'r*','MarkerSize',5)
        %         xlim([min(c),min(c)+50]);
        %         ylim([min(r),min(r)+50]);
        %         set(gca, 'YDir','reverse');
        %         hold on;
        %         plot(xfit,yfit,'k--');
        %         plot(Xmaxd,Ymaxd,'r.','MarkerSize',10);
        %         title(['Indx ',num2str(i)]);
        %     end
    end

end