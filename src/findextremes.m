function [Xmaxd,Ymaxd,Area_perm,Area,maxD]=findextremes(AUX,S0,i)
% findextremes - This function finds the furthest points of the diffusion cloud from the souce
% Input: AUX - is the matrix of the concentration cloud where values below
%               a concentration limit are zero
%        S0 - release site
%        i - index of the file to be displayed ( meant to check that the line is fitting well the cloud)
% Output: Xmaxd - index of the furthest point from release site in x-axis
%         Ymaxd - index of the furthest point from release site in y-axis
%         Area_perm - area under the curve of the cumulative diftribution
%                       function of the points in the perimeter of the cloud. 
%         Area - area under the curve of the cumulative diftribution
%                       function of all the points of the cloud. 
%         maxD - maximum distance to the release site.

    Ys=S0(1,1);
    Xs=S0(1,2);

    [r,c]=find(AUX);
    if range(r)==0 || range(c)==0
       k=1:length(r);
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

    [f,x_dis]=ecdf(dis_norm);
    [f_perm,x_disperm]=ecdf(dis_perimeter);
    
    %make sure the cumulative distribution of the distance goes from 0 to 1
    %and the for the area under the curve is estimated from 0-1
    limit=max(x_disperm)+0.01;
    x_disperm=[0; x_disperm;limit;18];
    f_perm=[0; f_perm;0;0];
    
    Area_perm=trapz(x_disperm,f_perm);
    Area=trapz(x_dis,f);
    
    if  i<10 
        
        figure
        plot(x_disperm,f_perm,'k-')
        title(['CDF perimeter Indx ',num2str(i)])
    end
    
     % Calculate linear regression
    if abs(Ys - Ymaxd) > abs(Xs - Xmaxd)
        P = polyfit([Ymaxd; Ys], [Xmaxd; Xs], 1);
        yfit = (min(r):max(r))';
        xfit = polyval(P, yfit);

    else
        P = polyfit([Xmaxd; Xs], [Ymaxd; Ys], 1);
        xfit = (min(c):max(c))';
        yfit = polyval(P, xfit);

    end
    
    %Display cloud and fitted line
    if i<10  
        figure
        plot(c,r,'.')
        hold on 
        plot(Xs,Ys,'r*','MarkerSize',5)
        xlim([min(c),min(c)+50]);
        ylim([min(r),min(r)+50]);
        set(gca, 'YDir','reverse');
        hold on;
        plot(xfit,yfit,'k--');
        plot(Xmaxd,Ymaxd,'r.','MarkerSize',10);
        title(['Indx ',num2str(i)]);
    end
end