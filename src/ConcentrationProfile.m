function ConcentrationProfile(file,molec,model,dx,Mx,S0,U,UT,Cb)

    % Transform units from SI to desired ones
    scaleC = convertFromSIunit(U);
    scaleT = convertFromSIunit(UT);
    Cb=Cb*scaleC; 
    C=load(file);
    Y=C.Data;
    t=C.t.*scaleC; 
    T=C.t(end) * scaleT;
    
    % Generate parameters
    L=((0:1:Mx-1)-(Mx-1)/2).*dx.*10^6; %Concentration profile length in µm
    Ls = [0.2, 0.5, 1, 2, 4]; %Distance in um
    L_idx=round(Ls./(dx*10^6))+S0(1,2);
    timepoints = T * [0.1/4, 1/4, 2/4, 3/4, 1];
    dis = [floor(L(1)), floor(L(1))/2, 0, ceil(L(end))/2, ceil(L(end))];
    tp=round(timepoints,2);
    th=Cb+Cb*0.05;
    
    % Extract and scale data
    if strcmp(model,'3D')==1
        Y1(:,:)=squeeze(Y(S0(1,1),:,S0(1,3),:));
        Y2(:,:)=squeeze(Y(S0(1,1),L_idx,S0(1,3),:));
    else
        Y1(:,:)=squeeze(Y(S0(1,1),:,:));
        Y2(:,:)=squeeze(Y1(S0(1,1),L_idx,:));
    end
    Y1(Y1<10^-4)=10^-4;
    Y1scale = Y1.* scaleC;
    Y2scale = Y2.* scaleC;
    Y2scale(Y2scale<th)=0;
    
    % Plot 2D surf concentration over distance and time 
    % create meshgrid
    [Ty,Lx]=meshgrid(t,L);
    
    figure 
    surf(Ty,Lx,Y1scale, 'EdgeColor', 'none');
    
    % Set plot properties
    set(gca, 'ColorScale', 'log', 'ZScale', 'log','XLim', [-0.08 tp(5)], 'XTick', tp, ...
        'YLim', [-4.8 4.5], 'YTick', dis);
    colormap jet;
    ax = gca;
    ax.Units = 'centimeters';
    pos=ax.Position;
    ax.Position = [pos(1) pos(2) 4 4];
    ax.XLabel.String = 'Distance (µm)';
    ax.YLabel.String = ['Time [' UT ']'];
    ax.ZLabel.String = [molec ' [' U ']'];
    
    % 3D Plot concentration profile over distance and time and basal
    % concentration
    
    %create meshgrid
    [L2y,T2x]=meshgrid(Ls,t);
    [Tx,Ly]=meshgrid(L,t);
    
    figure
    plot3(L2y,T2x,Y2scale(:,:),'-w','LineWidth',1);
    hold on
    surf(Ly,Tx,Y1scale', 'EdgeColor', 'none');
    colormap jet
    
    % Set plot properties
    ax2=gca;
    ax2.Units='centimeters';
    pos=ax2.Position;
    ax2.Position=[pos(1) pos(2) 4 4];
    set(gca, 'YLim', [t(1) t(end)], 'YTick', tp, 'XLim',[0 Ls(end)+0.1], ...
        'ColorScale','log', 'ZScale','log', 'Zlim', [10^-4,500]);
    caxis([10^-4 100]);
    
end