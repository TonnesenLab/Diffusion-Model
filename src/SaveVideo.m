%%% Create Video and save it %%%
function SaveVideo(Data,NameVideo,ECS,nt,rate,dt,originalTime,unitC,unitT,clim,Cbasal,fps,NameMask,MS)

    % Define scaling factors for time and concentration units
    scale=convertFromSIunit(unitC);
    scaleT=convertFromSIunit(unitT);
    
    % Convert basal concentration units from uM to mM amd then to desired
    % units
    Cb=Cbasal*10^-3*scale; 

    % Determine transparency limit
    transplim = max(Cb, clim(1));
    
    % Precompute constants and indices
    step = round(nt / rate);
    t(1)=0;
    t(2:step)=dt*rate*(1:step-1).*scaleT;
    idx = (1:round(size(Data,3)/step):size(Data,3));
    r = double((max(t) < 10));
    T=numel(originalTime);
    %idx2 = [2, round(0.33 * T), round(0.66 * T), T];
    %idx2 = [2, ceil(0.25 * T)+1, ceil(0.5 * T)+1, ceil(0.75 *T), T];
    %idx2 = [2, round(0.2 * T)+1, round(0.4 * T)+1, round(0.6 * T)+1, round(0.8 *T)+1, T];
    %idx2 = [21, 521]; %COV_1
    %idx2 = [5, 76,151,226]; 
    %idx2 = [2, 189,376,564]; 
    idx2 = [2, 201,401,601];% 0, 4, 8, 12 ms
    if strcmp(NameMask,'Agarose.tif')==1
        cl='white';

    else
        cl='gray';
    end
    
    % Display background ECS image
    figure('Color', 'w');
    ax = axes;
    h = axes;
    image(ax, ECS, 'CDataMapping', 'scaled');
    axis(ax, 'image');
    ax.CLim = [0, 255];
    colormap(ax, cl);
    set(ax, 'XTick', [], 'YTick', []);
    axpos=get(ax,'Position');
    set([ax,h],'Position',axpos);
    
    % Create video object
    writerObj = VideoWriter(NameVideo);
    writerObj.FrameRate = fps;
    open(writerObj);
    
    for i = idx
        IM = Data(:, :, i)*scale;
        TransparencyData = IM > transplim;
        hold(ax,'on'); 
        image(h,IM,'AlphaData', TransparencyData,'CDataMapping','scaled');
        axis(h,'image');
        set(h, 'ColorScale', 'log', 'XTick', [], 'YTick', []);
        h.Visible = 'off';
        h.CLim=clim;
        colormap(h,'jet');
        time = num2str(round(originalTime(i)*scaleT, r));
        ax.Title.String = ['Time = ', time, unitT];
        ax.Title.FontSize = 20;
        ax.Title.FontName = 'Arial';
        ax.Title.FontWeight = 'normal';
        
        cb = colorbar('FontSize', 12);
        cb.Ruler.Scale = 'log';
        cb.Label.String = unitC;
        cb.Label.FontSize = 14;
        cb.Position = [0.84, 0.13, 0.03, 0.78];

        axis image;
        drawnow;

        % Write frame to video
        writeVideo(writerObj, getframe(gcf)); 
    end
    close(gcf);
    % Save image montage if MS is 0
    if MS == 0
        for i = idx2
            IM = Data(:, :, i) * scale;
            TransparencyData = IM > transplim;

            figure('Color', 'w');
            ax = axes;
            h = axes;
            image(ax, ECS, 'CDataMapping', 'scaled');
            axis(ax, 'image');
            set(ax, 'XTick', [], 'YTick', []);
            ax.CLim = [0, 255];
            colormap(ax, cl);
            axpos=get(ax,'Position');
            set([ax,h],'Position',axpos);
            
            hold(ax,'on'); 
            image(h,IM,'AlphaData', TransparencyData,'CDataMapping','scaled');
            axis(h,'image');
            set(h, 'ColorScale', 'log', 'XTick', [], 'YTick', []);
            h.Visible = 'off';
            h.CLim=clim;
            colormap(h,'jet');
            
            time = num2str(round(originalTime(i) * scaleT, r));
            ax.Title.String = ['Time = ', time, unitT];
            ax.Title.FontSize = 20;
            ax.Title.FontName = 'Arial';
            ax.Title.FontWeight = 'normal';

            cb = colorbar('FontSize', 12);
            cb.Ruler.Scale = 'log';
            cb.Label.String = unitC;
            cb.Label.FontSize = 14;
            cb.Position = [0.84, 0.13, 0.03, 0.78];

            axis image;
            drawnow;
        end
    end

    % Close the video writer object
    close(writerObj);
end