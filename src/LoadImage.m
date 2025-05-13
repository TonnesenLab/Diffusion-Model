function [alpha,IM2,Mask]=LoadImage(app)
% LoadImage - This function loads the chosen image by the user with the
%               specified dimensions. It also displays the image
% Input: app - the app files
% Output: alpha - estimates the volume fraction of that image
%         IM2 - Image adjusted to the chosen size. 
%         Mask - Mask Image (grayscale) adjusted to the chosen size. The
%                   pixel values of the Mask determine the probability of diffusion. 

    Tab=app.Tab;
    aux=app.Aux.Value;
    
    if strcmp(Tab,'Multiple Simulations')==1
        S0 = [1, 1];
        Idx = [app.Y0_2.Value, app.X0_2.Value];
        My=app.My_2.Value;
        Mx=app.Mx_2.Value;
        Name=[app.imfolder,app.ImageDropDown_2.Value];
        NameMask=[app.imfolder,app.MaskDropDown_2.Value];
    else
        Idx = [app.idx1.Value, app.idx2.Value];
        if app.op ==4
            load([app.datafolder 'Source_GABArelease.mat']);
            S0(:,1)=SourceY;
            S0(:,2)=SourceX;
        else
            S0 = [app.Y0.Value, app.X0.Value];
        end
        My=app.My.Value;
        Mx=app.Mx.Value;
        Name=[app.imfolder,app.ImageDropDown.Value];
        NameMask=[app.imfolder,app.MaskDropDown.Value];
    end

    dim = [My-1, Mx-1];
    [Mask0, zstack] = OpenImage(NameMask, aux);
    

    if zstack 
        IM=imread(Name,app.Z0.Value);
        IM2=IM(Idx(1):Idx(1)+dim(1),Idx(2):Idx(2)+dim(2));
        Mask=Mask0(Idx(1):Idx(1)+dim(1),Idx(2):Idx(2)+dim(2),:);
        %image(Mask(:,:,ceil(Mz/2)),'CDataMapping','scaled');
    else
        IM=imread(Name);
        s=size(IM);
        if size(s,2)>2
            IM2=IM(Idx(1):Idx(1)+dim(1),Idx(2):Idx(2)+dim(2),:); 
        else
            IM2=IM(Idx(1):Idx(1)+dim(1),Idx(2):Idx(2)+dim(2));
        end
        Mask=Mask0(Idx(1):Idx(1)+dim(1),Idx(2):Idx(2)+dim(2));
        %image(Mask,'CDataMapping','scaled');
    end

    Mask_a=imread(NameMask);
    Mask_a=Mask_a(Idx(1):Idx(1)+dim(1),Idx(2):Idx(2)+dim(2));
    th = round((max(Mask_a(:)) - min(Mask_a(:))) / 2);
    [counts, ~] = imhist(Mask_a, 255);
    total=sum(counts);
    count_ECS=sum(counts(th+1:end));
    alpha= count_ECS / total;

    if strcmp(Tab,'Multiple Simulations')==1
        Sources=xlsread([app.datafolder app.SourceDropDown.Value]);% Y=Sources(:,1), X=Sources(:,2)
        app.Ax3.XTick=[];
        app.Ax3.YTick=[];
        image(app.Ax3,IM2,'CDataMapping','scaled');
        axis(app.Ax3,'image');
        colormap(app.Ax3,'gray');
        I0=app.idx0.Value;
        Iend=app.idxend.Value;
        hold(app.Ax3, 'on')
        plot(app.Ax3,Sources(I0:Iend,2),Sources(I0:Iend,1),'r.','Markersize',5);

    else
        app.Ax2.XTick=[];
        app.Ax2.YTick=[];
        image(app.Ax2,IM2); %'CDataMapping','scaled'
        %app.Ax2.CLim=[0,255];
        colormap(app.Ax2,'gray');
        axis(app.Ax2,'image');
        hold(app.Ax2, 'on') 
        if app.op==4
            plot(app.Ax2,S0(:,2),S0(:,1),'r.','Markersize',8);
        else
            plot(app.Ax2,S0(2),S0(1),'r.','Markersize',8);
        end
    end

end