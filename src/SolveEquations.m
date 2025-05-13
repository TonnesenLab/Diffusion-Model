function [lamda,dif_av,Data]=SolveEquations(op,rate,nt,dif,tp,K,SName,app,Cb,clim,NameMask)
% SoveEquations - This function chooses the rigth solver acoording to the
%                   user specifications

% Input: op - identifier for 2D, 3D, pseudo3D or multiple release
%               simulations
%        rate - iteration save rate. The concentration is save after the number 
%               of iterations indicated by the user.             
%        nt - number of iterations
%        dif - Free diffusion coefficient provided by the user
%        tp - duration of the release pulse
%        SName - file name introduced by user to save the data
%        app - app elements
%        Cb - basal concentraion
%        clim - concentration threshold to display the concentration cloud.
%        NameMask - Name and path of the SUSHI mask chosen by the user.

% Output: lamda - tortuosity of the image
%         dif_av - average diffusion coeffiecient after scaling.
%         Data - Matrix of the concentration of each pixel at each saved
%                   instant of time

    C0=app.Cin;
    ds2=app.ds.^2;
    S0 = app.S0;
    dt = app.dt;
    M = app.M;
    BC = app. bounds;
    aux=app.Aux.Value;
    freq=app.freq;
    
    %units
    unitC=app.unitconc;
    unitT=app.unittime;
    if op == 1 
        kapa=0;
        sc=0.90; %0.982
    elseif op == 2|| op == 4
        kapa=0.0003; % clearance in z direction 
        sc=0.9; %0.982 % scape coefficient at the boundaries. 
        C0=C0*0.2;

    elseif  op == 3
        kapa=K;
        sc=0.7; %0.97 for Agarose.tif
    end
    
    %Load diffusion matrices
    Mask=app.mask; 
    if range(Mask) ==0
        Mask=Mask.*aux;
    end
    IM2=app.IM;
    D=Mask .* dif;
    
    dif_av = mean(D(:), 'omitnan');
    lamda=sqrt(dif/dif_av);
    t(1)=0;
    t(2:nt/rate+1)=dt*rate*(1:nt/rate);

    
    % Solve the Equation for multiple simulation
    if strcmp(app.Tab,'Multiple Simulations')==1
        MS=1; 
        Source=xlsread([app.datafolder app.SourceDropDown.Value]);
        i0=app.idx0.Value;
        iend=app.idxend.Value;
        fps=app.fps_2.Value;

        for n=i0:iend
            n
            SN = sprintf('%s_%04d', SName, n);
            fname = fullfile(app.resfolder, [SN, '.mat']);
            S0=Source(n,[1,2]); % Adjust indices S0 1 = Y, 2 = X
            
            % Solve the Diffusion Equation and save the data and video
            if op==1 || op==2 
                if Mask(S0(1),S0(2))>0 
                    
                    [Data]=NumericalSolver2D(M,ds2,dt,nt,D,S0(1:2),C0,kapa,tp,rate,sc,Mask,Cb,BC);
                    save(fname,'Data','t');
                    if n==1 ||rem(n,120)==0

                        NameVideo=fullfile(app.resfolder, [SN, '.avi']);
                        SaveVideo(Data,NameVideo,IM2,nt,rate,dt,t,unitC,unitT,clim,Cb,fps,NameMask,MS);
                        save(fname,'Data','IM2','t');
                    end
                end
            end
        end
        
    %Solve the Equation for one simulation
    else
        SN=SName;
        fname=[app.resfolder,SN,'.mat'];

        if op==1 || op==2
            if size(D,3)>1
                D_new=D(:,:,round(size(D,3)/2));
            else
                D_new=D;
            end
            [Data]=NumericalSolver2D(M,ds2,dt,nt,D_new,S0(1:2),C0,kapa,tp,rate,sc,Mask,Cb, BC);
            save(fname,'Data','t','-v7.3');
        elseif op==3
            [Data]=NumericalSolver3D(M,ds2,dt,nt,D,S0(1:3),C0,kapa,tp,rate,sc,Mask,Cb);
            save(fname,'Data','t','-v7.3');
        elseif op==4
            Source_file=[app.datafolder 'Source_GABArelease.mat'];
            [Data]=NumericalSolver2D_multiple(M,ds2,dt,nt,D,Source_file,C0,kapa,tp,rate,sc,Mask,Cb, BC,freq);
            save(fname,'Data','t','-v7.3');
        end
    end
end