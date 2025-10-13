function [lamda,dif_av,Data]=SolveEquations(op,rate,nt,dif,tp,T,K,SName,app,Cb,clim,NameMask)
    C0=app.Cin;
    f=app.DiffusionModel;
    ds2=app.ds.^2;
    S0 = app.S0;
    dt = app.dt;
    M = app.M;
    BC = app. bounds;
    alpha=app.Alpha.Value;
    aux=app.Aux.Value;
    freq=app.freq;
    
    %units
    unitC=app.unitconc;
    unitT=app.unittime;
    if op == 1 
        kapa=0;
        sc=0.90; %0.982
    elseif op == 2|| op == 4
        %kapa=15e3; % clearance in z direction kapa depends on t kapa=15e3; 
        sc=0.9; %0.982 % scape coefficient at the boundaries. 
        kapa=375;
        C0=C0*0.2;
%         if strcmp(app.ImageDropDown.Value,'Agarose.tif')==1
%             x=app.Dif.Value * 10^-1;%dif value 10^-9 m2/s
%             C0=C0 * (0.353 * x^(-0.101));
%         end
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
    max(max(Mask))
    max(max(D))
    
    dif_av = mean(D(:), 'omitnan');
    lamda=sqrt(dif/dif_av);
    t(1)=0;
    t(2:nt/rate+1)=dt*rate*(1:nt/rate);

    
    % Solve the Equation for multiple simulation
    if strcmp(app.Tab,'Multiple Simulations')==1
        MS=1; 
        alpha=app.Alpha_2.Value;
        %IDX=xlsread(app.SourceDropDown.Value);
        Source=xlsread(app.SourceDropDown.Value);
        i0=app.idx0.Value;
        iend=app.idxend.Value;
        fps=app.fps_2.Value;

        for n=i0:iend
            n
            SN = sprintf('%s_%04d', SName, n);
            fname = fullfile(app.resfolder, [SN, '.mat']);
            S0=Source(n,[1,2]); % Adjust indices S0 1 = Y, 2 = X
            %S03=[Source(n,:),5];
            
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
    %         elseif op==3
    %             if Mask(S0(1),S0(2),S0(3))>0
        %             [Data]=NumericalSolver3D(M,ds2,dt,nt,D,S0(1:3),C0,kapa,tp,T,rate,f,alpha,CMAX,sc,Mask,Cb);
        %             save(fname,'Data','t','-v7.3'); 

    %                 Data3D1=squeeze(Data(:,:,round(M(3)/2),:));
    %                 NameVideo1=fullfile(app.resfolder, [SN, '_z2.avi']);
    %                 SaveVideo(Data3D1,NameVideo1,IM2,nt,rate,dt,t,unitC,unitT,clim,Cb,fps,NameMask,MS);    
    %                 save(fname,'Data','IM2','t','-v7.3');
    %             end
            end
        end
        
    %Solve the Equation for one simulation
    else
        SN=SName;
        fname=[app.resfolder,SN,'.mat'];

        if op==1 || op==2
            
            [Data]=NumericalSolver2D(M,ds2,dt,nt,D,S0(1:2),C0,kapa,tp,rate,sc,Mask,Cb, BC);
            save(fname,'Data','t','-v7.3');
        elseif op==3
            [Data]=NumericalSolver3D(M,ds2,dt,nt,D,S0(1:3),C0,kapa,tp,T,rate,f,alpha,sc,Mask,Cb, BC);
            save(fname,'Data','t','-v7.3');
        elseif op==4
            Source_file='Source_GABArelease_ROI1.mat';
            [Data]=NumericalSolver2D_multiple(M,ds2,dt,nt,D,Source_file,C0,kapa,tp,rate,sc,Mask,Cb, BC,freq);
            save(fname,'Data','t','-v7.3');
        end
    end
end