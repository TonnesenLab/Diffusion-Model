
function [rate,dt,nt,op,S0,M,ds,dif,tp,T,K,SName,Cb,clim,NameMask,unittime,unitconc, bounds,freq]=LoadParam(app)
    % Extract selected tab title
    Tab = app.TabGroup.SelectedTab.Title;

    if strcmp(Tab,'Multiple Simulations')==1
        % Extract values for multiple simulations
        unittime=app.Units_2.Value;
        unitconc=app.UnitsC_2.Value;
        SName=app.FileName_2.Value;
        NameMask=fullfile(app.imfolder, app.MaskDropDown_2.Value);
        ds=[app.dy_2.Value * 10^-6, app.dx_2.Value * 10^-6, app.dz_2.Value * 10^-6]; % in m
        K=app.Kappa_2.Value;
        M = [app.My_2.Value, app.Mx_2.Value, app.Mz_2.Value];
        dt=app.Timestep_2.Value * convertToSIunit(unittime);
        T=app.TotalDuration_2.Value * convertToSIunit(unittime);
        tp=app.Sourceduration_2.Value * convertToSIunit(unittime);
        S0 = [app.Y0_2.Value, app.X0_2.Value, app.Z0_2.Value];
        nt=T/dt;   
        rate=app.SaveRate_2.Value;
        Model = app.ModelDropDown_2.Value;
        dif=app.Dif_2.Value * 10^-10; % Diff in 3D
        Cb=0;
        freq=0;
        clim=[app.Clim1_2.Value, app.Clim2_2.Value];
        if strcmp(app.Boundaries2.Value,'scape') == 1
            bounds = 1;
        else
            bounds = 0;
        end
        
        
    else
        unittime=app.Units.Value;
        unitconc=app.UnitsC.Value;
        SName=app.FileName.Value;
        NameMask=fullfile(app.imfolder, app.MaskDropDown.Value);
        ds=[app.dy.Value * 10^-6, app.dx.Value * 10^-6, app.dz.Value * 10^-6];
        K=app.Kappa.Value;
        M = [app.My.Value, app.Mx.Value, app.Mz.Value];
        dt=app.Timestep.Value * convertToSIunit(unittime);
        T=app.TotalDuration.Value * convertToSIunit(unittime);
        tp=app.Sourceduration.Value * convertToSIunit(unittime);
        freq=app.Freq.Value;
        S0 = [app.Y0.Value, app.X0.Value, app.Z0.Value]; 
        nt=T/dt;
        rate=app.SaveRate.Value;
        Model = app.ModelDropDown.Value;
        dif=app.Dif.Value*10^-10;
        Cb=app.Cbasal.Value * 10^-3;
        clim=[app.Clim1.Value, app.Clim2.Value];
        if strcmp(app.Boundaries.Value,'scape') == 1
            bounds = 1;
        else
            bounds = 0;
        end
        
    end
    % Set operation type based on model
    if strcmp(Model,'2D')==1
        op=1;
        dif=dif * (2/3);
    elseif strcmp(Model,'2D multiple release')==1
        op=4;
        dif=dif;
    elseif strcmp(Model,'pseudo 3D')==1
        op=2;
        dif=dif;
    elseif strcmp(Model,'3D')==1
        op=3;
        dif=dif; 
    end
end
