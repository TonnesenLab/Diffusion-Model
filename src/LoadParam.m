function [rate,dt,nt,op,S0,M,ds,dif,tp,T,K,SName,Cb,clim,NameMask,unittime,unitconc, bounds,freq]=LoadParam(app)
% LoadParam - This function loads all the parameters the user has
%               introduced in the app. 
% Input: app - app files
% Output: rate - iteration save rate. The concentration is save after the
%               number of iterations indicated by the user. 
%       dt - time step for the iterations
%       nt - number of iterations
%       op - model identifier 2D, 3D, psuedo3D or pseudo 3D with multiple release sites
%       S0 - release site
%       ds - pixel size 
%       M - Image size
%       dif - Free Diffusion coefficient
%       tp - duration of the release pulse
%       T - duration of the simulation
%       K - Clearance due to uptake ( not scape or diffusion in z)
%       SName - file name chosen by user
%       Cb - basal concentraion
%       clim - Concentration limits
%       NameMask - name and path of the mask sushi image
%       unittime - units of time chosen
%       unitconc - units of concentration chosen
%       bounds - binary value that indicates toroidal (0) or scape (1) boundary
%               conditions
%       freq - Frequency of release in Hz

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
        dif_in=app.Dif_2.Value * 10^-10; % Diff in 3D
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
        dif_in=app.Dif.Value*10^-10;
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
        dif=dif_in * (2/3);
    elseif strcmp(Model,'2D multiple release')==1
        op=4;
        dif=dif_in;
    elseif strcmp(Model,'pseudo 3D')==1
        op=2;
        dif=dif_in;
    elseif strcmp(Model,'3D')==1
        op=3;
        dif=dif_in; 
    end
end
