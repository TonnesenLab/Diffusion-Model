
function [C0,maxC]=SourceType(app,op)
    dy=app.ds(1);
    dx=app.ds(2);
    dz=app.ds(3);
    d3=dx*dy*dz;
    dt=app.Timestep.Value*convertToSIunit(app.Units.Value);

if op==1
    TN=app.TN.Value; % TMA Transport number RTI
    F=96485.332;%s*A/mol (Faradays Constant)
    I=app.Current.Value*10^-12; % A
    Q=I*TN/(F);%mol/s
    C0=Q*dt/(d3);% mol/m3 = 10^-3 mol/L = 1 mM 
    
    Na=6.022*10^23; %Avogadro's number
    maxnump=app.Maxnump.Value*10^18*d3;%max number of particles per m3
    maxC=(maxnump*dt/(Na*d3)); % Concentration in mM
elseif op==2
    numP=app.Particles.Value;
    Na=6.022*10^23; %Avogadro's number
    C0=(numP/(Na*d3)); % Concentration in mM
    %C0=(numP*dt/(Na*d3)); % Concentration in mM
    app.UV.Value=round(d3*10^18,5);
    
    maxnump=app.Maxnump.Value*10^18*d3;%max number of particles in a voxel
    maxC=(maxnump/(Na*d3)); % Concentration in mM

end
end