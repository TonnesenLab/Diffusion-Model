%%% NUMERICAL METHOD TO SOLVE THE 2D Diffusion Equation for multiple release points at the same time %%% 
function [C_store]=NumericalSolver2D_multiple(M,ds2,dt,nt,Diff,Source_file,Q,kapa,tp,rate,sc,Mask,Cb, BC,freq)
    tic
    %Extract source coordinates
    Mx=M(2);
    My=M(1);
    dy2=ds2(1);
    dx2=ds2(2);
    
    load(Source_file);
    conjunto_indices = sub2ind([My, Mx], SourceY, SourceX);
    
    load('GABA_uptake_points2.mat');
    %uptake_indices = sub2ind([My, Mx], SourceY, SourceX);
    
 
    aux=1;
    pareja_indice=zeros(Mx*My,1);
    for i=1:My
        for j=1:Mx
            pareja_indice(aux) = sub2ind([My, Mx], i, j);
            aux=aux+1;
        end
    end
    %is_uptake=ismember(pareja_indice, uptake_indices);
    
    is_release=ismember(pareja_indice, conjunto_indices);

    %Initialize variables
    count=1;
    C=zeros(My,Mx);
    C_old=C;
    D=Diff;
    C_old(D>0)=Cb; % concentration in mM;
    C_store=zeros(My,Mx,round(nt/rate));
    AUX =(D~= 0);
    %kapa=kapa * dt;
    %Q = (Q /(alpha * dt)); % concentration in mM;
    Q1 = (Q / dt); % To compare with Rohishas data
    U =(Q1/10^4)*6; % uptake rate = 6 molecules / 100 ms
    k_freq=0;
    k_freq2=floor((1/freq)/dt);
    
    kup_freq=3000; % GABA uptake delay = 3 ms
    freq_up=10; %10 Hz
    kup_freq2=kup_freq+floor((1/freq_up)/dt); % 1 event every 0.1 s = 100 ms

    % Time Loop for FT euler method
    for k=1:nt
        if k==k_freq2 + 1
            k_freq=k_freq + floor((1/freq)/dt);
            k_freq2=k_freq2+floor((1/freq)/dt);
        end
        if k==kup_freq2 
            kup_freq=kup_freq + floor((1/freq_up)/dt);
            kup_freq2=kup_freq2+floor((1/freq_up)/dt);
        end
        
        % Apply source and uptake condition
        if k <= floor(k_freq + tp/dt)
            Q= Q1;
        else
            Q=0;
        end
        
        if k == floor(kup_freq)
            uptake=U;
        else
            uptake=0;
        end

        %Boundary conditions
        if BC == 1
            C(:,1)=sc * C_old(:,2);
            C(:,Mx)=sc * C_old(:,Mx-1);
            C(1,:)=sc * C_old(2,:);
            C(My,:)=sc * C_old(My-1,:);
            loopvalues_i = (2:My-1);
            loopvalues_j = (2:Mx-1);
        else
            loopvalues_i = (1:My);
            loopvalues_j = (1:Mx);
        end
        

        % Space loop for CS euler method
        for i=loopvalues_i
            for j=loopvalues_j
                aux=j+Mx*(i-1);
                %Skip if D = 0
                if D(i,j)==0 
                    continue
                end

                % Apply source  and uptake at specific locations 
                
                S = Q * is_release(aux);
                C_up = 0; % No uptake
                %C_up = uptake * is_uptake(aux);
                if BC == 0
                    % Calculate Cy  
                    if i == 1
                        Dmean1 = ((D(i+1,j) + D(i,j)) / 2) * (D(i+1,j) ~= 0);
                        Dmean2 = ((D(My,j) + D(i,j)) / 2 )* (D(My,j) ~= 0);
                        Cy = (Dmean1*(C_old(i + 1, j)-C_old(i, j)) + Dmean2*(C_old(My, j)-C_old(i, j))) / dy2;

                    elseif i == My
                        Dmean4 = ((D(1,j) + D(i,j)) / 2 )*  (D(1,j) ~= 0);
                        Dmean3 = ((D(i-1,j) + D(i,j)) / 2) * (D(i-1,j) ~= 0);
                        Cy = (Dmean4*(C_old(1, j)-C_old(i, j)) + Dmean3*(C_old(i - 1, j)-C_old(i, j))) / dy2;

                    else
                        Dmean1 = ((D(i+1,j) + D(i,j)) / 2) * (D(i+1,j) ~= 0);
                        Dmean3 = ((D(i-1,j) + D(i,j)) / 2) * (D(i-1,j) ~= 0);                   
                        Cy = (Dmean1*(C_old(i + 1, j)-C_old(i, j)) + Dmean3*(C_old(i - 1, j)-C_old(i, j))) / dy2;
                    end

                    % Calculate Cx 

                    if j == 1
                        Dmean5 = ((D(i,Mx) + D(i,j)) / 2 ) * (D(i,Mx) ~= 0);
                        Dmean6 = ((D(i,j+1) + D(i,j)) / 2) * (D(i,j+1) ~= 0);
                        Cx = (Dmean5*(C_old(i, Mx)- C_old(i, j)) + Dmean6*(C_old(i, j + 1) - C_old(i, j))) / dx2;

                    elseif j == Mx
                        Dmean8 = ((D(i,1) + D(i,j)) / 2) * (D(i,1) ~= 0);
                        Dmean7 = ((D(i,j-1) + D(i,j)) / 2) * (D(i,j-1) ~= 0);
                        Cx = (Dmean7*(C_old(i, j - 1)- C_old(i, j)) + Dmean8*(C_old(i, 1) - C_old(i, j))) / dx2;

                    else
                        Dmean6 = ((D(i,j+1) + D(i,j)) / 2) * (D(i,j+1) ~= 0);
                        Dmean7 = ((D(i,j-1) + D(i,j)) / 2) * (D(i,j-1) ~= 0);
                        Cx = (Dmean7*(C_old(i, j - 1)- C_old(i, j)) + Dmean6*(C_old(i, j + 1) - C_old(i, j))) / dx2;
                    end
                else
                    Dmean1 = ((D(i+1,j) + D(i,j)) / 2) *(D(i+1,j) ~= 0);
                    Dmean3 = ((D(i-1,j) + D(i,j)) / 2) * (D(i-1,j) ~= 0);                   
                    Cy = (Dmean1*(C_old(i + 1, j)-C_old(i, j)) + Dmean3*(C_old(i - 1, j)-C_old(i, j))) / dy2;
                    Dmean6 = ((D(i,j+1) + D(i,j)) / 2) * (D(i,j+1) ~= 0);
                    Dmean7 = ((D(i,j-1) + D(i,j)) / 2) * (D(i,j-1) ~= 0);
                    Cx = (Dmean7*(C_old(i, j - 1)- C_old(i, j)) + Dmean6*(C_old(i, j + 1) - C_old(i, j))) / dx2;
                end
                
                % Calculate new concentration
                C(i,j)=(Cx+Cy)*dt+S*dt+(1-kapa*Mask(i,j))*C_old(i, j)-C_up*dt;

                if C(i,j)<0
                    C(i,j)=0;
                end

            end
        end
        
        %Update C_old for the next time step
        C_old=C;
        
        %Save data if it's the first iteration or rate-th iteration
        if k==1 || rem(k,rate)==0
            C_store(:,:,count)= C(1:My,1:Mx);
            count=count+1
        end
        
    end
    toc
end