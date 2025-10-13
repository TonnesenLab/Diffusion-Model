function [C_store]=NumericalSolver3D(M,ds2,dt,nt,Diff,S0,Q,kapa,tp,T,rate,f,alpha,sc,Mask,Cb, BC)
    %Extract source coordinates
    Sy=S0(1);
    Sx=S0(2);
    Sz=S0(3);
    Mx=M(2);
    My=M(1);
    Mz=M(3);
    dy2=ds2(1);
    dx2=ds2(2);
    dz2=ds2(3);
    ds=sqrt(ds2);
    
    sc_z=sc*(ds2(1))/(ds(1)*ds(3));

    %Initialize variables
    count=1;

    AUX=zeros(My,Mx,Mz);
    AUX(Diff>0)=1;
    C_old=Cb * AUX;
    C_store=zeros(My,Mx,Mz,round(nt/rate));
    C=zeros(My,Mx,Mz);
    D=zeros(My,Mx);

    %Q = Q /(alpha * dt);
    Q=Q/dt;
    kapa=kapa*dt; % only for TMA simulation. 
    

    tic
 % Time Loop for FT euler method
    for k=1:nt 
        
        % Source On and Off after some time.
        if k <= floor(tp/dt)
            Q = Q;
        else
            Q = 0;
        end
        %Boundary conditions          
        %if BC == 1
            C(:, :, 1) = sc_z .* C_old(:, :, 2);
            C(:, :, Mz) = sc_z .* C_old(:, :, Mz - 1);
            C(:, 1, :) = sc .* C_old(:, 2, :);
            C(:, Mx, :) = sc .* C_old(:, Mx - 1, :);
            C(1, :, :) = sc .* C_old(2, :, :);
            C(My, :, :) = sc .* C_old(My - 1, :, :);
            loopvalues_i = (2:My-1);
            loopvalues_j = (2:Mx-1);
%         else
%             C(:, :, 1) = sc .* C_old(:, :, 2);
%             C(:, :, Mz) = sc .* C_old(:, :, Mz - 1);
%             loopvalues_i = (1:My);
%             loopvalues_j = (1:Mx);
%         end
    
        %Space loop for CS euler method
        for z = 2:Mz-1
            D = Diff(:,:,z);
            Ddown = Diff(:,:,z-1);
            Dup = Diff(:,:,z+1);
            Cold = C_old(:,:,z);
            Coldz_down = C_old(:,:,z-1);
            Coldz_up = C_old(:,:,z+1);
        
            for i = loopvalues_i
                for j = loopvalues_j

                    if D(i, j) == 0
                        continue
                    end

                    % Source On and Off depending on location
                    
                    S = Q * (Sz == z && Sy == i && Sx == j);
                    
                    Coldz = Cold(i,j);

                    % Calculate Cy 
                    Dmean1 = (D(i+1,j) + D(i,j)) / 2 * (D(i+1,j) ~= 0);
                    Dmean3 = (D(i-1,j) + D(i,j)) / 2 * (D(i-1,j) ~= 0);
                    
%                     if i == 1
%                         Dmean2 = (D(My,j) + D(i,j)) / 2 * (D(My,j) ~= 0);
%                         Cy = (Dmean1*(Cold(i + 1, j)-Coldz) + Dmean2*(Cold(My, j)-Coldz)) / dy2;
%                         
%                     elseif i == My
%                         Dmean4 = (D(1,j) + D(i,j)) / 2 * (D(1,j) ~= 0);
%                         Cy = (Dmean4*(Cold(1, j)-Coldz) + Dmean3*(Cold(i - 1, j)-Coldz)) / dy2;
%                         
%                     else
                        Cy = (Dmean1*(Cold(i + 1, j)-Coldz) + Dmean3*(Cold(i - 1, j)-Coldz)) / dy2;
%                     end

                    %Calculate Cx
                    Dmean6 = (D(i,j+1) + D(i,j)) / 2 * (D(i,j+1) ~= 0);
                    Dmean7 = (D(i,j-1) + D(i,j)) / 2 * (D(i,j-1) ~= 0);
                    
%                     if j == 1
%                         Dmean5 = (D(i,Mx) + D(i,j)) / 2 * (D(i,Mx) ~= 0);
%                         Cx = (Dmean5*(Cold(i, Mx)- Coldz) + Dmean6*(Cold(i, j + 1) - Coldz)) / dx2;
%                         
%                     elseif j == Mx
%                         Dmean8 = (D(i,1) + D(i,j)) / 2 * (D(i,1) ~= 0);
%                         Cx = (Dmean7*(Cold(i, j - 1)- Coldz) + Dmean8*(Cold(i, 1) - Coldz)) / dx2;
%                         
%                     else
                        Cx = (Dmean7*(Cold(i, j - 1)- Coldz) + Dmean6*(Cold(i, j + 1) - Coldz)) / dx2;
%                     end

                    %Calculate Cz
                    Dmean9 = (Ddown(i,j) + D(i,j)) / 2 * (Ddown(i,j) ~= 0);
                    Dmean10 = (Dup(i,j) + D(i,j)) / 2 * (Dup(i,j) ~= 0);
                    
                    Cz = (Dmean9*(Coldz_down(i,j)- Coldz) + Dmean10*(Coldz_up(i,j) - Coldz)) / dz2;

                    % Update concentration
                    C(i,j,z)=(Cx+Cy+Cz)*dt+S*dt+(1 - kapa*Mask(i,j,z))*Coldz;

                end
            end                       
        end

        % Update C_old for the next time step
        C_old = C;

        %Save the Data
        if k==1 || rem(k,rate)==0
            C_store(:,:,:,count)= C(1:My,1:Mx,1:Mz);
            count=count+1

%                 dlg=uiprogressdlg(f);
%                 dlg.ShowPercentage='on';
%                 dlg.Value=k/nt;
        end
        
    end
toc
end