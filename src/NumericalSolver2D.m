function [C_store]=NumericalSolver2D(M,ds2,dt,nt,Diff,S0,Q,kapa,tp,rate,sc,Mask,Cb, BC)
% NumericalSolver2D - This function solves in 2D and Psuedo3D the diffusion
%                      equation using forward euler method.
% Input: M - Image size
%        ds2 - pixel size squared                 
%        dt - time step for the iterations
%        nt - number of iterations
%        Diff - Diffusion matrix where the free diffusion coefficient is scaled 
%               by the diffusion probability for each pixel.
%        S0 - release site
%        Q - magnitude of the source 
%        kapa - Clearance due to uptake or diffusion in z
%        tp - duration of the release pulse
%        rate - iteration save rate. The concentration is save after the number of iterations indicated by the user.
%        sc - scape coefficient at the edges of the image to account for
%        Mask - mask sushi image
%        Cb - basal concentraion
%        BC - binary value that indicates toroidal (0) or scape (1) boundary
%               conditions
% Output: C_store - Matrix of the concentration of each pixel at each saved
%                   instant of time

 
    %Extract source coordinates
    Sy=S0(1);
    Sx=S0(2);
    Mx=M(2);
    My=M(1);
    dy2=ds2(1);
    dx2=ds2(2);

    %Initialize variables
    count=1;
    C=zeros(My,Mx);
    C_old=C;
    D=Diff;
 
    C_old(D>0)=Cb; % concentration in mM;
    C_store=zeros(My,Mx,round(nt/rate));
    %kapa=kapa * dt;
    %Q = (Q /(alpha * dt)); % concentration in mM;
    Q = (Q / dt); % To compare with Rohishas data
    tic
    
    % Time Loop for FT euler method
    for k=1:nt
        
        % Apply source condition
        if k <= floor(tp/dt)
            Q = Q; 
        else
            Q=0;
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
                %Skip if D = 0
                if D(i,j)==0 
                    continue
                end

                % Apply source condition
                S = Q*(Sy == i && Sx == j);

                % Calculate Cy  
                if i == 1
                    Dmean1 = (D(i+1,j) + D(i,j)) / 2 * (D(i+1,j) ~= 0);
                    Dmean2 = (D(My,j) + D(i,j)) / 2 * (D(My,j) ~= 0);
                    Cy = (Dmean1*(C_old(i + 1, j)-C_old(i, j)) + Dmean2*(C_old(My, j)-C_old(i, j))) / dy2;
                    
                elseif i == My
                    Dmean4 = (D(1,j) + D(i,j)) / 2 * (D(1,j) ~= 0);
                    Dmean3 = (D(i-1,j) + D(i,j)) / 2 * (D(i-1,j) ~= 0);
                    Cy = (Dmean4*(C_old(1, j)-C_old(i, j)) + Dmean3*(C_old(i - 1, j)-C_old(i, j))) / dy2;
                    
                else
                    Dmean1 = (D(i+1,j) + D(i,j)) / 2 * (D(i+1,j) ~= 0);
                    Dmean3 = (D(i-1,j) + D(i,j)) / 2 * (D(i-1,j) ~= 0);                   
                    Cy = (Dmean1*(C_old(i + 1, j)-C_old(i, j)) + Dmean3*(C_old(i - 1, j)-C_old(i, j))) / dy2;
                end
                
                % Calculate Cx 
                
                if j == 1
                    Dmean5 = (D(i,Mx) + D(i,j)) / 2 * (D(i,Mx) ~= 0);
                    Dmean6 = (D(i,j+1) + D(i,j)) / 2 * (D(i,j+1) ~= 0);
                    Cx = (Dmean5*(C_old(i, Mx)- C_old(i, j)) + Dmean6*(C_old(i, j + 1) - C_old(i, j))) / dx2;
                    
                elseif j == Mx
                    Dmean8 = (D(i,1) + D(i,j)) / 2 * (D(i,1) ~= 0);
                    Dmean7 = (D(i,j-1) + D(i,j)) / 2 * (D(i,j-1) ~= 0);
                    Cx = (Dmean7*(C_old(i, j - 1)- C_old(i, j)) + Dmean8*(C_old(i, 1) - C_old(i, j))) / dx2;
                    
                else
                    Dmean6 = (D(i,j+1) + D(i,j)) / 2 * (D(i,j+1) ~= 0);
                    Dmean7 = (D(i,j-1) + D(i,j)) / 2 * (D(i,j-1) ~= 0);
                    Cx = (Dmean7*(C_old(i, j - 1)- C_old(i, j)) + Dmean6*(C_old(i, j + 1) - C_old(i, j))) / dx2;
                end
                
                % Calculate new concentration
                C(i,j)=(Cx+Cy)*dt+S*dt+(1-kapa*Mask(i,j))*C_old(i, j);
                % C(i,j)=(Cx+Cy)*dt+S*dt+(1-kapa)*C_old(i, j);

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