function [P,Label]=DiffAnalysis(op,file,molec,model,ds,M,S0,U,UT,RTI,radius,ax,NumAnalysis,fit)

    [color, MOLEC] = colorcodefunction (model, molec);
    
    %Define axes units and labels 
    xlabstring=['Time [' UT ']'];
    ylabstring = [molec ' [' U ']'];
    str=[model ' Simulation'];
    p=gobjects(1,7);
    % Load Data and convert the units from SI to desired ones.
    scaleC = convertFromSIunit(U);
    scaleT = convertFromSIunit(UT);
    C=load(file);
    Y=C.Data .* scaleC; 
    t=C.t.*scaleC; 
    T=C.t(end) * scaleT; 

    % Get Measuring Points
    [MP,angles]=MeasuringPoints(radius,S0,M,ds,op);

    % Extract data for valid measuring points
    for i = 1:size(MP,1)
        if strcmp(model, '3D') == 1
            Y1(i,:) = squeeze(Y(MP(i, 1), MP(i, 2), MP(i, 3), :));
        else
            Y1(i,:) = squeeze(Y(MP(i, 1), MP(i, 2), :));
        end
    end
    
    % Find indices where the maximum value of Y1 is greater than 0
    valid_indices = max(Y1, [], 2) > 0;
    MPnew = MP(valid_indices, :);
    anglesnew = angles(valid_indices);
    anglesnew = rad2deg(anglesnew);
    for i = 1:size(MPnew,1)
        if strcmp(model,'3D')==1
            Ytotal(:,i) = squeeze(Y(MPnew(i,1),MPnew(i,2),MPnew(i,3),:));
            y=squeeze(Y(MPnew(i,1),MPnew(i,2),MPnew(i,3),:));
            [ymax,idx]=max(y);
            Ytotal_max(:,i)=ymax;
            Tmax(:,i)=t(idx);
            Imax(:,i)=idx;
        else
            
            Ytotal(:,i) = Y(MPnew(i,1),MPnew(i,2),:);
            y=squeeze(Y(MPnew(i,1),MPnew(i,2),:));
            [ymax,idx]=max(y);
            Ytotal_max(:,i)=ymax;
            Tmax(:,i)=t(idx);
            Imax(:,i)=idx;

        end
        if fit ==1
            x=t';
%             [~,index] = min(abs(x-50));
            [fitresult, gof] = createFit_decay(x(idx:end), Ytotal(idx:end,i));
            R(i)=gof.rsquare;
        end
    end
    
   
    % Save data
    if fit == 1
        save('diffcurve.mat','MPnew','MP','t','Ytotal','anglesnew','Ytotal_max','Tmax','Imax');
    else
        save('diffcurve.mat','MPnew','MP','t','Ytotal','anglesnew','Ytotal_max','Tmax','Imax');
    end

    % OP == 2 display diffusion curves
    if op==2

        R=num2str(radius*10^6);
        titlestr1={[molec,' diffusion at ',R,' µm']};
        %Plot data as mean +-std
        p(1)=plot(ax,t',Ytotal(:,1),'-','LineWidth',2,'Color',color);
        hold(ax,'on')
        plot(ax,t,Ytotal);
        %hold(ax,'on');
        %jbfill(t',Ymin,Ymax,color,color,0,0.3);
        %Indicate when the source is on
        % a=area(ax,[0,50],[10,10],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        % a.FaceAlpha=0.3;
        hold(ax,'on')
    end

    % Compare Models or molecules
    if op==3 || op==4
        if NumAnalysis~=1
            color=[0.6,0.6,0.6];
        end
        hold(ax,'on')

        %p(1)=jbfill(t,Ymax',Ymin',color,color,0,0.6);
        p(1)=plot(ax,t',Ytotal(:,1),'-','LineWidth',2,'Color',color);
        hold(ax,'on')
        plot(ax,t,Ytotal,'-','LineWidth',2,'Color',color);
        set(ax,'FontSize',14);
        if op==3
            titlestr1={[molec,' diffusion']};
        elseif op==4
            str=[molec ' Simulation'];
            titlestr1={[model,' diffusion']};
        end 
    end
    ax.XLim=[0 T];
    P=[p(1)];
    Label={str};
    
    % Display RTI experimental data
    Titles = {'Agar','Brain','Brain','Agar', 'Brain', 'Brain + Uptake'};
    LabelsRTI456 = {'TMA (Nicholson 2001)', 'dex70K (Nicholson 2001)', 'dex3K (Nicholson 2001)', 'BSA (Nicholson 2001)'};
    LabelsRTI123 = {'TMA (Hrabetová 2016)', 'TMA (Hrabetová 2016)', 'TMA (Kserr 2014)','_',  };
    fileNames123 = {'AgaroseData.xlsx', 'BrainData.xlsx','highestC.xlsx','lowestC.xlsx'};
    fileNames456 = {'Agar.txt', 'dex70K.txt', 'dex3K.txt', 'BSA.txt'};
    count=0;
    
    for i=1:6
        if RTI(i)==1
            titlestr1 = Titles{i};
            if RTI(i)== 1 || RTI(i)== 2
                C=xlsread(fileNames123{i});
                C(:,1)=C(:,1)-10;
                RTIlabel={str,LabelsRTI123{i}};
                ax.XLim=[0 140];
                ax.YLim=[0 1];
                idx=i+1;
            elseif RTI(i)== 3
                C1 = xlsread(fileNames123{i});
                C2 = xlsread(fileNames123{i+1});
                fun=fit(C2(:,1),C2(:,2),'smoothingspline');
                Cfit2=fun(C1(:,1));
                C(:,1)=C1(:,1);
                C(:,2)=C1(:,2);
                C(:,3)=Cfit2;
                RTIlabel={str,LabelsRTI123{i}};
                ax.XLim=[0 200];
                ax.YLim=[0 2];
                idx=i+1;
            elseif RTI(i)== 4 || RTI(i)== 5 || RTI(i)== 6
                fileID = fopen(fileNames456{MOLEC},'r');
                sizeC = [2 Inf];
                C=fscanf(fileID,'%f %f',sizeC)';
                fclose(fileID);
                RTIlabel={str,LabelsRTI456{MOLEC}};
                ax.XLim=[0 400];
                ax.YLim=[0 8];
                idx=4+MOLEC;
            end
            
            count=count+1;
            hold(ax,'on')
            p(idx)=plot(ax,C(:,1),C(:,2:end),'--','LineWidth',1.5,'Color',[0.2,0.2,0.2]);
            if count==3
                P=[p(1),p(2),p(3)];
                Label={str,'Agar','Brain'};
            elseif count == 2
                P=[p(1),p(idx)];
                Label={str,RTIlabel};
            end
        end 
    end
 
    % Set title
    ax.Title.String = titlestr1;
    ax.Title.FontSize = 12;
    ax.XLabel.String = xlabstring;
    ax.YLabel.String = ylabstring;
end