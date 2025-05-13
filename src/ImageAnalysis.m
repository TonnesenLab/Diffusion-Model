function [alpha,lamda]=ImageAnalysis(app)
% ImageAnalysis - This function gives an estimate of the volume fraction
%               and tortuosity of a given image. 
% Input: app - The app. content
% Output: alpha - volume fraction 
%         lamda - tortuosity

    Name=[app.imfolder,app.ImageDropDown_4.Value];
    C=imread(Name);
    J=imcomplement(C);
    Max_pvalue=255;% double(max(J(:)))
    imshow(J)
    AUX=imcrop;
    ECS=AUX<90;
    struct=AUX>=90;
    count_ECS=sum(ECS,'all');
    count_struct=sum(struct,'all');
    alpha=(count_ECS/count_struct)*100;

    imshow(C)
    AUX2=imcrop;
    lamda=sqrt(Max_pvalue/mean(AUX2(:)));

end