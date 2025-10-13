function [alpha,lamda]=ImageAnalysis(app)
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
alpha
imshow(C)
AUX2=imcrop;
lamda=sqrt(Max_pvalue/mean(AUX2(:)));
lamda

end