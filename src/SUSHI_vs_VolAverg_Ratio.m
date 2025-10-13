%Data 1 SUSHI
Data11=Data1.*10^3; % Convert Data to uM
Data11=Data11+0.001;
for i=1:length(Data11)
   aux=Data11(:,:,i);
   aux(Sushi400_400==0)=0;
   Data11(:,:,i)=aux;
end

%Data 2 Avrg. Volume
Data12=Data2.*10^3; % Convert Data to uM
Data22=Data12+0.001;
for i=1:length(Data22)
   aux2=Data22(:,:,i);
   aux2(Sushi400_400==0)=0;
   Data22(:,:,i)=aux2;
end

%Ratio = Sushi/V.Av

Data=Data11./Data22;
Data(isnan(Data))=0;
aux3=Data(:,:,3);

save('Sushi_VA_ratio_100Hz.mat','Data','t')