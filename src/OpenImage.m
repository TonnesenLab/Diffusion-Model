
function [Mask0,zstack]=OpenImage(NameMask,aux)

    Info=imfinfo(NameMask);
    NumberImages=length(Info);
    zstack=NumberImages > 1;

    if zstack
        C=zeros(Info(1).Height,Info(1).Width ,NumberImages);
        for i=1:NumberImages
            A=imread(NameMask, 'Index', i);
            C(:,:,i)=A;
        end 
        MaxPixelVal= max(C(:));
    else
        C=imread(NameMask);
        MaxPixelVal= max(C(:));
    end
    Mask0 = (double(C) .* aux) ./ double(MaxPixelVal);
    
end