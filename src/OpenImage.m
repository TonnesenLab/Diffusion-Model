function [Mask0,zstack]=OpenImage(NameMask,aux)
% OpenImage - This function reads the chosen sushi mask into a matrix and
%               transforms it into a 0-1 scale matrix
% Input: NameMask - Name and path of sushi mask
%        aux - scale value for homogeniezed field of view

% Output: Mask0 - Sushi Mask from 0-1 
%         zstack - binary value indicate if the image is a z-stack (1) or just
%               a single plane (0). 

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