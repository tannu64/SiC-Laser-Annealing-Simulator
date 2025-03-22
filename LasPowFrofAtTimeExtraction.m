function [LasNunber,LasPowProf]=LasPowFrofAtTimeExtraction(LasPro,tMesh)
%extract the laser power profile in time line
%LasNunber:the beam number of applied laser
%LasPowProf:laser power profile, unit:W/cm2

%cite function to get x related arrays
[tPointNumber,tAxis,tDelta]=tAxisExtraction(tMesh);
M=size(LasPro,1); %get the beam number of applied laser
LasNunber=M;   %number of laser beam
N=tPointNumber;  %number of x point
LasPowProf=zeros(M,N);  %creat laser power profile matrix with original values of 0
%fill the laser power profile matrix
for i=1:M
    %to tell if the profile is required as Gaussian
    if LasPro(i,5)==1
        tao=LasPro(i,2)/2.355;
        for j=1:N
            LasPowProf(i,j)=10^6*LasPro(i,3)*exp(-(tAxis(j)-3*tao-LasPro(i,4))^2/(2*(tao)^2))/(tao*(2*pi)^0.5);  
        end
    %to tell if the profile is required as trapezoid   
    elseif LasPro(i,5)==2
        %the trapezoid profile is defined with 3us incress and decress
        %process
        if LasPro(i,2)<3
            error("trapezoid laser pulsewidth must be over 3us")
        end
        for j=1:N
            if tAxis(j)<LasPro(i,4)
                LasPowProf(i,j)=0;
            elseif tAxis(j)<(LasPro(i,4)+3)
                LasPowProf(i,j)= (tAxis(j)-LasPro(i,4))/3*(10^6*LasPro(i,3)/LasPro(i,2));
            elseif tAxis(j)<(LasPro(i,4)+LasPro(i,2))
                LasPowProf(i,j)=10^6*LasPro(i,3)/LasPro(i,2);
            elseif tAxis(j)<(LasPro(i,4)+LasPro(i,2)+3)
                LasPowProf(i,j)= (LasPro(i,4)+LasPro(i,2)+3-tAxis(j))/3*(10^6*LasPro(i,3)/LasPro(i,2));
            else
                LasPowProf(i,j)=0;
            end
        end
    else
        error("LasPro's last colume defined wrong, it is only 1 or 2")
    end
end
end