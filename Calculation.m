
function [TMax,Tatx]=Calculation(MatPro,LasPro,RefAbs,xMesh,tMesh,ProfilAtDepth)
%extract the Temperature profile
%TMax:the array of highest temperature profile along depth(x)
%Tatx: the matrix of temperature profile along time at various depth(x)

%cite the function block
[Nx,Px,Dx,DxTrs]=xAxisExtraction(xMesh);
[Nt,Pt,Dt]=tAxisExtraction(tMesh);
[Cap,Den,Cond,LPR,a,R]=MatLasProtertyExtraction(MatPro,LasPro,RefAbs,xMesh);
[NL,PL]=LasPowFrofAtTimeExtraction(LasPro,tMesh);
%to tell whether the Dt is set smally enough
KMax=max(Cond);
CMin=min(Cap);
DMin=min(Den);
xMin=min(DxTrs);
r=KMax*Dt/(CMin*DMin*xMin^2);
if r>=0.5
    error("t mesh step is too big!! please compress it");
end
%creat profile matrix to save the temperature data at various depth
ND=size(ProfilAtDepth,2);
Tatx=zeros(ND,Nt);
%to tell whether the ProfileatDepth exceeds the x space
if min(ProfilAtDepth)<0||max(ProfilAtDepth)>max(Px)
    error("ProfilAtDepth exceeds the x space defined by xmesh, please correct it");
end
%creat array to save the hightest temperature along depth
TMax=zeros(1,Nx);  
%creat arrays to save the temperature profile of former moment and current  
T1=27*ones(1,Nx);    
T2=27*ones(1,Nx);
%initialize the rate of progress
CompletePercentage=0; 
for i=1:Nt   %start a loop on the timeline
    %show the rate of progress
    if round(i/Nt*100)~=CompletePercentage   %if the rate of progress changes, output new rate
        CompletePercentage=round(i/Nt*100);       
        fprintf("Complete Percentage: %d %%\n",CompletePercentage);
    end
    % calculate the temperature change of first x point 
    % calculate the powerdensity by thermalconductivity
    DP2=Cond(1)*(T1(2)-T1(1))/(DxTrs(1)*Dx(1)/10^4);   
    % calculate the powerdensity by laser irradiation
    DP3=0;
    for j=1:NL    %start a loop on the multiple laser
        DP3=DP3+PL(j,i)*(1-R(j))*a(j,1)*LPR(j,1);
    end
    % calculate the change of temperature 
    DP1=DP2+DP3;
    T2(1)=T1(1)+DP1*Dt/(Cap(1)*Den(1));
    % calculate the temperature change of inner x points
    for k=2:Nx-1  %x position loop
        DP2=Cond(k)*(T1(k+1)-T1(k))/(DxTrs(k)*Dx(k)/10^4)-Cond(k-1)*(T1(k)-T1(k-1))/(DxTrs(k-1)*Dx(k)/10^4);
        DP3=0;
        for j=1:NL
            DP3=DP3+PL(j,i)*(1-R(j))*a(j,k)*LPR(j,k);
        end
        DP1=DP2+DP3;
        T2(k)=T1(k)+DP1*Dt/(Cap(k)*Den(k));
    end
    % calculate the temperature change of last x point
    DP2=-Cond(Nx-1)*(T1(Nx)-T1(Nx-1))/(DxTrs(Nx-1)*Dx(Nx)/10^4);
    DP3=0;
    for j=1:NL  %laser loop
        DP3=DP3+PL(j,i)*(1-R(j))*a(j,Nx)*LPR(j,Nx);
    end
    DP1=DP2+DP3;
    T2(Nx)=T1(Nx)+DP1*Dt/(Cap(Nx)*Den(Nx));
    T1=T2;  
    %record the data we are going to stroe 
    NDC=1;
    for k=1:Nx
        if TMax(k)<T2(k)    %save the highest temperature profile
            TMax(k)=T2(k);  
        end
        if (-Dx(k)/2)<=(ProfilAtDepth(NDC)-Px(k))&&(ProfilAtDepth(NDC)-Px(k))<(Dx(k)/2)   %save the temperature profile at vaious depth
           Tatx(NDC,i)=T2(k);
           if NDC<ND
                NDC=NDC+1;   
           end
        end     
    end
end
end
