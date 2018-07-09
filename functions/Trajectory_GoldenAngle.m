function [Traj,DensityComp]=Trajectory_GoldenAngle(ntviews,nX);
%Generate 2D golden-angle radial trajectory and density copensation

% a=111.246117975;

%Gn = (1 + sqrt(5))/2;
%a=180/Gn;
a = 23.6281;

radian=mod((0:a:(ntviews-1)*a)*pi/180,2*pi);
Rho=[-floor(nX/2):floor(nX/2)];
Rho=Rho(1:nX);
Rho=Rho+0.5;
for ii=1:size(Rho,2)
    for jj=1:size(radian,2)
        Rho_temp=Rho/(nX/2)*0.5;
        X(ii,jj)=-Rho_temp(ii)*sin(radian(1,jj));
        Y(ii,jj)=Rho_temp(ii)*cos(radian(1,jj));
    end
end

Traj=X+i*Y;
DensityComp=sqrt(X.^2+Y.^2);
DensityComp=DensityComp/(max(DensityComp(:)));
return