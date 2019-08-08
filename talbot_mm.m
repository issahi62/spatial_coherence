function [It,x]=talbot_mm(z)
%
% Field behind binary amplitude grating at distance z
%
% z= propagation distance
% d = period of the grating
% f = fill factor
% Dx = size of the calculation window 
% N= number of the sampling points
% M = number of included diffraction orders
% lambda= wavelength
It =[];
lambda=633e-9;
N=4000;
M=237;
Dx=0.65e-3; %calculted by given parameter
f=0.1;
d=125*10^-6;

%z=[0.0494,0.0247,0.0165,0.01234,0.0082,0.0049]

k=2*pi/lambda;
x=linspace(-Dx/2,Dx/2,N);
Et=zeros(size(x));

for m=-M:M
% Add here the expressions for the complex amplitude coefficients
% tm and t0 that can be analytically solved from Eq. (3).

tm=1i/(2*pi*m)*(exp(-1i*2*pi*m*f)-1);
    if m==0
        tm=f;
    end
kxm=2*pi*m/d;
kzm=sqrt(k^2-kxm.^2);
Et=Et+tm*exp(1i*(kxm*x+kzm*z));
It=[(It), abs(Et).^2];
end

figure; 
plot(length(It), It.*10^2)
xlabel('x [mm]')
ylabel('Intensity (a.u)')