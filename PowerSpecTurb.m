function [TurbPhase] = PowerSpecTurb( Base,D, r0)
% ko=2*pi/Lambda;
Lo=1000; % Outer scale 
lo=0.001; % inner scale 
N = Base.H;
delx=1./(Base.H.*Base.dx);
kl= 5.92/lo/(2*pi);
kL= 1/Lo;

%generate spatial frequencies
fx=(-N/2: N/2-1)*delx;
[kx, ky]=meshgrid(fx);
% radial coordinate of wavenumbers 
knorm = sqrt(kx.^2 +ky.^2);

%Phase spectrum
PhaseSpec =  0.023 * r0^(-5/3) * exp (-(knorm./ kl).^2 ) ./(knorm.^2 + kL^2).^(11/6); %0 .23*(r0).^(-5/3)*knorm.^(-11/3);% Phase spectrum 
PhaseSpec(N/2 +1, N/2 +1) = 0;
C = ( randn(N) + 1i.*randn(N) );

%turbulent phase
TurbPhase = ifft2 (ifftshift(sqrt (PhaseSpec).*C*(delx))).*( Base.H ).^2 .*real(exp(-1i*pi*(Base.X + Base.Y))); %% ifftshift(ifft2(ifftshift(sqrt (PhaseSpec).*C*(delx))))*( Base.H ).^2 .*real(exp(1i*pi*(Base.X + Base.Y)));
phz_lo = zeros(size(TurbPhase)); 

%Add subharmonics to tubulent phase
NP=20;
for p = 1 : NP
    del_f = 1 / (NP^p*D);
    fx = (-1 : 1 ) * del_f;
    [fx, fy] = meshgrid(fx);
    f=sqrt(fx.^2+fy.^2);
    fm = 5.92/lo/(2*pi);
    f0 = 1/Lo;
    PSD_phi = 0.023*r0^(-5/3) * exp (-(f/fm).^2)./ (f.^2 + f0^2).^(11/6);
    PSD_phi(2,2) = 0;
    cn = (randn(3) + 1i*randn(3)) .* sqrt(PSD_phi)*del_f;
    SH = zeros(N);
    for ii = 1:9
        SH = SH + cn(ii)*exp(1i*2*pi*(fx(ii)*Base.X+fy(ii)*Base.Y));
    end
    phz_lo = phz_lo + SH;
end
TurbPhase = exp(1i * real( TurbPhase)); % final screen as a wavefront 
