%% Gaussian beam reflection from a thin film (diff. angles of incidence).
N=2^10; %range of numbers
L=1000; %electromagnetic field is defined in the interval (um)
a=5; %normalization of the illuminated area on the axis x (um)
r0b=3; %radius of the Gaussian beam (um)
lambda=0.6328; %wavelength (um)
hi=20*lambda; %optical path length, 20-film thickness
k0=2*pi/lambda; %wave vector at a given wavelength
gam=lambda/L;
al=20; %angle of incidence
z=636.62*0.0001; %number of diffraction lengths (dimensionless)

%% >Normalization--------------------------------------------------------%%
r0=r0b/a; %radius of the Gaussian beam (dimensionless)
Lb=L/a; %dimensionless width of the interval
DeltaX=Lb/N; %distance between points on the x axis (dimensionless)
DeltaKx=1/(N*DeltaX); %distance between the points in reciprocal space
p=-N/2:1:N/2-1; %points on the x axis
xx=p*DeltaX;
xxx=xx.*a;
m=-N/2:1:N/2-1; %points in reciprocal space

%% >Calculation of the reflection coefficients---------------------------%%
                            %permittivities
e1=1;      %air
e2=3.118;  %sapphire         
e3=13.138; %silicon
                            %refractive indices
n1=sqrt(e1); 
n2=sqrt(e2); 
n3=sqrt(e3); 
gg=(m*gam).^2;
cosQ1= (1-gg).^0.5;
cosQ2=(1-(e1/e2).*gg).^0.5;
cosQ3=(1-(e1/e3).*gg).^0.5;
h=hi.*n2; %optical thickness of the film
b=(2.*pi.*h./lambda).*cosQ2; %phase parameter

                %reflection coefficient for s-polarization (perpendicular)
% reflection coefficient on the border 1 - 2
r12s=(n1.*cosQ1-n2.*cosQ2)./(n1.*cosQ1+n2.*cosQ2);
% reflection coefficient on the border 2 - 3
r23s=(n2.*cosQ2-n3.*cosQ3)./(n2.*cosQ2+n3.*cosQ3); 
% amplitude reflection coefficient
Rs=(r12s+r23s.*exp(2.*1i.*b))./(1+r12s.*r23s.*exp(2.*1i.*b)); 

                %reflection coefficient for p-polarized (parallel)
% reflection coefficient on the border 1 - 2
r12p=(n2.*cosQ1-n1.*cosQ2)./(n2.*cosQ1+n1.*cosQ2);
% reflection coefficient on the border 2 - 3
r23p=(n3.*cosQ2-n2.*cosQ3)./(n3.*cosQ2+n2.*cosQ3); 
% amplitude reflection coefficient
Rp=(r12p+r23p.*exp(2.*1i.*b))./(1+r12p.*r23p.*exp(2.*1i.*b)); 
Q=acosd(cosQ1); %(h*180/pi)

%% >Analytical formula of the Gaussian beam------------------------------%%
ka2=(a*(2*pi)/lambda)^2;
AGZ=(1+4*z^2)^(-1/2)*exp(-xx.^2/(1+4*z^2));
FAG=z*ka2+atan(2*z)+xx.^2*(2*z/(4*z^2+1));
EGZ=AGZ.*exp(1i*FAG);

%% >Field of Gaussian beam with the added phase--------------------------%%
nnn=sin(al/180*pi)*L/lambda;
cosalpha0=nnn*lambda/L;
sintheta10=cosalpha0; 
FA0=2*pi*a/lambda*sintheta10.*xx; %slope of the Gaussian beam
A0r=exp(-(xx).^2./(r0).^2); %amplitude of the Gaussian beam at the start

%settings of polarization
HHi=0+1i;
E0=A0r.*exp(1i.*(FA0))/(2^0.5); 
E01=HHi.*A0r.*exp(1i.*(FA0))/(2^0.5); 

%on 1-st plane
Phaza=angle(E0); %phase of the incident beam
Int=abs(E0).^2; %intensity of the incident beam
Ampl=abs(E0); %amplitude of the incident beam
E00=trapz(Int); %integrated intensity of the incident beam
E001=trapz(Ampl); %integrated amplitude of the incident beam

%% >The expansion in the Fourier spectrum--------------------------------%%
F0=fftshift(fft(fftshift(E0)));        
F01=fftshift(fft(fftshift(E01)));  

%% >Field in front of the substrate--------------------------------------%%
zz=-0.5*z.*(m./Lb).^2; %phase of plane wave components
                %Fourier spectrum of the incident beam on the 2-nd plane
FZ=F0.*exp(1i*k0*a*(((k0.^2*a.^2)-(4*pi.^2*(m./Lb).^2)).^0.5)*z);
FZ1=F01.*exp(1i*k0*a*(((k0.^2*a.^2)-(4*pi.^2*(m./Lb).^2)).^0.5)*z);
                %Field of the incident beam on the 2-nd plane
EZ=fftshift(ifft(fftshift(FZ)));  
EZ1=fftshift(ifft(fftshift(FZ1))); 
Phaza1=angle(EZ); %phase of the incident beam on the 2-nd plane
Int1=abs(EZ).^2; %intensity of the incident beam on the 2-nd plane
Ampl1=abs(EZ); %amplitude of the incident beam on the 2-nd plane
EZZ=trapz(Int1); %integrated intens. of the incident beam on the 2-nd pl.

%% >Field passed into the substrate--------------------------------------%%
EZ1=fftshift(ifft2(fftshift(FZ1))); %field on the 2-nd plane (transmitted)
koef=(n3.*cosQ3)./(n1.*cosQ1);
Phaza4=angle(EZ1); %phase of the transmitted beam on the 2-nd plane
Int4=koef.*abs(EZ1).^2; %intensity of the transmitted beam on the 2-nd pl.
Ampl4=abs(EZ1); %аmplitude of the transmitted beam on the 2-nd plane
EZZ1=trapz(Int4); %integrated intens. of the transm. beam on the 2-nd pl.

%% >Field after reflection from the substrate----------------------------%%
                            %s-polarization
FZs=FZ.*Rs; %Fourier spectrum of the 2-nd plane s-polariztion (reflected)
EZs=fftshift(ifft(fftshift(FZs))); %field on the 2-nd plane (reflected)
Phaza2=angle(EZs); %phase of the reflected beam on the 2-nd plane
Int2=abs(EZs).^2; %intensity of the incident beam on the 2-nd plane
Ampl2=abs(EZs); %аmplitude of the incident beam on the 2-nd plane
EZZs=trapz(Int2); %integrated intens. of the incident beam on the 2-nd pl.
                            %p-polarization
FZp=FZ1.*Rp; %Fourier spectrum of the 2-nd plane p-polariztion (reflected)
EZp=fftshift(ifft(fftshift(FZp))); %field on the 2-nd plane (reflected)
Phaza3=angle(EZp); %phase of the reflected beam on the 2-nd plane
Int3=abs(EZp).^2; %intensity of the incident beam on the 2-nd plane
Ampl3=abs(EZp); %аmplitude of the incident beam on the 2-nd plane
EZZp=trapz(Int3); %integrated intens. of the incident beam on the 2-nd pl.
Int4=Int2+Int3; %Intensity s&p-polarization

%% >Angles---------------------------------------------------------------%%
Hi= EZs./EZp;
Az=-angle(Hi).*1/2;
Az=Az.*180/pi;
El=atand((abs(Hi)-1)./(abs(Hi)+1));
Amp=(abs(EZp).^2+abs(EZs).^2).^0.5;

%% >Visualization--------------------------------------------------------%%
%{
subplot(6,1,1), plot(xxx,Int1,'LineWidth',1);
title('Intensity in front of the substrate:');
subplot(6,1,2), plot(xxx,Int2,'LineWidth',1);
title('Intensity after reflection from the substrate (s-polarization):');
subplot(6,1,3), plot(xxx,Int3,'LineWidth',1);
title('Intensity after reflection from the substrate (p-polarization):');
subplot(6,1,4), plot(xxx,Int4,'LineWidth',1);
title('Intensity after reflection from the substrate (s&p-polarization):');
subplot(6,1,5), plot(m,abs(F0),'LineWidth',1);
title('The expansion in the Fourier spectrum:');
subplot(6,1,6), plot(Q,abs(Rs).^2,Q,abs(Rp).^2,'LineWidth',1);
title('Amplitude reflection coefficient (rs,rp):');
%}
%{
figure
subplot(5,1,1), plot(xxx,Az);
title('The azimuth of the polarization ellipse:');
subplot(5,1,2), plot(xxx,El);
title('Angle of Ellipticity of the polarization ellipse:');
subplot(5,1,3), plot(xxx,Amp);
title('The amplitude of the polarization ellipse:');
subplot(5,1,4), plot(xxx,(Amp.^2./(El.^2+1)).^0.5);
title('Semi-major axis of the polarization ellipse:');
subplot(5,1,5), plot(xxx,((Amp.^2).*(El.^2)./(El.^2+1)).^0.5);
title('Minor axis of the polarization ellipse:');
%}
subplot(6,1,1), plot(xxx,Int1,'LineWidth',1.5);
title('Интенсивность перед подложкой:');
subplot(6,1,2), plot(xxx,Int2,'LineWidth',1.5);
title('Интенсивность после отражения от подложки(s-поляризация):');
subplot(6,1,3), plot(xxx,Int3,'LineWidth',1.5);
title('Интенсивность после отражения от подложки(p-поляризация):');
subplot(6,1,4), plot(xxx,Int4,'LineWidth',1.5);
title('Интенсивность после отражения от подложки(s&p-поляризация):');
subplot(6,1,5), plot(m,abs(F0),'LineWidth',1.5);
title('Разложение в Фурье спектр:');
subplot(6,1,6), plot(Q,abs(Rs).^2,Q,abs(Rp).^2,'LineWidth',1.5);
title('Амплитуда коээфициентов отражения (rs,rp):');
figure
plot(xxx,Az,'LineWidth',3,'Color','red');
xlabel('X');
ylabel('Градус');
title('Азимут эллипса поляризации:');
grid on
figure
plot(xxx,El,'LineWidth',3,'Color','green');
xlabel('X');
ylabel('Градус');
title('Угол эллиптичности поляризации:');
grid on

