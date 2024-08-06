%% Truck model(cabin, seat and boot) 13 Dof + biodynamic model with 12 DoF
tic;close all;clear all;
ntire=6;           %NUmber of tires [int]
vel=80/3.6;        %Speed [m/s]
Lt=1000;           %Track length [m]
time=Lt/vel;       %Time of analysis [s]
flag_print=1;      %Flag to print graphs         


b1=5.18;           %Dist. from front susp. to vehicle's CG[m]
b2=0.62;           %Dist. from inter. susp. to vehicle's CG [m]
b3=1.97;           %Dist. from rear susp. to vehicle's CG [m]
b4=6.78;           %Dist. from front cabin support to vehicle's CG [m]
b5=4.68;           %Dist. from rear cabin support to vehicle's CG [m]
Lx1=b1;            %Dist. left front susp. to vehicle's CG [m]
Lx2=Lx1;           %Dist. rigth front susp. to vehicle's CG [m]
Lx3=-b2;           %Dist. left interm. susp. to vehicle's CG [m]
Lx4=Lx3;           %Dist. rigth interm. susp. to vehicle's CG [m]
Lx5=-b3;           %Dist. left rear susp. to vehicle's CG [m]
Lx6=Lx5;           %Dist. rigth rear susp. to vehicle's CG [m]
wb_intermediario=Lx1-Lx4;   %wheelbase for interm. susp. [m]  => front susp. to interm. susp. 1 to 2
wb=Lx1-Lx6;        %wheelbase [m] => front susp. to rear suspension  1 to 3

%Referênce to tires
%Ub(1,:)=Track for front left tire   %Ub(2,:)=Track for front rigth tire
%Ub(3,:)=Track for left interm. tire %Ub(4,:)=%Track for rigth interm. tire
%Ub(5,:)=%Track for rear left tire   %Ub(6,:)=%Track for rear rigth tire

%%Road profile generation according to ISO 8608
T=Lt/vel;            %Simulation time
dt1=0.0005;           %User specified time interval for analysis
nmin=0.011;          %minimum spactial frequency [ciclo/m]
nmax=2.83;           %maximum spatial frequency [ciclo/m]
wmin=2*pi*vel*nmin;  %minimum angular frequency [rad/s]
wmax=2*pi*vel*nmax;  %maximum angular frequency [rad/s]
fmin=vel*nmin;       %minimum frequency [Hz]
fmax=vel*nmax;       %maximum frequency [Hz]
dt0=(1/(10*fmax));   %minimum integration time interval to represent accordingly the track [s]
dt=min(dt0,dt1);     %actual time interval for analysis to be possible to generate road waviness
np=round(T/dt);      %number of time steps for analysis
t=linspace(0.0,T,np);%vector of time instants of the analysis

%Generation of Road profile by ISO 8608:2016
tdef1=(Lx1-Lx3)/vel;%time interval from first tire axis to the intermediate tire axis
tdef2=(Lx1-Lx5)/vel;%time interval from intermediate tire axis and the last tire axis
n0=0.1;             %Spatial frequency by ISO8608 [cycle/m]
no0=n0;
road=['A','B','C','D','E','F','G','H'];
cn=1;                %rad profile class; 1=A, 2=B, 3=C, 4=D, 5=E; (1=smooth, 5=very irregular)
Gn0=(4^(cn+1))*10^-6;%Wave number by the ISO standard according to the road uneaviness [m^3]
ww=2.0;              %road profile exponent for the spectrald density
nls=1000;             %number of spectral lines use dfor irregularities generation
w=linspace(wmin,wmax,nls); %discrete circular frequencies [rad/s]
nwave=linspace(nmin,nmax,nls); %discrete spatial frequencies  [1/m]
deltaw=(wmax-wmin)/nls;deltaf=(fmax-fmin)/nls;deltan=(nmax-nmin)/nls; %frequency resolution 
f=linspace(fmin,fmax,nls); %discrete frequencies [Hz]
fs=1/dt;df=1/T;dn=df/vel;  %Sampling frequencies
%Initialize road profiles and corresponding derivatives
ub1=zeros(1,np);ub2=ub1;ub3=ub1;ub4=ub1;ub5=ub1;ub6=ub1;
ubp1=ub1;ubp2=ub1;ubp3=ub1;ubp4=ub1;ubp5=ub1;ubp6=ub1;

%Generates the road profile in all tires, along time, according to
%velocity, correlation,...
[Ub]=road_profile(vel,ww,no0,Gn0,nls,tdef1,tdef2,T,np,wmin,wmax,ntire,flag_print);
%Update the road profile adding any previous irregularities, e.g., 'bump' 
ub1=ub1+Ub(1,:);ub2=ub2+Ub(2,:);ub3=ub3+Ub(3,:);ub4=ub4+Ub(4,:);ub5=ub5+Ub(5,:);ub6=ub6+Ub(6,:);
%Derivatives of the road profile
ubp1=[0.0 diff(ub1)/dt]; ubp2=[0.0 diff(ub2)/dt]; ubp3=[0.0 diff(ub3)/dt]; ubp4=[0.0 diff(ub4)/dt];ubp5=[0.0 diff(ub5)/dt]; ubp6=[0.0 diff(ub6)/dt];

%Print road profile
if flag_print == 1
    figure; %Road profile(and derivatives) in time scale
    subplot(2,1,1),plot(t, ub1(:),'-r',t ,ub2(:),'-g',t ,ub3(:),'-b',t, ub4(:),'-k',t, ub5(:),'-.k',t, ub6(:),'-.b');title('Perfil da Pista');           legend('ub1 [m]','ub2 [m]','ub3 [m]','ub4 [m]','ub5 [m]','ub6 [m]');    xlabel('Time [s]');ylabel('[m]');
    subplot(2,1,2),plot(t,ubp1(:),'-r',t,ubp2(:),'-g',t,ubp3(:),'-b',t,ubp4(:),'-k',t,ubp5(:),'-.k',t,ubp6(:),'-.b');title('Derivada do Perfil da Pista');legend('ubp1 [m]','ubp2 [m]','ubp3 [m]','ubp4 [m]','ubp5 [m]','ubp6 [m]');xlabel('Time [s]');ylabel('[m]');
    figure;%Road profile (and derivatives) in space scale
    subplot(2,1,1),plot(t*vel, ub1(:),'-r',t*vel ,ub2(:),'-g',t*vel ,ub3(:),'-b',t*vel, ub4(:),'-k',t*vel, ub5(:),'-.k',t*vel, ub6(:),'-.b');title('Perfil da Pista');           legend('ub1 [m]','ub2 [m]','ub3 [m]','ub4 [m]','ub5 [m]','ub6 [m]');    xlabel('Distância [m]');ylabel('[m]');
    subplot(2,1,2),plot(t*vel,ubp1(:),'-r',t*vel,ubp2(:),'-g',t*vel,ubp3(:),'-b',t*vel,ubp4(:),'-k',t*vel,ubp5(:),'-.k',t*vel,ubp6(:),'-.b');title('Derivada do Perfil da Pista');legend('ubp1 [m]','ubp2 [m]','ubp3 [m]','ubp4 [m]','ubp5 [m]','ubp6 [m]');xlabel('Distância [m]');ylabel('[m/m]');
end

