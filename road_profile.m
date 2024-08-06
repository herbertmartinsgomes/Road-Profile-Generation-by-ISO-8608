%GERAÇÃO DO PERFIL DE PISTA
function [Ub]=road_profile(vel,ww,no0,Gn0,nls,tdef1,tdef2,T,np,wmin,wmax,ntire,iflag)
    iflag=1;      %prints or not graphs
    s=rng;rng('default');
    w=linspace(wmin,wmax,nls)'; deltaw=(wmax-wmin)/nls;
    ctf=0.01;%    %Cut-off frequency (1/m) for correlation between parallel road tracks
    deltat=T/np;z=1*sqrt(-1);
    nmin=wmin/(2*pi*vel);nmax=wmax/(2*pi*vel);fmax=vel*nmax;
    fs=1/deltat;
    Ub=zeros(ntire,np);Gw=zeros(1,nls);
    %pontosextras=floor(0/deltat);
    pontosextras=floor(tdef2/deltat);
    f=(linspace(wmin/(2*pi),wmax/(2*pi),nls))';
    trilha=zeros((np+pontosextras),ntire);
    for k=1:nls
        no=w(k)/(2*pi*vel);         %número de onda
        Gn=Gn0*((no/no0))^(-ww);    %valor de Gn   
        Gw(k)=Gn*(1.0/(2.0*pi*vel));
        gamma=((ctf^2)/(ctf^2+(f(k)/vel).^2));
        Gwc(k)=gamma*Gw(k);
           
        Gwt(k)=Gw(k)*exp(z*w(k)*tdef1);
        Gwct(k)=gamma.*Gwt(k); 
        Gwf(k)=Gw(k)*exp(z*w(k)*tdef2);
        Gwcf(k)=gamma.*Gwf(k); 
    end
        phase=2.0*pi*rand(nls,2);%phi1=0;
    t=linspace(0,T,np+pontosextras); %tempo ficticio para levar em conta a defasagem entre pneus dianteiros e traseiros
    
    n=np+pontosextras;
    trilha(:,1)=sum(repmat(sqrt(abs(Gw')*deltaw),1,n).*sin((w*t)+repmat(phase(:,1),1,n)),1)';
    trilha(:,2)=sum(repmat(sqrt(abs((Gw-Gwc)')*deltaw),1,n).*sin((w*t)+repmat(phase(:,1),1,n)),1)' + sum( repmat(sqrt(abs(Gwc')*deltaw),1,n).*sin((w*t)+repmat(phase(:,2),1,n)),1)';
    
    Ub(1,:)=trilha(pontosextras+1:pontosextras+np,1);   %Trilha do pneu dianteiro esquerdo
    Ub(2,:)=trilha(pontosextras+1:pontosextras+np,2);   %Trilha do pneu dianteiro direito
    
    Ub(3,:)=trilha(pontosextras-round((tdef1/tdef2)*pontosextras)+1:pontosextras+np-round((tdef1/tdef2)*pontosextras),1);        %Trilha do pneu intermediario esquerdo
    Ub(4,:)=trilha(pontosextras-round((tdef1/tdef2)*pontosextras)+1:pontosextras+np-round((tdef1/tdef2)*pontosextras),2);        %Trilha do pneu intermediario direito
    
    Ub(5,:)=trilha(1:np,1);   %Trilha do pneu traseiro  esquerdo 
    Ub(6,:)=trilha(1:np,2);   %Trilha do pneu traseiro  direito

%     no=w./(2*pi*vel);         %número de onda
%     Gn=Gn0*(no./no0).^(-ww);    %valor de Gn   
%     Gw=Gn.*(1.0/(2.0*pi*vel));        
%     t=linspace(0,T,np+pontosextras); %tempo ficticio para levar em conta a defasagem entre pneus dianteiros e traseiros    
%     n=np+pontosextras;    
%     phase=repmat(2.0*pi*rand(nls,1),1,n);    
%     trilha(:,1)=(sqrt(abs(Gw)*deltaw)'*sin((w*t)+phase))';  
%     Ub(1,:)=trilha(pontosextras+1:pontosextras+np,1);   %Trilha do pneu 
%     
    if iflag==1        
        %Densidade Espectral unicaudal pelo periodograma
        df=1/(t(np));
        nfft=round(np/2+1);window=round(0.1*length(Ub(1,:)));noverlap=round(0.1*window);
        [hh, ff]=pwelch(Ub(1,:),window,noverlap,nfft,fs,'onesided');
        hh=(2*hh.*df*(np/2+1)/nfft);
        nesp2=ff./vel; hh=(vel/df)*(hh);
        nesp=linspace(nmin,nmax,np/2+1); nesp=nesp';      %número de onda conforme a velocidade do veículo
        Gn=Gn0*(nesp.*(1/no0)).^(-ww);
        i=find(nesp2>=nmin,1,'first');j=find(nesp2>=nmax,1,'first');
        nespn=logspace(log(0.001),10,2000);   %número de onda conforme a velocidade do veículo
        nespn=nespn';
        GnA=2*(4^(1+1))*10^-6*(nespn.*(1/no0)).^(-ww);
        GnB=2*(4^(2+1))*10^-6*(nespn.*(1/no0)).^(-ww);
        GnC=2*(4^(3+1))*10^-6*(nespn.*(1/no0)).^(-ww);
        GnD=2*(4^(4+1))*10^-6*(nespn.*(1/no0)).^(-ww);
        GnE=2*(4^(5+1))*10^-6*(nespn.*(1/no0)).^(-ww);
        GnF=2*(4^(6+1))*10^-6*(nespn.*(1/no0)).^(-ww);
        GnG=2*(4^(7+1))*10^-6*(nespn.*(1/no0)).^(-ww);
        GnH=2*(4^(8+1))*10^-6*(nespn.*(1/no0)).^(-ww);
    
        loglog(nespn,GnA,'k','LineWidth',2);hold on;
        loglog(nespn,GnB,'k','LineWidth',2);
        loglog(nespn,GnC,'k','LineWidth',2);
        loglog(nespn,GnD,'k','LineWidth',2);
        loglog(nespn,GnE,'k','LineWidth',2);
        loglog(nespn,GnF,'k','LineWidth',2);
        loglog(nespn,GnG,'k','LineWidth',2);   
        loglog(nesp2(i:j),hh(i:j),'b','LineWidth',2);
    
        title('Road Class');xlabel('Spatial Frequency \it{n}\rm [cycles/m]');ylabel('\it{G_n}\rm(\it{n}\rm) Displacement Power Spectral Density  [m³]');  
        grid on;% grid minor;
        ylim([1E-8 1.1]);xlim([0.0049 10]);
        text(0.6,5e-7/1.1,'\bf{A}');
        text(0.6,20e-7/1.1,'\bf{B}');
        text(0.6,80e-7/1.1,'\bf{C}');
        text(0.6,320e-7/1.1,'\bf{D}');
        text(0.6,1280e-7/1.1,'\bf{E}');
        text(0.6,5120e-7/1.1,'\bf{F}');
        text(0.6,20480e-7/1.1,'\bf{G}');
        text(0.6,81920e-7/1.1,'\bf{H}');
        hold on;
        set(gca, 'XColor', [0.05 0.05 0.05]);
        %----
        rng(s);
    end
end
