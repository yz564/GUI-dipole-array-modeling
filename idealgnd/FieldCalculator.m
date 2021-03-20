function FieldCalculator(freq,xstart,ystart,height,xlength,ylength,xstep,ystep,Magnitude,Phase)
eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;
n0=376.8194;
c = 3e8;
k0 = 2*pi*freq*(mu0*eps0)^0.5;
unit=1e-3;
Dipoles=importdata('ReconstructedDipole.csv');
Dipoles.textdata(1,:)=[];
Type_d=Dipoles.textdata(:,1);
[DipoleNum,nnn]=size(Dipoles.data);

X=xstart:xstep:xstart+xlength;
Y=ystart:ystep:ystart+ylength;
numx=size(X,2);
numy=size(Y,2);
N=numx*numy;

x=zeros(1,N);
y=zeros(1,N);
z=zeros(1,N);

for i=1:numx
    for j=1:numy
        x(j+(i-1)*numy) = X(i);
        y(j+(i-1)*numy) = Y(j);
        z(j+(i-1)*numy) = height;
    end
end

Dx =Dipoles.data(:,1)*unit;
Dy =Dipoles.data(:,2)*unit;
Dz =Dipoles.data(:,3)*unit;
%Magnitude =Dipoles.data(:,4);
%Phase =Dipoles.data(:,5);
% D=Magnitude.*exp(1i*(Phase/180));
Pz=zeros(DipoleNum,1);
Px=zeros(DipoleNum,1);
Py=zeros(DipoleNum,1);
Mz=zeros(DipoleNum,1);
Mx=zeros(DipoleNum,1);
My=zeros(DipoleNum,1);

Ez=zeros(DipoleNum,N);
Ex=zeros(DipoleNum,N);
Ey=zeros(DipoleNum,N);
Hz=zeros(DipoleNum,N);
Hx=zeros(DipoleNum,N);
Hy=zeros(DipoleNum,N);

%Pz
for n=1:DipoleNum
    dx=Dx(n)*ones(1,N);
    dy=Dy(n)*ones(1,N);
    dz=Dz(n)*ones(1,N);
    r1=((x-dx).^2+(y-dy).^2+(z-dz).^2).^0.5;
    fr1=exp(-1j*k0.*r1)./r1;
    g1r1=(3./((k0.*r1).^2)+1j.*3./(k0.*r1)-1).*fr1;
    g2r1=(2./((k0.*r1).^2)+1j.*2./(k0.*r1)).*fr1;
    g3r1=(1./(k0.*r1)+1j).*fr1;
    KE=-1j*k0*n0/4/pi;
    KH=k0/4/pi;
    
    if strcmpi(Type_d(n),'Px')==1
        Px(n)= Magnitude(n)*exp(1i*(Phase(n)/180*pi));
        Ex(n,:)=KE*(-((y-dy).^2+(z-dz).^2)./r1./r1.*g1r1+g2r1)*Px(n);
        Ey(n,:)=KE*((x-dx).*(y-dy)./r1./r1.*g1r1)*Px(n);
        Ez(n,:)=KE*((z-dz).*(x-dx)./r1./r1.*g1r1)*Px(n);
        Hx(n,:)=zeros(N, 1)*Px(n);
        Hy(n,:)=KH*(-(z-dz)./r1.*g3r1)*Px(n);
        Hz(n,:)=KH*((y-dy)./r1.*g3r1)*Px(n);
    elseif strcmpi(Type_d(n),'Py')==1
        Py(n)= Magnitude(n)*exp(1i*(Phase(n)/180*pi));
        Ex(n,:)=KE*((x-dx).*(y-dy)./r1./r1.*g1r1)*Py(n);
        Ey(n,:)=KE*(-((z-dz).^2+(x-dx).^2)./r1./r1.*g1r1+g2r1)*Py(n);
        Ez(n,:)=KE*((y-dy).*(z-dz)./r1./r1.*g1r1)*Py(n);
        Hx(n,:)=KH*((z-dz)./r1.*g3r1)*Py(n);
        Hy(n,:)=zeros(N, 1)*Py(n);
        Hz(n,:)=KH*(-(x-dx)./r1.*g3r1)*Py(n);
    elseif strcmpi(Type_d(n),'Pz')==1
        Pz(n)= Magnitude(n)*exp(1i*(Phase(n)/180*pi));
        Ex(n,:)=KE*((z-dz).*(x-dx)./r1./r1.*g1r1)*Pz(n);
        Ey(n,:)=KE*((y-dy).*(z-dz)./r1./r1.*g1r1)*Pz(n);
        Ez(n,:)=KE*((-((y-dy).^2+(x-dx).^2))./r1./r1.*g1r1+g2r1)*Pz(n);
        Hx(n,:)=KH*(-(y-dy)./r1.*g3r1)*Pz(n);
        Hy(n,:)=KH*((x-dx)./r1.*g3r1)*Pz(n); 
        Hz(n,:)=zeros(N, 1)*Pz(n);
    elseif strcmpi(Type_d(n),'Mx')==1
        Mx(n)= Magnitude(n)*exp(1i*(Phase(n)/180*pi));
        Ex(n,:)=zeros(N, 1)*Mx(n);
        Ey(n,:)=KE*(-(z-dz)./r1.*g3r1).*k0*Mx(n);
        Ez(n,:)=KE*((y-dy)./r1.*g3r1).*k0*Mx(n);
        Hx(n,:)=KH*((-((y-dy).^2+(z-dz).^2))./r1./r1.*g1r1+g2r1).*k0*Mx(n);
        Hy(n,:)=KH*((x-dx).*(y-dy)./r1./r1.*g1r1).*k0*Mx(n);
        Hz(n,:)=KH*((z-dz).*(x-dx)./r1./r1.*g1r1).*k0*Mx(n);
    elseif strcmpi(Type_d(n),'My')==1
        My(n)= Magnitude(n)*exp(1i*(Phase(n)/180*pi));
        Ex(n,:)=KE*((z-dz)./r1.*g3r1).*k0*My(n);  
        Ey(n,:)=zeros(N, 1)*My(n);
        Ez(n,:)=KE*(-(x-dx)./r1.*g3r1).*k0*My(n);
        Hx(n,:)=KH*((x-dx).*(y-dy)./r1./r1.*g1r1).*k0*My(n);
        Hy(n,:)=KH*(-((x-dx).^2+(z-dz).^2)./r1./r1.*g1r1+g2r1).*k0*My(n);  
        Hz(n,:)=KH*((y-dy).*(z-dz)./r1./r1.*g1r1).*k0*My(n);  
    elseif strcmpi(Type_d(n),'Mz')==1 
        Mz(n)= Magnitude(n)*exp(1i*(Phase(n)/180*pi));
        Ex(n,:)=KE*(-(y-dy)./r1.*g3r1).*k0*Mz(n);
        Ey(n,:)=KE*((x-dx)./r1.*g3r1).*k0*Mz(n);
        Ez(n,:)=zeros(N, 1)*Mz(n);
        Hx(n,:)=KH*((z-dz).*(x-dx)./r1./r1.*g1r1).*k0*Mz(n);
        Hy(n,:)=KH*((y-dy).*(z-dz)./r1./r1.*g1r1).*k0*Mz(n);
        Hz(n,:)=KH*(-((x-dx).^2+(y-dy).^2)./r1./r1.*g1r1+g2r1).*k0*Mz(n);
    end
end
Hx0=sum(Hx);
Hy0=sum(Hy);
Hz0=sum(Hz);
Ex0=sum(Ex);
Ey0=sum(Ey);
Ez0=sum(Ez);

%%
    Ex_mag=abs(Ex0);
    Ex_phase=180*angle(Ex0)/pi;
    Ex_Mag=zeros(numy,numx);
    Ex_Phase=zeros(numy,numx);
    
    for i=1:numy
        for j=1:numx
            Ex_Mag(numy-j+1,i)=Ex_mag(j+numx*(i-1));
            Ex_Phase(numy-j+1,i)=Ex_phase(j+numx*(i-1));
        end
    end
    figure(1)
    subplot(2, 3, 1);
    imagesc(Ex_Mag);
    title('Ex Mag', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 3, 4);
    imagesc(Ex_Phase);
    title('Ex Phase', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet')
    %%
    Ey_mag=abs(Ey0);
    Ey_phase=180*angle(Ey0)/pi;
    Ey_Mag=zeros(numy,numx);
    Ey_Phase=zeros(numy,numx);
    
    for i=1:numy
        for j=1:numx
            Ey_Mag(numy-j+1,i)=Ey_mag(j+numx*(i-1));
            Ey_Phase(numy-j+1,i)=Ey_phase(j+numx*(i-1));
        end
    end
    figure(1)
    subplot(2, 3, 2);
    imagesc(Ey_Mag);
    title('Ey Mag', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 3, 5);
    imagesc(Ey_Phase);
    title('Ey Phase', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet')
    %%
    Ez_mag=abs(Ez0);
    Ez_phase=180*angle(Ez0)/pi;
    Ez_Mag=zeros(numy,numx);
    Ez_Phase=zeros(numy,numx);
    
    for i=1:numy
        for j=1:numx
            Ez_Mag(numy-j+1,i)=Ez_mag(j+numx*(i-1));
            Ez_Phase(numy-j+1,i)=Ez_phase(j+numx*(i-1));
        end
    end
    figure(1)
    subplot(2, 3, 3);
    imagesc(Ez_Mag);
    title('Ez Mag', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 3, 6);
    imagesc(Ez_Phase);
    title('Ez Phase', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet')
    %%
    Hx_mag=abs(Hx0);
    Hx_phase=180*angle(Hx0)/pi;
    Hx_Mag=zeros(numy,numx);
    Hx_Phase=zeros(numy,numx);
    
    for i=1:numy
        for j=1:numx
            Hx_Mag(numy-j+1,i)=Hx_mag(j+numx*(i-1));
            Hx_Phase(numy-j+1,i)=Hx_phase(j+numx*(i-1));
        end
    end
    figure(2)
    subplot(2, 3, 1);
    imagesc(Hx_Mag);
    title('Hx Mag', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 3, 4);
    imagesc(Hx_Phase);
    title('Hx Phase', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet')
    %%
    Hy_mag=abs(Hy0);
    Hy_phase=180*angle(Hy0)/pi;
    Hy_Mag=zeros(numy,numx);
    Hy_Phase=zeros(numy,numx);
    
    for i=1:numy
        for j=1:numx
            Hy_Mag(numy-j+1,i)=Hy_mag(j+numx*(i-1));
            Hy_Phase(numy-j+1,i)=Hy_phase(j+numx*(i-1));
        end
    end
    figure(2)
    subplot(2, 3, 2);
    imagesc(Hy_Mag);
    title('Hy Mag', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 3, 5);
    imagesc(Hy_Phase);
    title('Hy Phase', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet')
%         %%
    Hz_mag=abs(Hz0);
    Hz_phase=180*angle(Hz0)/pi;
    Hz_Mag=zeros(numy,numx);
    Hz_Phase=zeros(numy,numx);
    
    for i=1:numy
        for j=1:numx
            Hz_Mag(numy-j+1,i)=Hz_mag(j+numx*(i-1));
            Hz_Phase(numy-j+1,i)=Hz_phase(j+numx*(i-1));
        end
    end
    figure(2)
    subplot(2, 3, 3);
    imagesc(Hz_Mag);
    title('Hz Mag', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 3, 6);
    imagesc(Hz_Phase);
    title('Hz Phase', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet')
end