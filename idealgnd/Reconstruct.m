function varargout = Reconstruct(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Reconstruct_OpeningFcn, ...
                   'gui_OutputFcn',  @Reconstruct_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function Reconstruct_OpeningFcn(hObject, eventdata, handles, varargin)

clc

handles.output = hObject;

guidata(hObject, handles);

function varargout = Reconstruct_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function editOutput2_Callback(hObject, eventdata, handles)

function editOutput2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editOutput1_Callback(hObject, eventdata, handles)

function editOutput1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox1_Callback(hObject, eventdata, handles)

function checkbox2_Callback(hObject, eventdata, handles)

function checkbox3_Callback(hObject, eventdata, handles)

function checkbox4_Callback(hObject, eventdata, handles)

function checkbox5_Callback(hObject, eventdata, handles)

function checkbox6_Callback(hObject, eventdata, handles)

function editScanCenter_Callback(hObject, eventdata, handles)

function editScanCenter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)

function edit8_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)

function edit9_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit11_Callback(hObject, eventdata, handles)

function edit11_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Model_Callback(hObject, eventdata, handles)

switch get(handles.Model,'Value')
    case 2
        handles.model=1; 
    case 3
        handles.model=2; 
    case 4
        handles.model=3; 
    otherwise
end
handles.output = hObject;
guidata(hObject, handles);

function Model_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Unit_Callback(hObject, eventdata, handles)
switch get(handles.Unit,'Value')
    case 2
        handles.unit=1; 
    case 3
        handles.unit=2; 
    otherwise
end
handles.output = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Unit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function freq_Callback(hObject, eventdata, handles)

function freq_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FINISH_Callback(hObject, eventdata, handles)

set(0,'defaultfigurecolor','w');
unit = 1e-3; %[m]

Dipole=importdata('DipoleInfo.csv');
DipoleLocation_0=unit*Dipole.data;
[DipoleNum,xyz]=size(DipoleLocation_0);
DipoleType=zeros(DipoleNum,6);
Type_d=Dipole.textdata;

for nn=1:DipoleNum
    if strcmpi(Type_d(nn),'Px')==1
        DipoleType(nn,:)=[1,0,0,0,0,0];
    elseif strcmpi(Type_d(nn),'Py')==1
        DipoleType(nn,:)=[0,1,0,0,0,0];
    elseif strcmpi(Type_d(nn),'Pz')==1
        DipoleType(nn,:)=[0,0,1,0,0,0];
    elseif strcmpi(Type_d(nn),'Mx')==1
        DipoleType(nn,:)=[0,0,0,1,0,0];
    elseif strcmpi(Type_d(nn),'My')==1
        DipoleType(nn,:)=[0,0,0,0,1,0];
    elseif strcmpi(Type_d(nn),'Mz')==1
        DipoleType(nn,:)=[0,0,0,0,0,1];
    end
end

% DipoleNum=1;
% DipoleType=[0,0,0,1,0,0]; % [Px, Py, Pz, Mx, My, Mz]
%[dipoleX, dipoleY, dipoleZ]

F1=get(handles.checkbox1,'Value');
F2=get(handles.checkbox2,'Value');
F3=get(handles.checkbox3,'Value');
F4=get(handles.checkbox4,'Value');
F5=get(handles.checkbox5,'Value');
F6=get(handles.checkbox6,'Value');
Unit=handles.unit;
SolutionModel=handles.model;
if SolutionModel==1
    GND=0;
else 
    GND=1;
end

fieldOpt = [F1 F2 F3 F4 F5 F6]; % [useEx, useEy, useEz, useHx, useHy, useHz]

%Scan plane information

ScanCenter = str2num(get(handles.editScanCenter,'String'));
scanCenter = unit*[ScanCenter(1),ScanCenter(2)];
xNumScan = str2num(get(handles.edit6,'String'));    
yNumScan = str2num(get(handles.edit7,'String'));    
xSpcScan = unit*str2num(get(handles.edit8,'String'));
ySpcScan = unit*str2num(get(handles.edit9,'String'));
scanHeight = unit*ScanCenter(3);

eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;
n0=376.8194;
c = 3e8;
freq =1e9*str2num(get(handles.freq,'String')); %[Hz]
k0 = 2*pi*freq*(mu0*eps0)^0.5;
c_m_dipole=1j*n0*2*pi*freq/c;
% Near Field Scanning Plane
scanXYZ=zeros(yNumScan, xNumScan, 3);
cen_x = scanCenter(1);
cen_y = scanCenter(2);

if SolutionModel==3
Data_H1 = textread('H1.fld');
Data_E1 = textread('E1.fld');
Data_H2 = textread('H2.fld');
Data_E2 = textread('E2.fld');
Ex0=Data_E1(:,4)+1i*Data_E1(:,5)-Data_E2(:,4)-1i*Data_E2(:,5);
Ey0=Data_E1(:,6)+1i*Data_E1(:,7)-Data_E2(:,6)-1i*Data_E2(:,7);
Ez0=Data_E1(:,8)+1i*Data_E1(:,9)+Data_E2(:,8)+1i*Data_E2(:,9);
Hx0=Data_H1(:,4)+1i*Data_H1(:,5)+Data_H2(:,4)+1i*Data_H2(:,5);
Hy0=Data_H1(:,6)+1i*Data_H1(:,7)+Data_H2(:,6)+1i*Data_H2(:,7);
Hz0=Data_H1(:,8)+1i*Data_H1(:,9)-Data_H2(:,8)-1i*Data_H2(:,9);
else
Data_H1 = textread('H1.fld');
Data_E1 = textread('E1.fld');
Ex0=Data_E1(:,4)+1i*Data_E1(:,5);
Ey0=Data_E1(:,6)+1i*Data_E1(:,7);
Ez0=Data_E1(:,8)+1i*Data_E1(:,9);
Hx0=Data_H1(:,4)+1i*Data_H1(:,5);
Hy0=Data_H1(:,6)+1i*Data_H1(:,7);
Hz0=Data_H1(:,8)+1i*Data_H1(:,9);
end

if (mod(xNumScan, 2) == 0)    
    left_num_x = xNumScan/2;
    right_num_x = xNumScan/2;
    x_start = cen_x - (left_num_x - 1 + 0.5) * xSpcScan;
    x_end = cen_x + (right_num_x - 1 + 0.5) * xSpcScan; 
    line_x = x_start:xSpcScan:x_end;  
else
    left_num_x = (xNumScan - 1)/2;
    right_num_x = (xNumScan - 1)/2;
    x_start = cen_x - left_num_x * xSpcScan;
    x_end = cen_x + right_num_x * xSpcScan;
    line_x = x_start:xSpcScan:x_end;
end

if (mod(yNumScan, 2) == 0)
    left_num_y = yNumScan/2;
    right_num_y = yNumScan/2;
    y_start = cen_y - (left_num_y - 1 + 0.5) * ySpcScan;
    y_end = cen_y + (right_num_y - 1 + 0.5) * ySpcScan;
    line_y = y_end:-ySpcScan:y_start;    
else
    left_num_y = (yNumScan - 1)/2;
    right_num_y = (yNumScan - 1)/2;
    y_start = cen_y - left_num_y * ySpcScan;
    y_end = cen_y + right_num_y * ySpcScan;
    line_y = y_end:-ySpcScan:y_start;
end
    
for county = 1:1:yNumScan
    scanXYZ(county, :, 1) = line_x;
end 
for countx = 1:1:xNumScan
    scanXYZ(:, countx, 2) = line_y;
end
scanXYZ(:,:,3) = scanHeight;

% Input data: fields on the scanning plane
Ex = zeros(yNumScan, xNumScan);
Ey = zeros(yNumScan, xNumScan);
Ez = zeros(yNumScan, xNumScan);
Hz = zeros(yNumScan, xNumScan);
Hx = zeros(yNumScan, xNumScan);
Hy = zeros(yNumScan, xNumScan);

ScanField = [];  
% Data_H = textread('E_NearField_Height4mm.fld');
% Hx0 = Data_H(:,4) + 1j .* Data_H(:,5);
% Hy0 = Data_H(:,6) + 1j .* Data_H(:,7);
% Hz0 = Data_H(:,8) + 1j .* Data_H(:,9);

if fieldOpt(1)==1
%     Data_Ex= importdata('Ex.csv');
%     Ex0 = Data_Ex.data(:,2) + 1j .* Data_Ex.data(:,3);
    Ex_mag=abs(Ex0);
    Ex_phase=180*angle(Ex0)/pi;
    Ex_Mag=zeros(yNumScan,xNumScan);
    Ex_Phase=zeros(yNumScan,xNumScan);
    
    for i=1:yNumScan
        for j=1:xNumScan
            Ex(yNumScan-j+1,i)=Ex0(j+xNumScan*(i-1));
            Ex_Mag(yNumScan-j+1,i)=Ex_mag(j+xNumScan*(i-1));
            Ex_Phase(yNumScan-j+1,i)=Ex_phase(j+xNumScan*(i-1));
        end
    end
    figure(11)
    subplot(2, 1, 1);
    imagesc(Ex_Mag);
    title('Ex Mag (scan)', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 1, 2);
    imagesc(Ex_Phase);
    title('Ex Phase (scan)', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet');
    ScanEx = reshape(Ex,yNumScan*xNumScan, 1);
    ScanField = [ScanField;ScanEx];  
end

if fieldOpt(2)==1
%     Data_Ey= importdata('Ey.csv');
%     Ey0 = Data_Ey.data(:,2) + 1j .* Data_Ey.data(:,3);
    Ey_mag=abs(Ey0);
    Ey_phase=180*angle(Ey0)/pi;
    Ey_Mag=zeros(yNumScan,xNumScan);
    Ey_Phase=zeros(yNumScan,xNumScan);
    
    for i=1:yNumScan
        for j=1:xNumScan
            Ey(yNumScan-j+1,i)=Ey0(j+xNumScan*(i-1));
            Ey_Mag(yNumScan-j+1,i)=Ey_mag(j+xNumScan*(i-1));
            Ey_Phase(yNumScan-j+1,i)=Ey_phase(j+xNumScan*(i-1));
        end
    end
    figure(12)
    subplot(2, 1, 1);
    imagesc(Ey_Mag);
    title('Ey Mag (scan)', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 1, 2);
    imagesc(Ey_Phase);
    title('Ey Phase (scan)', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet');
    ScanEy = reshape(Ey,yNumScan*xNumScan, 1);
    ScanField = [ScanField;ScanEy];  
end

if fieldOpt(3)==1
%     Data_Ez= importdata('Ez.csv');
%     Ez0 = Data_Ez.data(:,2) + 1j .* Data_Ez.data(:,3);
    Ez_mag=abs(Ez0);
    Ez_phase=180*angle(Ez0)/pi;
    Ez_Mag=zeros(yNumScan,xNumScan);
    Ez_Phase=zeros(yNumScan,xNumScan);
    
    for i=1:yNumScan
        for j=1:xNumScan
            Ez(yNumScan-j+1,i)=Ez0(j+xNumScan*(i-1));
            Ez_Mag(yNumScan-j+1,i)=Ez_mag(j+xNumScan*(i-1));
            Ez_Phase(yNumScan-j+1,i)=Ez_phase(j+xNumScan*(i-1));
        end
    end
    figure(13)
    subplot(2, 1, 1);
    imagesc(Ez_Mag);
    title('Ez Mag (scan)', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 1, 2);
    imagesc(Ez_Phase);
    title('Ez Phase (scan)', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet');
    ScanEz = reshape(Ez,yNumScan*xNumScan, 1);
    ScanField = [ScanField;ScanEz];    
end

if fieldOpt(4)==1
%     Data_Hx= importdata('Hx.csv');
%     Hx0 = Data_Hx.data(:,2) + 1j .* Data_Hx.data(:,3);
    Hx_mag=abs(Hx0);
    Hx_phase=180*angle(Hx0)/pi;
    Hx_Mag=zeros(yNumScan,xNumScan);
    Hx_Phase=zeros(yNumScan,xNumScan);
    
    for i=1:yNumScan
        for j=1:xNumScan
            Hx(yNumScan-j+1,i)=Hx0(j+xNumScan*(i-1));
            Hx_Mag(yNumScan-j+1,i)=Hx_mag(j+xNumScan*(i-1));
            Hx_Phase(yNumScan-j+1,i)=Hx_phase(j+xNumScan*(i-1));
        end
    end
 figure(14)
    subplot(2, 2, 1);
    imagesc(Hx_Mag);
    title('Hx Mag [A/m]', 'fontsize', 12);
    set (gca,'YDir','normal');
    xticks([5 10 15 20])
    xticklabels({'-5','0','5','10'})
    xlabel('X axis [mm]')
    yticks([5 10 15 20])
    yticklabels({'-5','0','5','10'})
    ylabel('Y axis [mm]')  
    colorbar;
    subplot(2, 2, 2);
    imagesc(Hx_Phase);
    title('Hx Phase [deg]', 'fontsize', 12);
    set (gca,'YDir','normal');
    xticks([5 10 15 20])
    xticklabels({'-5','0','5','10'})
    xlabel('X axis [mm]')
    yticks([5 10 15 20])
    yticklabels({'-5','0','5','10'})
    ylabel('Y axis [mm]')  
    colorbar;
    colormap('Jet');
    ScanHx = reshape(Hx,yNumScan*xNumScan, 1);
    ScanField = [ScanField;ScanHx];    
end

if fieldOpt(5)==1
%     Data_Hy= importdata('Hy.csv');
%     Hy0 = Data_Hy.data(:,2) + 1j .* Data_Hy.data(:,3);
    Hy_mag=abs(Hy0);
    Hy_phase=180*angle(Hy0)/pi;
    Hy_Mag=zeros(yNumScan,xNumScan);
    Hy_Phase=zeros(yNumScan,xNumScan);
    
    for i=1:yNumScan
        for j=1:xNumScan
            Hy(yNumScan-j+1,i)=Hy0(j+xNumScan*(i-1));
            Hy_Mag(yNumScan-j+1,i)=Hy_mag(j+xNumScan*(i-1));
            Hy_Phase(yNumScan-j+1,i)=Hy_phase(j+xNumScan*(i-1));
        end
    end
    figure(14)
    subplot(2, 2, 3);
    imagesc(Hy_Mag);
    title('Hy Mag [A/m]', 'fontsize', 12);
    set (gca,'YDir','normal');
    xticks([5 10 15 20])
    xticklabels({'-5','0','5','10'})
    xlabel('X axis [mm]')
    yticks([5 10 15 20])
    yticklabels({'-5','0','5','10'})
    ylabel('Y axis [mm]')   
    
    colorbar;
    subplot(2, 2, 4);
    imagesc(Hy_Phase);
    title('Hy Phase [deg]', 'fontsize', 12);
    set (gca,'YDir','normal'); 
    xticks([5 10 15 20])
    xticklabels({'-5','0','5','10'})
    xlabel('X axis [mm]')
    yticks([5 10 15 20])
    yticklabels({'-5','0','5','10'})
    ylabel('Y axis [mm]')   
    colorbar;
    colormap('Jet');
    ScanHy = reshape(Hy,yNumScan*xNumScan, 1);
    ScanField = [ScanField;ScanHy];    
end


if fieldOpt(6)==1
%     Data_Hz= importdata('Hz.csv');
%     Hz0 = Data_Hz.data(:,2) + 1j .* Data_Hz.data(:,3);
    Hz_mag=abs(Hz0);
    Hz_phase=180*angle(Hz0)/pi;
    Hz_Mag=zeros(yNumScan,xNumScan);
    Hz_Phase=zeros(yNumScan,xNumScan);
    for i=1:yNumScan
        for j=1:xNumScan
            Hz(yNumScan-j+1,i)=Hz0(j+xNumScan*(i-1));
            Hz_Mag(yNumScan-j+1,i)=Hz_mag(j+xNumScan*(i-1));
            Hz_Phase(yNumScan-j+1,i)=Hz_phase(j+xNumScan*(i-1));
        end
    end
    figure(16)
    subplot(2, 1, 1);
    imagesc(Hz_Mag);
    title('Hz Mag (scan)', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    subplot(2, 1, 2);
    imagesc(Hz_Phase);
    title('Hz Phase (scan)', 'fontsize', 12);
    set (gca,'YDir','normal');
    colorbar;
    colormap('Jet');
    ScanHz = reshape(Hz,yNumScan*xNumScan, 1);
    ScanField = [ScanField;ScanHz];
end

Constant=ones(DipoleNum,1);


% Scan plane matrix reshape
[plane_ynum,plane_xnum,plane_depth] = size(scanXYZ);
pointsnum = plane_ynum * plane_xnum;
x = reshape(scanXYZ(:,:,1), pointsnum, 1);
y = reshape(scanXYZ(:,:,2), pointsnum, 1);
z = reshape(scanXYZ(:,:,3), pointsnum, 1);

KE=-1j*k0*n0/4/pi;
KH=k0/4/pi;
FieldNum=sum(fieldOpt);
Tmatrix=zeros(FieldNum*pointsnum,DipoleNum);
for i=1:DipoleNum
    m=find(DipoleType(i,:));
    dipoleX=DipoleLocation_0(i,1);
    dipoleY=DipoleLocation_0(i,2);
    dipoleZ=DipoleLocation_0(i,3);
    dx =ones(pointsnum, 1) * dipoleX;
    dy =ones(pointsnum, 1) * dipoleY;
    dz =ones(pointsnum, 1) * dipoleZ;
   
    if GND==1
    r1=((x-dx).^2+(y-dy).^2+(z-dz).^2).^0.5;
    r2=((x-dx).^2+(y-dy).^2+(z+dz).^2).^0.5;
    fr1=exp(-1j*k0.*r1)./r1;
    fr2=exp(-1j*k0.*r2)./r2;
    g1r1=(3./((k0.*r1).^2)+1j.*3./(k0.*r1)-1).*fr1;
    g2r1=(2./((k0.*r1).^2)+1j.*2./(k0.*r1)).*fr1;
    g3r1=(1./(k0.*r1)+1j).*fr1;
    g1r2=(3./((k0.*r2).^2)+1j*3./(k0.*r2)-1).*fr2;
    g2r2=(2./((k0.*r2).^2)+1j*2./(k0.*r2)).*fr2;
    g3r2=(1./(k0.*r2)+1j).*fr2;
    
    T_Ex_Px=KE*(-((y-dy).^2+(z-dz).^2)./r1./r1.*g1r1+g2r1 + ((y-dy).^2+(z+dz).^2)./r2./r2.*g1r2-g2r2);
    T_Ey_Px=KE*((x-dx).*(y-dy)./r1./r1.*g1r1 - (x-dx).*(y-dy)./r2./r2.*g1r2);
    T_Ez_Px=KE*((z-dz).*(x-dx)./r1./r1.*g1r1 - (z+dz).*(x-dx)./r2./r2.*g1r2);
    T_Hx_Px=zeros(pointsnum, 1);
    T_Hy_Px=KH*(-(z-dz)./r1.*g3r1 + (z+dz)./r2.*g3r2);
    T_Hz_Px=KH*((y-dy)./r1.*g3r1 - (y-dy)./r2.*g3r2);

    T_Ex_Py=KE*((x-dx).*(y-dy)./r1./r1.*g1r1 - (x-dx).*(y-dy)./r2./r2.*g1r2);
    T_Ey_Py=KE*(-((z-dz).^2+(x-dx).^2)./r1./r1.*g1r1+g2r1 + ((z+dz).^2+(x-dx).^2)./r2./r2.*g1r2-g2r2);
    T_Ez_Py=KE*((y-dy).*(z-dz)./r1./r1.*g1r1 - (y-dy).*(z+dz)./r2./r2.*g1r2);
    T_Hx_Py=KH*((z-dz)./r1.*g3r1 - (z+dz)./r2.*g3r2);
    T_Hy_Py=zeros(pointsnum, 1);
    T_Hz_Py=KH*(-(x-dx)./r1.*g3r1 + (x-dx)./r2.*g3r2);

    T_Ex_Pz=KE*((z-dz).*(x-dx)./r1./r1.*g1r1+(z+dz).*(x-dx)./r2./r2.*g1r2);
    T_Ey_Pz=KE*((y-dy).*(z-dz)./r1./r1.*g1r1+(y-dy).*(z+dz)./r2./r2.*g1r2);
    T_Ez_Pz=KE*((-((y-dy).^2+(x-dx).^2))./r1./r1.*g1r1+g2r1+(-((y-dy).^2+(x-dx).^2))./r2./r2.*g1r2+g2r2);
    T_Hx_Pz=KH*(-(y-dy)./r1.*g3r1-(y-dy)./r2.*g3r2);
    T_Hy_Pz=KH*((x-dx)./r1.*g3r1+(x-dx)./r2.*g3r2); 
    T_Hz_Pz=zeros(pointsnum, 1);

    T_Ex_Mx=zeros(pointsnum, 1);
    T_Ey_Mx=KE*(-(z-dz)./r1.*g3r1-(z+dz)./r2.*g3r2).*k0;
    T_Ez_Mx=KE*((y-dy)./r1.*g3r1+(y-dy)./r2.*g3r2).*k0;
    T_Hx_Mx=KH*((-((y-dy).^2+(z-dz).^2))./r1./r1.*g1r1+g2r1+(-((y-dy).^2+(z+dz).^2))./r2./r2.*g1r2+g2r2).*k0;
    T_Hy_Mx=KH*((x-dx).*(y-dy)./r1./r1.*g1r1+(x-dx).*(y-dy)./r2./r2.*g1r2).*k0;
    T_Hz_Mx=KH*((z-dz).*(x-dx)./r1./r1.*g1r1+(z+dz).*(x-dx)./r2./r2.*g1r2).*k0;

    T_Ex_My=KE*((z-dz)./r1.*g3r1+(z+dz)./r2.*g3r2).*k0;  
    T_Ey_My=zeros(pointsnum, 1);
    T_Ez_My=KE*(-(x-dx)./r1.*g3r1-(x-dx)./r2.*g3r2).*k0;
    T_Hx_My=KH*((x-dx).*(y-dy)./r1./r1.*g1r1+(x-dx).*(y-dy)./r2./r2.*g1r2).*k0;
    T_Hy_My=KH*(-((x-dx).^2+(z-dz).^2)./r1./r1.*g1r1+g2r1-((x-dx).^2+(z+dz).^2)./r2./r2.*g1r2+g2r2).*k0;  
    T_Hz_My=KH*((y-dy).*(z-dz)./r1./r1.*g1r1+(y-dy).*(z+dz)./r2./r2.*g1r2).*k0;  

    T_Ex_Mz=KE*(-(y-dy)./r1.*g3r1 + (y-dy)./r2.*g3r2).*k0;
    T_Ey_Mz=KE*((x-dx)./r1.*g3r1 - (x-dx)./r2.*g3r2).*k0;
    T_Ez_Mz=zeros(pointsnum, 1);
    T_Hx_Mz=KH*((z-dz).*(x-dx)./r1./r1.*g1r1 - (z+dz).*(x-dx)./r2./r2.*g1r2).*k0;
    T_Hy_Mz=KH*((y-dy).*(z-dz)./r1./r1.*g1r1 - (y-dy).*(z+dz)./r2./r2.*g1r2).*k0;
    T_Hz_Mz=KH*(-((x-dx).^2+(y-dy).^2)./r1./r1.*g1r1+g2r1 + ((x-dx).^2+(y-dy).^2)./r2./r2.*g1r2-g2r2).*k0;
    else
    r1=((x-dx).^2+(y-dy).^2+(z-dz).^2).^0.5;
    fr1=exp(-1j*k0.*r1)./r1;
    g1r1=(3./((k0.*r1).^2)+1j.*3./(k0.*r1)-1).*fr1;
    g2r1=(2./((k0.*r1).^2)+1j.*2./(k0.*r1)).*fr1;
    g3r1=(1./(k0.*r1)+1j).*fr1;
    T_Ex_Px=KE*(-((y-dy).^2+(z-dz).^2)./r1./r1.*g1r1+g2r1);
    T_Ey_Px=KE*((x-dx).*(y-dy)./r1./r1.*g1r1);
    T_Ez_Px=KE*((z-dz).*(x-dx)./r1./r1.*g1r1);
    T_Hx_Px=zeros(pointsnum, 1);
    T_Hy_Px=KH*(-(z-dz)./r1.*g3r1);
    T_Hz_Px=KH*((y-dy)./r1.*g3r1);

    T_Ex_Py=KE*((x-dx).*(y-dy)./r1./r1.*g1r1);
    T_Ey_Py=KE*(-((z-dz).^2+(x-dx).^2)./r1./r1.*g1r1+g2r1);
    T_Ez_Py=KE*((y-dy).*(z-dz)./r1./r1.*g1r1);
    T_Hx_Py=KH*((z-dz)./r1.*g3r1);
    T_Hy_Py=zeros(pointsnum, 1);
    T_Hz_Py=KH*(-(x-dx)./r1.*g3r1);

    T_Ex_Pz=KE*((z-dz).*(x-dx)./r1./r1.*g1r1);
    T_Ey_Pz=KE*((y-dy).*(z-dz)./r1./r1.*g1r1);
    T_Ez_Pz=KE*((-((y-dy).^2+(x-dx).^2))./r1./r1.*g1r1+g2r1);
    T_Hx_Pz=KH*(-(y-dy)./r1.*g3r1);
    T_Hy_Pz=KH*((x-dx)./r1.*g3r1); 
    T_Hz_Pz=zeros(pointsnum, 1);

    T_Ex_Mx=zeros(pointsnum, 1);
    T_Ey_Mx=KE*(-(z-dz)./r1.*g3r1).*k0;
    T_Ez_Mx=KE*((y-dy)./r1.*g3r1).*k0;
    T_Hx_Mx=KH*((-((y-dy).^2+(z-dz).^2))./r1./r1.*g1r1+g2r1).*k0;
    T_Hy_Mx=KH*((x-dx).*(y-dy)./r1./r1.*g1r1).*k0;
    T_Hz_Mx=KH*((z-dz).*(x-dx)./r1./r1.*g1r1).*k0;

    T_Ex_My=KE*((z-dz)./r1.*g3r1).*k0;  
    T_Ey_My=zeros(pointsnum, 1);
    T_Ez_My=KE*(-(x-dx)./r1.*g3r1).*k0;
    T_Hx_My=KH*((x-dx).*(y-dy)./r1./r1.*g1r1).*k0;
    T_Hy_My=KH*(-((x-dx).^2+(z-dz).^2)./r1./r1.*g1r1+g2r1).*k0;  
    T_Hz_My=KH*((y-dy).*(z-dz)./r1./r1.*g1r1).*k0;  

    T_Ex_Mz=KE*(-(y-dy)./r1.*g3r1).*k0;
    T_Ey_Mz=KE*((x-dx)./r1.*g3r1).*k0;
    T_Ez_Mz=zeros(pointsnum, 1);
    T_Hx_Mz=KH*((z-dz).*(x-dx)./r1./r1.*g1r1).*k0;
    T_Hy_Mz=KH*((y-dy).*(z-dz)./r1./r1.*g1r1).*k0;
    T_Hz_Mz=KH*(-((x-dx).^2+(y-dy).^2)./r1./r1.*g1r1+g2r1).*k0;    
    end
    T=[];        
        if m==1
            if fieldOpt(1)==1
                T=[T;T_Ex_Px];
            end
            if fieldOpt(2)==1
                T=[T;T_Ey_Px];
            end    
            if fieldOpt(3)==1
                T=[T;T_Ez_Px];
            end
            if fieldOpt(4)==1
                T=[T;T_Hx_Px];
            end 
            if fieldOpt(5)==1
                T=[T;T_Hy_Px];
            end
            if fieldOpt(6)==1
                T=[T;T_Hz_Px];
            end     
            Tmatrix(:,i)=T;
        elseif m==2
            if fieldOpt(1)==1
                T=[T;T_Ex_Py];
            end
            if fieldOpt(2)==1
                T=[T;T_Ey_Py];
            end    
            if fieldOpt(3)==1
                T=[T;T_Ez_Py];
            end
            if fieldOpt(4)==1
                T=[T;T_Hx_Py];
            end 
            if fieldOpt(5)==1
                T=[T;T_Hy_Py];
            end
            if fieldOpt(6)==1
                T=[T;T_Hz_Py];
            end     
            Tmatrix(:,i)=T;
        elseif m==3
            if fieldOpt(1)==1
                T=[T;T_Ex_Pz];
            end
            if fieldOpt(2)==1
                T=[T;T_Ey_Pz];
            end    
            if fieldOpt(3)==1
                T=[T;T_Ez_Pz];
            end
            if fieldOpt(4)==1
                T=[T;T_Hx_Pz];
            end 
            if fieldOpt(5)==1
                T=[T;T_Hy_Pz];
            end
            if fieldOpt(6)==1
                T=[T;T_Hz_Pz];
            end     
            Tmatrix(:,i)=T;
        elseif m==4
            if fieldOpt(1)==1
                T=[T;T_Ex_Mx];
            end
            if fieldOpt(2)==1
                T=[T;T_Ey_Mx];
            end    
            if fieldOpt(3)==1
                T=[T;T_Ez_Mx];
            end
            if fieldOpt(4)==1
                T=[T;T_Hx_Mx];
            end 
            if fieldOpt(5)==1
                T=[T;T_Hy_Mx];
            end
            if fieldOpt(6)==1
                T=[T;T_Hz_Mx];
            end     
            Tmatrix(:,i)=T;            
            Constant(i,1)=c_m_dipole;   
        elseif m==5
            if fieldOpt(1)==1
                T=[T;T_Ex_My];
            end
            if fieldOpt(2)==1
                T=[T;T_Ey_My];
            end    
            if fieldOpt(3)==1
                T=[T;T_Ez_My];
            end
            if fieldOpt(4)==1
                T=[T;T_Hx_My];
            end 
            if fieldOpt(5)==1
                T=[T;T_Hy_My];
            end
            if fieldOpt(6)==1
                T=[T;T_Hz_My];
            end     
            Tmatrix(:,i)=T;            
            Constant(i,1)=c_m_dipole;
        elseif m==6 
            if fieldOpt(1)==1
                T=[T;T_Ex_Mz];
            end
            if fieldOpt(2)==1
                T=[T;T_Ey_Mz];
            end    
            if fieldOpt(3)==1
                T=[T;T_Ez_Mz];
            end
            if fieldOpt(4)==1
                T=[T;T_Hx_Mz];
            end 
            if fieldOpt(5)==1
                T=[T;T_Hy_Mz];
            end
            if fieldOpt(6)==1
                T=[T;T_Hz_Mz];
            end     
            Tmatrix(:,i)=T;  
            Constant(i,1)=c_m_dipole; 
        end
end

DipoleQuantity0 = ((Tmatrix')*Tmatrix)\(Tmatrix')*ScanField; %M dipole unit: A*m*m
Dipole_Quantity = DipoleQuantity0.*Constant; %M dipole unit: V*m;  Constant=1j*n0*2*pi*freq/c;
Dipole_Magnitude_HFSS=abs(Dipole_Quantity);
Dipole_Phase_HFSS=180*angle(Dipole_Quantity)/pi;% Should be 'Dipole_Phase=180*angle(Dipole_Quantity0)/pi'
Dipole_Magnitude0=abs(DipoleQuantity0);
Dipole_Phase0=180*angle(DipoleQuantity0)/pi;
if Unit==1
Dipole_Magnitude=abs(Dipole_Quantity);
Dipole_Phase=Dipole_Phase_HFSS;
else
Dipole_Magnitude=abs(DipoleQuantity0);
Dipole_Phase=Dipole_Phase0;
end
set(handles.editOutput2,'String',num2str(Dipole_Magnitude'));
set(handles.edit11,'String',num2str(Dipole_Phase'));




%Write CSV File
Type=[];
for i=1:DipoleNum
    if DipoleType(i,1)==1
        Type=[Type;'Px'];
    elseif DipoleType(i,2)==1
        Type=[Type;'Py'];
    elseif DipoleType(i,3)==1
        Type=[Type;'Pz'];        
    elseif DipoleType(i,4)==1
        Type=[Type;'Mx'];
    elseif DipoleType(i,5)==1
        Type=[Type;'My']; 
    elseif DipoleType(i,6)==1
        Type=[Type;'Mz'];  
    end
end
Magnitude=Dipole_Magnitude;
Phase_unit_Degree=Dipole_Phase;
Location_X_mm=1000*DipoleLocation_0(:,1);
Location_Y_mm=1000*DipoleLocation_0(:,2);
Location_Z_mm=1000*DipoleLocation_0(:,3);
%C={Type,Location_X_mm,Location_Y_mm,Location_Z_mm,Magnitude_Punit_Am_Munit_Vm,Phase_unit_Degree};
if DipoleNum==1
    Type={Type};
    T=table(Type,Location_X_mm,Location_Y_mm,Location_Z_mm,Magnitude,Phase_unit_Degree)
else
    T=table(Type,Location_X_mm,Location_Y_mm,Location_Z_mm,Magnitude,Phase_unit_Degree)
end
writetable(T,'ReconstructedDipole.csv');


%Write HFSS Script
fp = fopen('HFSS_Script.txt','w');
index_dipole = 1;
HistoryID = 20;
for i = 1:DipoleNum
            fprintf(fp,['$begin' char(39) 'IncHDWave' num2str(index_dipole) char(39),'\n']);
            fprintf(fp,['ID=',num2str(index_dipole+HistoryID),'\n']);            
            fprintf(fp,['BoundType=' char(39) 'Hertzian-Dipole Incident Wave' char(39) '\n']);
            fprintf(fp,'ParentBndID=-1\n');
            fprintf(fp,'IsCartesian=true\n'); 
            fprintf(fp,['kX='  char(39) '0' char(39) '\n']);
            fprintf(fp,['kY='  char(39) '0' char(39) '\n']);
            fprintf(fp,['kZ='  char(39) '1' char(39) '\n']);
            fprintf(fp,['OriginX='  char(39) num2str(DipoleLocation_0(i,1)) 'meter' char(39) '\n']);
            fprintf(fp,['OriginY='  char(39) num2str(DipoleLocation_0(i,2)) 'meter' char(39) '\n']);
            fprintf(fp,['OriginZ='  char(39) num2str(DipoleLocation_0(i,3)) 'meter' char(39) '\n']);
            fprintf(fp,['SphereRadius='  char(39) '1e-4meter' char(39) '\n']);
            if DipoleType(i,1)==1   
                fprintf(fp,['EoX=' char(39) '1' char(39) '\n']);
                fprintf(fp,['EoY=' char(39) '0' char(39) '\n']);
                fprintf(fp,['EoZ=' char(39) '0' char(39) '\n']);
                fprintf(fp,'IsElectricDipole=true\n');%electric dipole
            elseif DipoleType(i,2)==1   
                fprintf(fp,['EoX=' char(39) '0' char(39) '\n']);
                fprintf(fp,['EoY=' char(39) '1' char(39) '\n']);
                fprintf(fp,['EoZ=' char(39) '0' char(39) '\n']);
                fprintf(fp,'IsElectricDipole=true\n');%electric dipole
            elseif DipoleType(i,3)==1   
                fprintf(fp,['EoX=' char(39) '0' char(39) '\n']);
                fprintf(fp,['EoY=' char(39) '0' char(39) '\n']);
                fprintf(fp,['EoZ=' char(39) '1' char(39) '\n']);
                fprintf(fp,'IsElectricDipole=true\n');%electric dipole    
            elseif DipoleType(i,4)==1  
                fprintf(fp,['EoX=' char(39) '1' char(39) '\n']);
                fprintf(fp,['EoY=' char(39) '0' char(39) '\n']);
                fprintf(fp,['EoZ=' char(39) '0' char(39) '\n']);
                fprintf(fp,'IsElectricDipole=false\n');%magnetic dipole
            elseif DipoleType(i,5)==1  
                fprintf(fp,['EoX=' char(39) '0' char(39) '\n']);
                fprintf(fp,['EoY=' char(39) '1' char(39) '\n']);
                fprintf(fp,['EoZ=' char(39) '0' char(39) '\n']);
                fprintf(fp,'IsElectricDipole=false\n');%magnetic dipole
            elseif DipoleType(i,6)==1  
                fprintf(fp,['EoX=' char(39) '0' char(39) '\n']);
                fprintf(fp,['EoY=' char(39) '0' char(39) '\n']);
                fprintf(fp,['EoZ=' char(39) '1' char(39) '\n']);
                fprintf(fp,'IsElectricDipole=false\n');%magnetic dipole
            end  
            fprintf(fp,['$end' char(39) 'IncHDWave' num2str(index_dipole) char(39),'\n']);
            index_dipole=index_dipole+1;
end
        
index_dipole=1;
for i = 1:DipoleNum   
            fprintf(fp,['SourceEntry(ID=',num2str(index_dipole+HistoryID),', Index=0, Terminal=false, Terminated=false, Magnitude='...
                        char(39) num2str(Dipole_Magnitude_HFSS(i)) char(39),', phase='...
                        char(39) num2str(Dipole_Phase_HFSS(i)) 'deg' char(39) ')\n']);% e dipole 
            index_dipole=index_dipole+1;            

end
fclose(fp);

if SolutionModel==1
   FieldCalculator(freq,x_start,y_start,scanHeight,x_end-x_start,y_end-y_start,xSpcScan,ySpcScan,Dipole_Magnitude0,Dipole_Phase0); 
end

if SolutionModel==2
    FieldCalculator_GND(freq,x_start,y_start,scanHeight,x_end-x_start,y_end-y_start,xSpcScan,ySpcScan,Dipole_Magnitude0,Dipole_Phase0);
end

%Write CST Script
DipoleLocation=DipoleLocation_0;
P_len = 0.2;  % mm is unit, before it's 0.05, now 0.2 make it big, 0.5 is good
offset=0.2;
M_len=0.2; % mm is unit, before it's 0.05, now 0.2 make it big,0.5 is good
P_I = 1e3/P_len;
M_thick=0.005;
edge=0.002;
M_I=1e6/(M_len*M_len); %  for M dipole, I*dS=1

fp = fopen('CST_Geometry_Script.txt','w');

fprintf(fp,'Option Explicit\n');
fprintf(fp,'Sub Main\n');

if SolutionModel==2
fprintf(fp,'With Brick\n');
fprintf(fp,'.Reset\n');
fprintf(fp,'.Name ("GND")\n');
fprintf(fp,'.Component ("Ground")\n');
fprintf(fp,'.Material ("PEC")\n');
fprintf(fp,'.Xrange (-20, 20)\n');
fprintf(fp,'.Yrange (-20, 20)\n');
fprintf(fp,'.Zrange (-0.035, 0)\n');
fprintf(fp,'.Create\n');
fprintf(fp,'End With\n');
end

fprintf(fp,'With Units\n');
fprintf(fp,'.Geometry ("mm")\n');
fprintf(fp,'.Frequency ("ghz")\n');
fprintf(fp,'.Time ("s")\n');
fprintf(fp,'.TemperatureUnit ("kelvin")\n');
fprintf(fp,'End With\n');
fprintf(fp,'With Background\n');
fprintf(fp,'.ResetBackground\n');
fprintf(fp,'.XminSpace "0.0"\n');
fprintf(fp,'.XmaxSpace "0.0"\n');
fprintf(fp,'.YminSpace "0.0"\n');
fprintf(fp,'.YmaxSpace "0.0"\n');
fprintf(fp,'.ZminSpace "0.0"\n');
fprintf(fp,'.ZmaxSpace "0.0"\n');
fprintf(fp,'.ApplyInAllDirections "False"\n');
fprintf(fp,'End With\n');
fprintf(fp,'With Material\n');
fprintf(fp,'.Reset\n');
fprintf(fp,'.FrqType "all"\n');
fprintf(fp,'.ChangeBackgroundMaterial\n');
fprintf(fp,'End With\n');
fprintf(fp,'Dim temp As String\n');

for i = 1:DipoleNum
    if DipoleType(i,1)==1
    fprintf(fp,'temp ="With DiscretePort" + vbCrLf + _\n');
	fprintf(fp,'" .Reset" + vbCrLf + _\n');	
	fprintf(fp,['" .PortNumber "+CStr(',num2str(i),')+vbCrLf + _\n']); 
	fprintf(fp,'" .Type ""Current"" "+ vbCrLf + _\n');
	fprintf(fp,['" .Current "+CStr(',num2str(P_I),')+ vbCrLf + _\n']);
    fprintf(fp,['" .SetP1 ""False"", " + CStr(',num2str(-P_len/2+DipoleLocation(i,1)*1000+offset),')+", "+CStr(',num2str(DipoleLocation(i,2)*1000),')+", "+CStr(',num2str(DipoleLocation(i,3)*1000),') + vbCrLf + _\n']);
	fprintf(fp,['" .SetP2 ""False"", " + CStr(',num2str(P_len/2+DipoleLocation(i,1)*1000+offset),')+", "+CStr(',num2str(DipoleLocation(i,2)*1000),')+", "+CStr(',num2str(DipoleLocation(i,3)*1000),') + vbCrLf + _\n']);
	fprintf(fp,'" .Monitor ""True"" "+ vbCrLf + _\n');
	fprintf(fp,'" .Create" + vbCrLf + _\n');
	fprintf(fp,'"End With"\n');
	fprintf(fp,['AddtoHistory("define port:"+CStr(',num2str(i),'), temp)\n']);
    
     fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Px',num2str(i),'_','1")\n']);
         fprintf(fp,'.Component ("Px")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(-P_len/2+DipoleLocation(i,1)*1000+offset-edge),',', num2str(-P_len/2+DipoleLocation(i,1)*1000+offset),')\n']);
         fprintf(fp,['.Yrange (',num2str(-(edge/2)+DipoleLocation(i,2)*1000),',', num2str((edge/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(DipoleLocation(i,3)*1000),',', num2str(DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
         
         
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Px',num2str(i),'_','2")\n']);
         fprintf(fp,'.Component ("Px")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(P_len/2+DipoleLocation(i,1)*1000+offset),',', num2str(P_len/2+DipoleLocation(i,1)*1000+offset+edge),')\n']);
         fprintf(fp,['.Yrange (',num2str(-(edge/2)+DipoleLocation(i,2)*1000),',', num2str((edge/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(DipoleLocation(i,3)*1000),',', num2str(DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
         
         
    elseif DipoleType(i,2)==1
    fprintf(fp,'temp ="With DiscretePort" + vbCrLf + _\n');
	fprintf(fp,'" .Reset" + vbCrLf + _\n');	
	fprintf(fp,['" .PortNumber "+CStr(',num2str(i),')+vbCrLf + _\n']); 
	fprintf(fp,'" .Type ""Current"" "+ vbCrLf + _\n');
	fprintf(fp,['" .Current "+CStr(',num2str(P_I),')+ vbCrLf + _\n']);
    fprintf(fp,['" .SetP1 ""False"", " + CStr(',num2str(DipoleLocation(i,1)*1000),')+", "+CStr(',num2str(-P_len/2+DipoleLocation(i,2)*1000+offset),')+", "+CStr(',num2str(DipoleLocation(i,3)*1000),') + vbCrLf + _\n']);
	fprintf(fp,['" .SetP2 ""False"", " + CStr(',num2str(DipoleLocation(i,1)*1000),')+", "+CStr(',num2str(P_len/2+DipoleLocation(i,2)*1000+offset),')+", "+CStr(',num2str(DipoleLocation(i,3)*1000),') + vbCrLf + _\n']);
	fprintf(fp,'" .Monitor ""True"" "+ vbCrLf + _\n');
	fprintf(fp,'" .Create" + vbCrLf + _\n');
	fprintf(fp,'"End With"\n');
	fprintf(fp,['AddtoHistory("define port:"+CStr(',num2str(i),'), temp)\n']);  
    
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Py',num2str(i),'_','1")\n']);
         fprintf(fp,'.Component ("Py")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(-(edge/2)+DipoleLocation(i,1)*1000),',', num2str((edge/2)+DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(-P_len/2+DipoleLocation(i,2)*1000+offset-edge),',', num2str(-P_len/2+DipoleLocation(i,2)*1000+offset),')\n']);
         fprintf(fp,['.Zrange (',num2str(DipoleLocation(i,3)*1000),',', num2str(DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
         
         
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Py',num2str(i),'_','2")\n']);
         fprintf(fp,'.Component ("Py")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(-(edge/2)+DipoleLocation(i,1)*1000),',', num2str((edge/2)+DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(P_len/2+DipoleLocation(i,2)*1000+offset-edge),',', num2str(P_len/2+DipoleLocation(i,2)*1000+offset),')\n']);
         fprintf(fp,['.Zrange (',num2str(DipoleLocation(i,3)*1000),',', num2str(DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
    
    
    
    elseif DipoleType(i,3)==1
    fprintf(fp,'temp ="With DiscretePort" + vbCrLf + _\n');
	fprintf(fp,'" .Reset" + vbCrLf + _\n');	
	fprintf(fp,['" .PortNumber "+CStr(',num2str(i),')+vbCrLf + _\n']); 
	fprintf(fp,'" .Type ""Current"" "+ vbCrLf + _\n');
	fprintf(fp,['" .Current "+CStr(',num2str(P_I),')+ vbCrLf + _\n']);
    fprintf(fp,['" .SetP1 ""False"", " + CStr(',num2str(DipoleLocation(i,1)*1000),')+", "+CStr(',num2str(DipoleLocation(i,2)*1000),')+", "+CStr(',num2str(-P_len/2+DipoleLocation(i,3)*1000+offset),') + vbCrLf + _\n']);
	fprintf(fp,['" .SetP2 ""False"", " + CStr(',num2str(DipoleLocation(i,1)*1000),')+", "+CStr(',num2str(DipoleLocation(i,2)*1000),')+", "+CStr(',num2str(P_len/2+DipoleLocation(i,3)*1000+offset),') + vbCrLf + _\n']);
	fprintf(fp,'" .Monitor ""True"" "+ vbCrLf + _\n');
	fprintf(fp,'" .Create" + vbCrLf + _\n');
	fprintf(fp,'"End With"\n');
	fprintf(fp,['AddtoHistory("define port:"+CStr(',num2str(i),'), temp)\n']);    
    
             
    
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Pz',num2str(i),'_','1")\n']);
         fprintf(fp,'.Component ("Pz")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(DipoleLocation(i,1)*1000),',', num2str(DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(-(edge/2)+DipoleLocation(i,2)*1000),',', num2str((edge/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(-P_len/2+DipoleLocation(i,3)*1000+offset-edge),',', num2str(-P_len/2+DipoleLocation(i,3)*1000+offset),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
         
         
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Pz',num2str(i),'_','2")\n']);
         fprintf(fp,'.Component ("Pz")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(DipoleLocation(i,1)*1000),',', num2str(DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(-(edge/2)+DipoleLocation(i,2)*1000),',', num2str((edge/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(P_len/2+DipoleLocation(i,3)*1000+offset-edge),',', num2str(P_len/2+DipoleLocation(i,3)*1000+offset),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
    
    
    elseif DipoleType(i,4)==1
    fprintf(fp,'temp ="With DiscretePort" + vbCrLf + _\n');
	fprintf(fp,'" .Reset" + vbCrLf + _\n');	
	fprintf(fp,['" .PortNumber "+CStr(',num2str(i),')+vbCrLf + _\n']); 
	fprintf(fp,'" .Type ""Current"" "+ vbCrLf + _\n');
	fprintf(fp,['" .Current "+CStr(',num2str(M_I),')+ vbCrLf + _\n']);
	fprintf(fp,['" .SetP2 ""False"", " + CStr(',num2str(DipoleLocation(i,1)*1000-offset),')+", "+CStr(',num2str((M_len/2)+DipoleLocation(i,2)*1000),')+", "+CStr(',num2str((M_len/2)+DipoleLocation(i,3)*1000),') + vbCrLf + _\n']);
	fprintf(fp,['" .SetP1 ""False"", " + CStr(',num2str(DipoleLocation(i,1)*1000-offset),')+", "+CStr(',num2str((M_len/2)+DipoleLocation(i,2)*1000),')+", "+CStr(',num2str(-(M_len/2)+DipoleLocation(i,3)*1000),') + vbCrLf + _\n']);
	fprintf(fp,'" .Monitor ""True"" "+ vbCrLf + _\n');
	fprintf(fp,'" .Create" + vbCrLf + _\n');
	fprintf(fp,'"End With"\n');
	fprintf(fp,['AddtoHistory("define port:"+CStr(',num2str(i),'), temp)\n']);
         % y=0.1 to -0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Mx',num2str(i),'_','1")\n']);
         fprintf(fp,'.Component ("Mx")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(DipoleLocation(i,1)*1000-M_thick-offset),',', num2str(DipoleLocation(i,1)*1000+M_thick-offset),')\n']);
         fprintf(fp,['.Yrange (',num2str((M_len/2)+DipoleLocation(i,2)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(-(M_len/2)+DipoleLocation(i,3)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
         % z=-0.1 to 0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Mx',num2str(i),'_','2")\n']);
         fprintf(fp,'.Component ("Mx")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(DipoleLocation(i,1)*1000-M_thick-offset),',', num2str(DipoleLocation(i,1)*1000+M_thick-offset),')\n']);
         fprintf(fp,['.Yrange (',num2str(-(M_len/2)+DipoleLocation(i,2)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(-(M_len/2)+DipoleLocation(i,3)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
      % y=-0.1 to 0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Mx',num2str(i),'_','3")\n']);
         fprintf(fp,'.Component ("Mx")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(DipoleLocation(i,1)*1000-M_thick-offset),',', num2str(DipoleLocation(i,1)*1000+M_thick-offset),')\n']);
         fprintf(fp,['.Yrange (',num2str(-(M_len/2)+DipoleLocation(i,2)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(+(M_len/2)+DipoleLocation(i,3)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');        
    elseif DipoleType(i,5)==1
    fprintf(fp,'temp ="With DiscretePort" + vbCrLf + _\n');
	fprintf(fp,'" .Reset" + vbCrLf + _\n');	
	fprintf(fp,['" .PortNumber "+CStr(',num2str(i),')+vbCrLf + _\n']); 
	fprintf(fp,'" .Type ""Current"" "+ vbCrLf + _\n');
	fprintf(fp,['" .Current "+CStr(',num2str(M_I),')+ vbCrLf + _\n']);
	fprintf(fp,['" .SetP2 ""False"", " + CStr(',num2str((M_len/2)+DipoleLocation(i,1)*1000),')+", "+CStr(',num2str(DipoleLocation(i,2)*1000-offset),')+", "+CStr(',num2str((M_len/2)+DipoleLocation(i,3)*1000),') + vbCrLf + _\n']);
	fprintf(fp,['" .SetP1 ""False"", " + CStr(',num2str(-(M_len/2)+DipoleLocation(i,1)*1000),')+", "+CStr(',num2str(DipoleLocation(i,2)*1000-offset),')+", "+CStr(',num2str((M_len/2)+DipoleLocation(i,3)*1000),') + vbCrLf + _\n']);
	fprintf(fp,'" .Monitor ""True"" "+ vbCrLf + _\n');
	fprintf(fp,'" .Create" + vbCrLf + _\n');
	fprintf(fp,'"End With"\n');
	fprintf(fp,['AddtoHistory("define port:"+CStr(',num2str(i),'), temp)\n']);
         % z=0.1 to -0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("My',num2str(i),'_','1")\n']);
         fprintf(fp,'.Component ("My")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(-(M_len/2)+DipoleLocation(i,1)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(DipoleLocation(i,2)*1000-M_thick-offset),',', num2str(DipoleLocation(i,2)*1000+M_thick-offset),')\n']);
         fprintf(fp,['.Zrange (',num2str((M_len/2)+DipoleLocation(i,3)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
         % x=-0.1 to 0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("My',num2str(i),'_','2")\n']);
         fprintf(fp,'.Component ("My")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(-(M_len/2)+DipoleLocation(i,1)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(DipoleLocation(i,2)*1000-M_thick-offset),',', num2str(DipoleLocation(i,2)*1000+M_thick-offset),')\n']);
         fprintf(fp,['.Zrange (',num2str(-(M_len/2)+DipoleLocation(i,3)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
      % z=-0.1 to 0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("My',num2str(i),'_','3")\n']);
         fprintf(fp,'.Component ("My")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(+(M_len/2)+DipoleLocation(i,1)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(DipoleLocation(i,2)*1000-M_thick-offset),',', num2str(DipoleLocation(i,2)*1000+M_thick-offset),')\n']);
         fprintf(fp,['.Zrange (',num2str(-(M_len/2)+DipoleLocation(i,3)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,3)*1000),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');         
    elseif DipoleType(i,6)==1
    fprintf(fp,'temp ="With DiscretePort" + vbCrLf + _\n');
	fprintf(fp,'" .Reset" + vbCrLf + _\n');	
	fprintf(fp,['" .PortNumber "+CStr(',num2str(i),')+vbCrLf + _\n']); 
	fprintf(fp,'" .Type ""Current"" "+ vbCrLf + _\n');
	fprintf(fp,['" .Current "+CStr(',num2str(M_I),')+ vbCrLf + _\n']);
	fprintf(fp,['" .SetP1 ""False"", " + CStr(',num2str((M_len/2)+DipoleLocation(i,1)*1000),')+", "+CStr(',num2str((M_len/2)+DipoleLocation(i,2)*1000),')+", "+CStr(',num2str(DipoleLocation(i,3)*1000-offset),') + vbCrLf + _\n']);
	fprintf(fp,['" .SetP2 ""False"", " + CStr(',num2str(-(M_len/2)+DipoleLocation(i,1)*1000),')+", "+CStr(',num2str((M_len/2)+DipoleLocation(i,2)*1000),')+", "+CStr(',num2str(DipoleLocation(i,3)*1000-offset),') + vbCrLf + _\n']);
	fprintf(fp,'" .Monitor ""True"" "+ vbCrLf + _\n');
	fprintf(fp,'" .Create" + vbCrLf + _\n');
	fprintf(fp,'"End With"\n');
	fprintf(fp,['AddtoHistory("define port:"+CStr(',num2str(i),'), temp)\n']);
         % y=0.1 to -0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Mz',num2str(i),'_','1")\n']);
         fprintf(fp,'.Component ("Mz")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(-(M_len/2)+DipoleLocation(i,1)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str((M_len/2)+DipoleLocation(i,2)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(DipoleLocation(i,3)*1000-M_thick-offset),',', num2str(DipoleLocation(i,3)*1000+M_thick-offset),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
         % x=-0.1 to 0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Mz',num2str(i),'_','2")\n']);
         fprintf(fp,'.Component ("Mz")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(-(M_len/2)+DipoleLocation(i,1)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(-(M_len/2)+DipoleLocation(i,2)*1000),',', num2str(-(M_len/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(DipoleLocation(i,3)*1000-M_thick-offset),',', num2str(DipoleLocation(i,3)*1000+M_thick-offset),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');
      % y=-0.1 to 0.1
         fprintf(fp,'With Brick\n');
         fprintf(fp,'.Reset\n');
         fprintf(fp,['.Name ("Mz',num2str(i),'_','3")\n']);
         fprintf(fp,'.Component ("Mz")\n');
         fprintf(fp,'.Material ("PEC")\n');
         fprintf(fp,['.Xrange (',num2str(+(M_len/2)+DipoleLocation(i,1)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,1)*1000),')\n']);
         fprintf(fp,['.Yrange (',num2str(-(M_len/2)+DipoleLocation(i,2)*1000),',', num2str(+(M_len/2)+DipoleLocation(i,2)*1000),')\n']);
         fprintf(fp,['.Zrange (',num2str(DipoleLocation(i,3)*1000-M_thick-offset),',', num2str(DipoleLocation(i,3)*1000+M_thick-offset),')\n']);
         fprintf(fp,'.Create\n');
         fprintf(fp,'End With\n');      
    end
end
fprintf(fp,'With Solver\n');
fprintf(fp,'.FrequencyRange ( 1.0, 2.0 )\n');
fprintf(fp,'.Method "Hexahedral"\n');
fprintf(fp,'.SteadyStateLimit "-30.0"\n');
fprintf(fp,'.NormingImpedance "50"\n');
fprintf(fp,'.SteadyStateDurationType "Number of pulses"\n');
fprintf(fp,'.NumberOfPulseWidths "300"\n');
fprintf(fp,'.SteadyStateDurationTime "1.06667e-008"\n');
fprintf(fp,'.SteadyStateDurationTimeAsDistance "3200"\n');
fprintf(fp,'.Start\n');
fprintf(fp,'End With\n');
fprintf(fp,'End Sub\n'); 
fclose(fp);

fp = fopen('CST_DipoleValue_Script.txt','w');

fprintf(fp,'Option Explicit\n');
fprintf(fp,'Sub Main\n');
fprintf(fp,'With CombineResults\n');
fprintf(fp,'.Reset\n');
fprintf(fp,'.SetMonitorType ("frequency")\n');
fprintf(fp,'.EnableAutomaticLabeling (False)\n');
fprintf(fp,'.SetLabel ("Code combination")\n');
fprintf(fp,'.SetNone\n');
fprintf(fp,'End With\n');
for i=1:DipoleNum
    fprintf(fp,['CombineResults.SetExcitationValues ("port",',num2str(i),',1,Cstr(',num2str(Dipole_Magnitude0(i)),'),',num2str(Dipole_Phase0(i)),')\n']);
end
fprintf(fp,'CombineResults.Run\n');
fprintf(fp,'End Sub\n');
fclose(fp);