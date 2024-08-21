%% Simulasi 3D untuk CFL=0.5
clc;clear all

% Inisialisasi Daerah Asal
Lx = 10; % Ukuran domain spasial x
Ly = 10; % Ukuran domain spasial y
dx = 0.1; % Langkah diskritisasi pada sumbu x
dy = dx; % Langkah diskritisasi pada sumbu y
nx = fix(Lx/dx);
ny = fix(Ly/dy);
x = linspace(0,Lx,nx);
y = linspace(0, Ly, ny);
T = 10; % Total waktu simulasi

% Field Variable
u = zeros(nx,ny);
um1 = u; % w pada waktu n-1
up1 = u; % w pada waktu n+1

% Parameter Lanjutan
CFL = 0.5;  % dengan CFL = c.dt/dx
c = 1;
dt = CFL*dx/c;

% Time Stepping Loop
t = 0;
while t<T
    % Kondisi Batas Reflektif
    u(:,[1,end])=0;
    u([1 end],:)=0;

    % Solusi
    t = t+dt;
    um1 = u;
    u = up1; % Menyimpan array saat ini dan sebelumnya

    % source 
    u(50,50)=dt^2*20*sin(30*pi*t/20);
    for i=2:nx-1
        for j=2:ny-1
            up1(i,j)=2*u(i,j)-um1(i,j)+CFL^2*(u(i+1,j) + u(i,j+1) - 4*u(i,j) + u(i-1,j) + u(i,j-1));
        end
    end

    % Visualisasi pada langkah tertentu
    clf;
    subplot(2,1,1);
    imagesc(x,y, u');
    colorbar;
    caxis([-0.02 0.02]);
    title(sprintf('t=%.2f',t));
    subplot(2,1,2);
    mesh(x,y,u');
    colorbar;
    caxis([-0.02 0.02]);
    axis([0 Lx 0 Ly -0.05 0.05]);
    shg;
    pause(0.01);
end

% Simulasi 2D untuk CFL=0.5
% Inisialisasi Derah Asal
Lx = 10; % Ukuran domain spasial x
Ly = 10; % Ukuran domain spasial y
dx = 0.1; % Langkah diskritisasi pada sumbu x
dy = dx; % Langkah diskritisasi pada sumbu y
nx = fix(Lx/dx);
ny = fix(Ly/dy);
x = linspace(0,Lx,nx);
y = linspace(0, Ly, ny);
T = 10; % Total waktu simulasi

% Field Variable
u = zeros(nx,ny);
um1 = u; % w pada waktu n-1
up1 = u; % w pada waktu n+1

% Parameter Lanjutan
CFL = 0.5;  % dengan CFL = c.dt/dx
c = 1;
dt = CFL*dx/c;

%Time Stepping Loop
t = 0;
while t<T

    % Kondisi Batas Reflektif
    u(:,[1,end])=0;
    u([1 end],:)=0;

    % Solusi
    t = t+dt;
    um1 = u;
    u = up1; % Menyimpan array saat ini dan sebelumnya

    % source 
    u(50,50)=dt^2*20*sin(30*pi*t/20);
    for i=2:nx-1
        for j=2:ny-1
            up1(i,j)=2*u(i,j)-um1(i,j)+CFL^2*(u(i+1,j) + u(i,j+1) - 4*u(i,j) + u(i-1,j) + u(i,j-1));
        end
    end

    % Visualisasi pada langkah tertentu
    clf;
    subplot(2,1,1);
    plot(x,u');
    title(sprintf('t=%.2f',t));
    subplot(2,1,2);
    plot(y,u');
    title(sprintf('t=%.2f',t))
    shg;
    pause(0.01);
end

%% Simulasi 3D untuk CFL=0.8
clc;clear all

% Inisialisasi Daerah Asal
Lx = 10; % Ukuran domain spasial x
Ly = 10; % Ukuran domain spasial y
dx = 0.1; % Langkah diskritisasi pada sumbu x
dy = dx; % Langkah diskritisasi pada sumbu y
nx = fix(Lx/dx);
ny = fix(Ly/dy);
x = linspace(0,Lx,nx);
y = linspace(0, Ly, ny);
T = 10; % Total waktu simulasi

% Field Variable
u = zeros(nx,ny);
um1 = u; % w pada waktu n-1
up1 = u; % w pada waktu n+1

% Parameter Lanjutan
CFL = 0.8;  % dengan CFL = c.dt/dx
c = 1;
dt = CFL*dx/c;

% Time stepping loop
t = 0;
while t<T

    % Kondisi Batas Reflektif
    u(:,[1,end])=0;
    u([1 end],:)=0;

    % Solusi
    t = t+dt;
    um1 = u;
    u = up1; % Menyimpan array saat ini dan sebelumnya

    % source 
    u(50,50)=dt^2*20*sin(30*pi*t/20);
    for i=2:nx-1
        for j=2:ny-1
            up1(i,j)=2*u(i,j)-um1(i,j)+CFL^2*(u(i+1,j) + u(i,j+1) - 4*u(i,j) + u(i-1,j) + u(i,j-1));
        end
    end

    % Visualisasi pada langkah tertentu
    clf;
    subplot(2,1,1);
    imagesc(x,y, u');
    colorbar;
    caxis([-0.02 0.02]);
    title(sprintf('t=%.2f',t));
    subplot(2,1,2);
    mesh(x,y,u');
    colorbar;
    caxis([-0.02 0.02]);
    axis([0 Lx 0 Ly -0.05 0.05]);
    shg;
    pause(0.01);
end

% Simulasi untuk CFL=0.8
% Inisialisasi Daerah Asal
Lx = 10; % Ukuran domain spasial x
Ly = 10; % Ukuran domain spasial y
dx = 0.1; % Langkah diskritisasi pada sumbu x
dy = dx; % Langkah diskritisasi pada sumbu y
nx = fix(Lx/dx);
ny = fix(Ly/dy);
x = linspace(0,Lx,nx);
y = linspace(0, Ly, ny);
T = 10; % Total waktu simulasi

% Field Variable
u = zeros(nx,ny);
um1 = u; % w pada waktu n-1
up1 = u; % w pada waktu n+1

% Parameter Lanjutan
CFL = 0.8;  % dengan CFL = c.dt/dx
c = 1;
dt = CFL*dx/c;

% Time stepping loop
t = 0;
while t<T

    % Kondisi Batas Reflektif
    u(:,[1,end])=0;
    u([1 end],:)=0;

    % Solusi
    t = t+dt;
    um1 = u;
    u = up1; % Menyimpan array saat ini dan sebelumnya

    % source 
    u(50,50)=dt^2*20*sin(30*pi*t/20);
    for i=2:nx-1
        for j=2:ny-1
            up1(i,j)=2*u(i,j)-um1(i,j)+CFL^2*(u(i+1,j) + u(i,j+1) - 4*u(i,j) + u(i-1,j) + u(i,j-1));
        end
    end

    % Visualisasi pada langkah tertentu
    clf;
    subplot(2,1,1);
    plot(x,u');
    title(sprintf('t=%.2f',t));
    subplot(2,1,2);
    plot(y,u');
    title(sprintf('t=%.2f',t))
    shg;
    pause(0.01);
end