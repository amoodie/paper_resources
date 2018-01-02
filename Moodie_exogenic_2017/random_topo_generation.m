% Andrew J. Moodie
% construct synthetic topography to play with in GMT and show how filtering works.
clear all; 
close all;

n = 2000; % size of DEM
base = 20; % base elevation
t = ones(n,n) * base;

%% MAKE DEM
% add y ridge
riy = 1000;
offset = 150;
conpick = 900;
ss = t(:, conpick)' + (sin( 0.094 * (1:n) - 1.1) + 1) * 40 * ((n - (n + 250 - conpick) + 1))/2000;
for i = 1:n;
    tt = t(:,i)' + (sin( 0.094 * (1:n) - 1.1) + 1) * 40 * ((n - (n + 250 - i) + 1))/2000;
    if offset < n
        t(riy:riy+60, i) =  tt(600:660);
        t(riy:riy+60, i) =  ss(600:660); % comment this line for a solid bar ridge
        t(riy-90:riy-30, offset) =  tt(600:660);
        offset = offset + 1;
    end
    if offset == n
       offset = 1; 
    end
end
t(riy-90:riy-30, 1:450) = repmat(base, 61, 450);
k = t(riy-90:riy+60, :); % retain to replace after bulges
       
% apply a regional tilting across area of interest
for i = 1:n;
    t(:,i) = t(:,i)' + abs(sin( 0.0017 * (1:n) -0) + 0) * 100;
end

% add pseudo-random roughness elements
% add large elements
b.n = 200; % number to add
b.x = randi([100 1400], 1, b.n); % rand x coordinate (within range)
b.y = randi([100 1400],1,b.n); % rand y coordinate (within range)
b.d = randi([300 600],1,b.n); % rand diameter (within range)
b.e = randi([200 400],1,b.n); % elevation multiplier 
for i = 1:b.n;
    b.s = fspecial('gaussian', b.d(i), (b.d(i)/6.5)) * b.e(i)^2; % build element to randomly gen specs
    t(b.x(i):(b.x(i)+b.d(i)-1), b.y(i):(b.y(i)+b.d(i)-1)) = t(b.x(i):(b.x(i)+b.d(i)-1) , b.y(i):(b.y(i)+b.d(i)-1)) + (b.s); % replace the topography at x and y with the element
end

% add small elements (noise)
z = 1.8*randn(n); % random noise
h = fspecial('average',7); % smooth morphology
Z = conv2(z, h); % smooth it a bit 
t = t + Z(1:n,1:n); % add the noise to the DEM

% add medium elements
md.n = 10000; % number of mounds
md.dmi = 30; % min diam
md.dma = 250; % max diam
md.emi = 200; % min elevation multiplier
md.ema = 700; % max elevation multiplier
md.d = randi([md.dmi md.dma], 1, md.n); % diam of mounds
md.c = randi([(1+md.dma) (n-md.dma)], md.n, 2); % coordinates of mounds; X, Y in columns; represents the topleft corner
md.e = randi([md.emi md.ema], 1, md.n); % elevation multipliers
for mm = 1:md.n;
    md.s = fspecial('gaussian', md.d(mm), (md.d(mm)/5)) * md.e(mm); % mound structure
    t(md.c(mm,1):md.c(mm,1)+md.d(mm)-1, md.c(mm,2):md.c(mm,2)+md.d(mm)-1) = t(md.c(mm,1):md.c(mm,1)+md.d(mm)-1, md.c(mm,2):md.c(mm,2)+md.d(mm)-1) + md.s; % replace
end

t = rot90(t, -1); % rotate the entire DEM (just for visualizing)

%% VIEW DEM
% 2D map using view
figure('Visible','on'); % change to on for figure
surf(t,'EdgeColor','None');
view(2);

%% TO NETCDF
% https://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2
grdwrite2(1:n, 1:n, t, 'high_valley.nc')
