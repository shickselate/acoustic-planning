clear; clc;

% --- 1. Grid setup ---
Nx = 96; Ny = 96; Nz = 96; dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% --- 2. Medium ---
medium.sound_speed = 1500;
medium.density = 1000;

% --- 3. Time array (60 µs) ---
kgrid.makeTime(medium.sound_speed, [], 60e-6);
t = kgrid.t_array; Nt = length(t);

% --- 4. Source plane ---
source.p_mask = zeros(Nx, Ny, Nz);
% source.p_mask(:, :, 1) = 1;

% Circular aperture
[X, Y] = meshgrid(1:Nx, 1:Ny);
cx = Nx/2; cy = Ny/2;
mask2D = sqrt((X - cx).^2 + (Y - cy).^2) <= Nx * 0.45;
source.p_mask(:, :, 1) = mask2D';

% --- 5. Focus point ---
focus = [Nx/2, Ny/2, 40];  % (x, y, z)
f0 = 500e3;  % 500 kHz
c = medium.sound_speed;

% --- 6. Window function (Tukey-like) ---
win_x = myTukeywin(Nx, 0.6);
win_y = myTukeywin(Ny, 0.6);
apod = win_x * win_y';  % 2D apodisation window

% --- 7. Create time-delayed signals ---
active = find(source.p_mask);
source.p = zeros(length(active), Nt);

for n = 1:length(active)
    [ix, iy, iz] = ind2sub(size(source.p_mask), active(n));
    dist = sqrt((ix - focus(1))^2 + (iy - focus(2))^2 + (iz - focus(3))^2);
    delay = (dist * dx) / c;
    amp = apod(ix, iy) * 1e5;  % apodised + scaled
    source.p(n, :) = amp * sin(2 * pi * f0 * (t - delay));
end

% --- 8. Sensor across full grid ---
sensor.mask = ones(Nx, Ny, Nz);

% --- 9. Run GPU simulation ---
sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
    'DataCast', 'single', ...
    'PMLInside', false, ...
    'PMLSize', 10, ...
    'PlotSim', false, ...
    'BinaryName', 'kspaceFirstOrder-CUDA.exe');

% --- 10. Visualisation ---
p_max = max(sensor_data, [], 2);
p_map = reshape(p_max, Nx, Ny, Nz);

% Axial
figure;
imagesc(squeeze(p_map(:, :, focus(3))));
axis image; colormap hot; colorbar;
title('Axial Slice at Focal Depth (Z)');

% Sagittal
figure;
imagesc(squeeze(p_map(focus(1), :, :))');
axis image; colormap hot; colorbar;
title('Sagittal Slice (Y-Z view)');

% Coronal
figure;
imagesc(squeeze(p_map(:, focus(2), :))');
axis image; colormap hot; colorbar;
title('Coronal Slice (X-Z view)');

% 3D isosurface
thresh = 0.8 * prctile(p_map(:), 99);  % adaptive threshold
fv = isosurface(p_map, thresh);
figure;
p = patch(fv);
isonormals(p_map, p);
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1]); view(3); camlight; lighting gouraud;
title('3D Isosurface of Focused Pressure');

% --- 11. Custom Tukey window function ---
function w = myTukeywin(N, r)
    if r <= 0
        w = ones(N,1);
    elseif r >= 1
        w = hann(N);
    else
        w = zeros(N,1);
        L = floor(r*(N-1)/2);
        for n = 1:N
            if n <= L
                w(n) = 0.5 * (1 + cos(pi * ((2*n)/(r*(N-1)) - 1)));
            elseif n > L && n <= N - L
                w(n) = 1;
            else
                w(n) = 0.5 * (1 + cos(pi * ((2*n)/(r*(N-1)) - (2/r) + 1)));
            end
        end
    end
end
