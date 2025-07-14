clear; clc;

% --- 1. Grid setup ---
Nx = 96; Ny = 96; Nz = 96; dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% --- 2. Medium properties ---
medium.sound_speed = 1500;
medium.density = 1000;

% --- 3. Time array (80 µs) ---
kgrid.makeTime(medium.sound_speed, [], 80e-6);
t = kgrid.t_array; Nt = length(t);

% --- 4. Define focal point ---
focus = [Nx/2, Ny/2, 40];  % [x, y, z]

% --- 5. Create sensor mask at source plane (z = 1) ---
% sensor.mask = zeros(Nx, Ny, Nz);
% sensor.mask(:, :, 1) = 1;  % source plane

% --- Circular transducer mask: 44 mm diameter at z = 0 ---
radius_mm = 22;
[X, Y] = meshgrid(1:Nx, 1:Ny);
cx = Nx / 2; cy = Ny / 2;

mask2D = sqrt((X - cx).^2 + (Y - cy).^2) <= radius_mm;
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(:, :, 1) = mask2D';


sensor.record = {'p'};

% --- 6. Point source at focus ---
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(focus(1), focus(2), focus(3)) = 1;

f0 = 500e3;
num_cycles = 4;
source_freq = toneBurst(1/kgrid.dt, f0, num_cycles);

source.p = source_freq;

% --- 7. Forward simulation: record at source plane ---
sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
    'DataCast', 'single', ...
    'PMLInside', false, ...
    'PMLSize', 10, ...
    'PlotSim', false, ...
    'BinaryName', 'kspaceFirstOrder-CUDA.exe');

% --- 8. Time reversal: use recorded signals as new source ---
source_tr.p_mask = sensor.mask;
source_tr.p = fliplr(sensor_data.p);  % reverse in time
source_tr.p_mode = 'additive';

sensor_tr.mask = ones(Nx, Ny, Nz);  % observe full field

sensor_data_tr = kspaceFirstOrder3DG(kgrid, medium, source_tr, sensor_tr, ...
    'DataCast', 'single', ...
    'PMLInside', false, ...
    'PMLSize', 10, ...
    'PlotSim', false, ...
    'BinaryName', 'kspaceFirstOrder-CUDA.exe');

% --- 9. Visualisation ---
p_max = max(sensor_data_tr, [], 2);
p_map = reshape(p_max, Nx, Ny, Nz);

% Axial
figure;
imagesc(squeeze(p_map(:, :, focus(3))));
axis image; colormap hot; colorbar;
title('Axial Slice at Focal Depth (Z)');

% 3D isosurface
thresh = 0.8 * prctile(p_map(:), 99);
fv = isosurface(p_map, thresh);
figure;
p = patch(fv);
isonormals(p_map, p);
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1]); view(3); camlight; lighting gouraud;
title('3D Isosurface of Time-Reversed Focus');

% --- 10. (Optional) Line plot along Z axis ---
line_profile = squeeze(p_map(focus(1), focus(2), :));
figure;
plot((1:Nz)*dx*1e3, line_profile);
xlabel('Depth (mm)'); ylabel('Max Pressure');
title('Axial Beam Profile from Time-Reversed Signal');


% Sagittal (Y-Z) slice through centre
figure;
imagesc(squeeze(p_map(focus(1), :, :))');  % transpose for correct orientation
axis image; colormap hot; colorbar;
xlabel('Y'); ylabel('Z');
title('Sagittal Slice (Y-Z view)');



% --- 11. Extract arrival times from transducer plane (z = 0) ---
sensor_mask = sensor.mask(:, :, 1);
[ix, iy] = find(sensor_mask);
num_pts = length(ix);

arrival_time = zeros(Nx, Ny);  % initialise
for n = 1:num_pts
    p_trace = sensor_data.p(n, :);  % nth recorded signal
    [~, max_idx] = max(p_trace);    % index of peak
    arrival_time(ix(n), iy(n)) = kgrid.t_array(max_idx);
end

% --- 12. Convert arrival times to lens thickness map ---
delay = arrival_time - min(arrival_time(:));  % in seconds
c0 = medium.sound_speed;
c_lens = 2500;  % speed of sound in lens material (e.g. plastic)

% Lens thickness equation from group delay
thickness = delay .* (c0 * c_lens) ./ (c_lens - c0);  % in metres

% --- 13. Visualise and scale thickness map ---
thickness_mm = thickness * 1e3;
thickness_mm = thickness_mm - min(thickness_mm(:));  % zero baseline

% Scale to max height (e.g. 5 mm max)
thickness_mm = thickness_mm / max(thickness_mm(:)) * 20;

figure;
imagesc(thickness_mm); axis image; colorbar;
title('Lens Thickness Map (mm)');


% CLEAN SCRIPT: Crop lens to 44 mm circular region and export as PNG + 3D surface

% --- CONFIGURATION ---
dx = 1e-3;              % Grid spacing [m]
crop_mm = 48;           % Desired lens diameter [mm]
crop_px = crop_mm / (dx * 1e3);  % Convert mm to pixels (assuming dx = 1 mm)
assert(mod(crop_px,1)==0, 'Crop size must be an integer number of pixels');
crop_px = round(crop_px);

% --- Load or define thickness_mm (replace with your own if needed) ---
% Example placeholder:
% thickness_mm = rand(96, 96) * 5;  % Random dummy data for testing
% Replace this with your actual lens thickness map

% 1. Center crop the square around the lens
[Nx, Ny] = size(thickness_mm);
cx = round(Nx / 2); cy = round(Ny / 2);
half = crop_px / 2;
x_idx = round(cx - half + 1 : cx + half);
y_idx = round(cy - half + 1 : cy + half);

lens_crop = thickness_mm(x_idx, y_idx);

% 2. Apply circular mask
[Xc, Yc] = meshgrid(1:crop_px, 1:crop_px);
rc = crop_px / 2;
mask_circle = sqrt((Xc - rc).^2 + (Yc - rc).^2) <= rc;

lens_circular = lens_crop;
lens_circular(~mask_circle) = 0;

% 3. Export grayscale PNG height map
norm_lens = (lens_circular - min(lens_circular(:))) ...
           / (max(lens_circular(:)) - min(lens_circular(:)));
img_out = uint8(255 * norm_lens);
imwrite(img_out, 'lens_heightmap_44mm_circle.png');

figure;
imshow(img_out);
title('Circular Lens Height Map (PNG, 44 mm)');

% 4. 3D surface plot
[x_mm, y_mm] = meshgrid((1:crop_px) * dx * 1e3);

figure;
surf(x_mm, y_mm, lens_circular');
shading interp;
colormap turbo;
colorbar;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Lens Thickness (mm)');
title('3D Surface of 44 mm Circular Lens');
view(45, 30); axis equal;



% --- Save lens height map as grayscale PNG ---
% height_img = uint8(255 * mat2gray(thickness_mm));  % scale to [0, 255]

% --- Manual normalisation to [0, 255] without Image Toolbox ---
min_val = min(thickness_mm(:));
max_val = max(thickness_mm(:));
norm_img = (thickness_mm - min_val) / (max_val - min_val);  % scale to 0–1
height_img = uint8(255 * norm_img);  % scale to 0–255

imwrite(height_img, 'lens_heightmap.png');

figure;
imshow(height_img);
title('2D Height Map of Lens (PNG)');


;


% --- 3D surface plot of lens ---
[xs, ys] = meshgrid((1:Nx) * dx * 1e3, (1:Ny) * dx * 1e3);  % mm scale

figure;
surf(xs, ys, thickness_mm');
shading interp;
colormap turbo;
colorbar;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Lens Thickness (mm)');
title('3D Surface Plot of Acoustic Lens');
view(45, 30); axis tight;

