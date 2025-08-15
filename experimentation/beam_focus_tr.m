%% FOCUSED ULTRASOUND (k-Wave 1.4) — Time-Reversal + Clean Visuals + Lens Map
% Stephen Hicks · 500 kHz · 44 mm circular aperture · z_f = 70 mm · dx = 0.5 mm

clear; clc;
addpath(genpath('k-Wave'));

%% -------------------- CONFIG --------------------
c0      = 1500;          % m/s (water)
rho0    = 1000;          % kg/m^3
f0      = 500e3;         % Hz
z_f_mm  = 70;            % mm focal depth (set 40 if you wish)
D_mm    = 44;            % mm aperture diameter
dx      = 0.5e-3;        % m grid step  (0.5 mm → ~6 PPW at 500 kHz)
ppw     = (c0/f0) / dx;  % points per wavelength
assert(ppw >= 6, 'PPW too low (%.2f). Increase grid resolution (smaller dx).', ppw);

% Domain size (keep tight but with margin for PML)
x_span_mm = 100; y_span_mm = 100; z_span_mm = 120;     % physical span (mm)
Nx = ceil(x_span_mm/(dx*1e3)); Ny = ceil(y_span_mm/(dx*1e3)); Nz = ceil(z_span_mm/(dx*1e3));
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% Medium
medium.sound_speed = c0;
medium.density     = rho0;

% Time: long enough to reach focus and back a bit
t_end = 2 * (z_f_mm/1000)/c0;    % s
kgrid.makeTime(c0, [], t_end);
Nt = length(kgrid.t_array);

% Aperture mask (z = first plane)
D_px   = round(D_mm/(dx*1e3));
R_px   = D_px/2;
[Xp, Yp] = meshgrid(1:Nx, 1:Ny);
cx = Nx/2; cy = Ny/2;
circ = sqrt((Xp - cx).^2 + (Yp - cy).^2) <= R_px;

% Strong lateral apodisation (no toolboxes)
tuk = myTukeywin(Nx, 0.7) * myTukeywin(Ny, 0.7)';  % 2D window
apod2D = single(tuk .* circ);

% Focal voxel
fz = round(z_f_mm/(dx*1e3));
fx = round(Nx/2); fy = round(Ny/2);

%% -------------------- TIME-REVERSAL: forward step --------------------
% Point source at focus emits a short tone burst
num_cycles = 5;
tb = toneBurst(1/kgrid.dt, f0, num_cycles);                     % row vector

source = struct();
source.p_mask = false(Nx,Ny,Nz);
source.p_mask(fx,fy,fz) = true;
source.p = tb;                                                  % same for single-point mask

% Record pressure on the aperture plane
sensor = struct();
sensor.mask = false(Nx,Ny,Nz);
sensor.mask(:,:,1) = circ';                                     % circular aperture at z=1
sensor.record = {'p'};

sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor);  % CPU solver (safe)
% GPU alt (if installed): kspaceFirstOrder3DG(...)

%% -------------------- TIME-REVERSAL: playback step --------------------
% Time reverse and play back from the aperture (optionally apodised)
source_tr = struct();
source_tr.p_mask = sensor.mask;
p_rec = sensor_data.p;                                          % [numActive x Nt]
p_rev = fliplr(p_rec);                                          % time reverse

% Apply apodisation across the aperture elements
act = find(sensor.mask);
[ax, ay, ~] = ind2sub(size(sensor.mask), act);
w = apod2D(sub2ind([Nx Ny], ax, ay));
source_tr.p = p_rev .* single(w);                               % weight each channel
source_tr.p_mode = 'additive';

% Observe full field
sensor_tr = struct();
sensor_tr.mask = true(Nx,Ny,Nz);

sd_tr = kspaceFirstOrder3D(kgrid, medium, source_tr, sensor_tr);

%% -------------------- VISUALISATION (time-gated) --------------------
% Time gate around expected focus arrival
t_focus = (z_f_mm/1000)/c0;                    % s
gate = (kgrid.t_array >= (t_focus - 15e-6)) & (kgrid.t_array <= (t_focus + 15e-6));

sd_resh = reshape(sd_tr, [], Nt);
p_gate_max = max(sd_resh(:, gate), [], 2);     % max within gate only
p_map = reshape(p_gate_max, Nx, Ny, Nz);

% Axial: XZ through centre Y
figure; imagesc((1:Nz)*dx*1e3, (1:Nx)*dx*1e3, squeeze(p_map(:, round(Ny/2), :)));
axis image; colormap hot; colorbar;
xlabel('Z (mm)'); ylabel('X (mm)'); title('Axial Slice (X–Z), time-gated');

% Coronal: YZ through centre X
figure; imagesc((1:Nz)*dx*1e3, (1:Ny)*dx*1e3, squeeze(p_map(round(Nx/2), :, :))');
axis image; colormap hot; colorbar;
xlabel('Z (mm)'); ylabel('Y (mm)'); title('Coronal Slice (Y–Z), time-gated');

% Focal spot: XY at z = fz
figure; imagesc((1:Nx)*dx*1e3, (1:Ny)*dx*1e3, p_map(:,:,fz)');
axis image; colormap hot; colorbar;
xlabel('X (mm)'); ylabel('Y (mm)'); title(sprintf('Focal Spot (X–Y @ Z = %d mm)', z_f_mm));

% Axial line profile on axis
figure; plot((1:Nz)*dx*1e3, squeeze(p_map(fx,fy,:)),'LineWidth',1.5);
xlabel('Z (mm)'); ylabel('Time-gated max pressure'); grid on;
title('On-axis axial profile');

%% -------------------- LENS THICKNESS MAP (optional) --------------------
% Compute arrival-time map on aperture from the FORWARD step (sensor_data.p).
% Using time index of peak arrival per element.
arrival_time = nan(Nx,Ny);
act = find(sensor.mask);
[ax, ay, ~] = ind2sub(size(sensor.mask), act);
for k = 1:numel(act)
    ptrace = sensor_data.p(k,:);                % pressure vs time for element k
    [~, idx] = max(ptrace);
    arrival_time(ax(k), ay(k)) = kgrid.t_array(idx);
end

% Clean: set outside aperture to NaN
arrival_time(~circ) = NaN;

% Convert delay to thickness for a lens material (phase plate)
c_lens = 2500;                                  % m/s (e.g., plastic)
delay = arrival_time - nanmin(arrival_time(:)); % seconds; min is zero
thickness = delay .* (c0 * c_lens) ./ (c_lens - c0);   % metres

% Make mm map, zero baseline, crop to 48 mm circular if you want
thickness_mm = thickness * 1e3;
thickness_mm = thickness_mm - nanmin(thickness_mm(:));
thickness_mm(~circ) = 0;                        % zero outside aperture

% Optional: crop to a specific circular diameter (e.g., 48 mm)
crop_mm = 48;
crop_px = round(crop_mm/(dx*1e3));
half   = crop_px/2;
x_idx  = round(cx - half + 1 : cx + half);
y_idx  = round(cy - half + 1 : cy + half);
lens_crop = thickness_mm(x_idx, y_idx);

[Xc,Yc] = meshgrid(1:crop_px, 1:crop_px);
rc = crop_px/2;
mask_circle = sqrt((Xc - rc).^2 + (Yc - rc).^2) <= rc;
lens_circ = lens_crop;
lens_circ(~mask_circle) = 0;

% Height-map export (16-bit PNG/TIFF not requiring Image Toolbox)
img16 = uint16( 65535 * (lens_circ / max(lens_circ(:) + eps)) );
imwrite(img16, 'lens_heightmap_48mm_circle.png');   % many tools accept PNG16

% 3D surface preview
[xmm, ymm] = meshgrid((1:crop_px)*dx*1e3);
figure; surf(xmm, ymm, lens_circ'); shading interp; colormap turbo; colorbar;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Lens thickness (mm)');
title('3D Surface — 48 mm circular lens (time-reversal derived)');
view(45,30); axis equal tight;

%% -------------------- Helper: Tukey window --------------------
function w = myTukeywin(N, r)
    if r <= 0, w = ones(N,1); return; end
    if r >= 1
        n = (0:N-1)'; w = 0.5*(1 - cos(2*pi*n/(N-1))); return; % Hann
    end
    w = zeros(N,1);
    L = floor(r*(N-1)/2);
    for n = 1:N
        if n <= L
            w(n) = 0.5*(1 + cos(pi*((2*n)/(r*(N-1)) - 1)));
        elseif n > L && n <= N - L
            w(n) = 1;
        else
            w(n) = 0.5*(1 + cos(pi*((2*n)/(r*(N-1)) - (2/r) + 1)));
        end
    end
end
