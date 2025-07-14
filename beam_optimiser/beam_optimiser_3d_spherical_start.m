clear; clc;

% --- Parameters ---
Nx = 44; Ny = 44; Nz = 128; dx = 1e-3;
f0 = 500e3;
c = 1500;
focus = [Nx/2, Ny/2, 30];
iterations = 2;
amplitude = 1e5;
window_r = 0.7;
learning_rate = 0.1;


% --- k-Wave setup ---
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);
medium.sound_speed = c;
medium.density = 1000;
kgrid.makeTime(c, [], 60e-6);
t = kgrid.t_array;
Nt = length(t);

% --- Create source mask ---
source_mask = zeros(Nx, Ny, Nz);
source_mask(:, :, 1) = 1;
source.p_mask = source_mask;
active = find(source_mask);
num_sources = length(active);

% --- Apodisation ---
wx = myTukeywin(Nx, window_r);
wy = myTukeywin(Ny, window_r);
apod = wx * wy';


% --- Smarter spherical phase initialisation ---
[X, Y] = meshgrid(1:Nx, 1:Ny);
Z = ones(size(X));
dists = sqrt((X - focus(1)).^2 + (Y - focus(2)).^2 + (focus(3)).^2);  % approx. to focus
delays = dists * dx / c;
phases = 2 * pi * f0 * delays;
phases = mod(phases, 2*pi);  % wrap phase to [0, 2pi]

% --- Visualise Initial Spherical Wavefront Phase Map ---
figure;
imagesc(phases);  % unwrapped phase
axis image; colormap turbo; colorbar;
title('Initial Phase Map (Unwrapped)');
xlabel('X'); ylabel('Y');

% Also visualise as a wrapped phase map (0 to 2π)
figure;
imagesc(mod(phases, 2*pi));
axis image; colormap hsv; colorbar;
title('Initial Phase Map (Wrapped to [0, 2π])');
xlabel('X'); ylabel('Y');

% Optional: 3D surface plot
[Xp, Yp] = meshgrid(1:Nx, 1:Ny);
figure;
surf(Xp, Yp, phases');
shading interp; colormap turbo; colorbar;
title('3D Surface of Initial Phase');
xlabel('X'); ylabel('Y'); zlabel('Phase (radians)');
view(45, 30); axis tight;



% --- Target mask and FSLR log ---
target_mask = false(Nx, Ny, Nz);
target_mask(focus(1)-1:focus(1)+1, focus(2)-1:focus(2)+1, focus(3)-1:focus(3)+1) = true;
target_idx = find(target_mask);
non_target_idx = find(~target_mask);
FSLR_log = zeros(iterations, 1);

% --- Sensor everywhere ---
sensor.mask = ones(Nx, Ny, Nz);

% --- Iterative optimisation ---
for iter = 1:iterations
    if iter == 1
        best_phases = phases;
        best_p_map = zeros(Nx, Ny, Nz);
        best_iter = 1;
    end

    % Build time-domain signals
    source.p = zeros(num_sources, Nt);
    for n = 1:num_sources
        [ix, iy, iz] = ind2sub(size(source_mask), active(n));
        delay = (sqrt((ix - focus(1))^2 + (iy - focus(2))^2 + (iz - focus(3))^2) * dx) / c;
        total_phase = 2 * pi * f0 * delay - phases(ix, iy);  % phase subtracted
        amp = apod(ix, iy) * amplitude;
        source.p(n, :) = amp * sin(2 * pi * f0 * t - total_phase);
    end

    % Run simulation
    sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, ...
        'DataCast', 'single', 'PMLInside', false, 'PMLSize', 10, ...
        'BinaryName', 'kspaceFirstOrder-CUDA.exe', 'PlotSim', false);

    % Analyse pressure
    p_max = max(sensor_data, [], 2);
    p_map = reshape(p_max, Nx, Ny, Nz);

    % Compute FSLR
    P_focus = mean(p_map(target_idx));
    P_sidelobe = mean(p_map(non_target_idx));
    FSLR_log(iter) = 20 * log10(P_focus / P_sidelobe);
    fprintf('Iter %d: FSLR = %.2f dB\n', iter, FSLR_log(iter));

    % Plot axial slice
    figure(10); clf;
    imagesc(squeeze(p_map(:, :, focus(3))));
    axis image; colormap hot; colorbar;
    title(sprintf('Axial Slice, Iter %d — FSLR = %.2f dB', iter, FSLR_log(iter)));
    drawnow;

    % Save best result
    if iter == 1 || FSLR_log(iter) > max(FSLR_log(1:iter-1))
        best_phases = phases;
        best_p_map = p_map;
        best_iter = iter;
    end

    % Update phase map
    grad = zeros(size(phases));
    region = p_map(focus(1)-1:focus(1)+1, focus(2)-1:focus(2)+1, focus(3)-1:focus(3)+1);
    contribution = mean(region(:));
    for n = 1:num_sources
        [ix, iy, iz] = ind2sub(size(source_mask), active(n));
        grad(ix, iy) = grad(ix, iy) + contribution;
    end
    grad = grad / max(grad(:));
    phases = phases + learning_rate * pi * (grad - 0.5);
end

% --- Save final result ---
save('best_beam_spherical_start.mat', 'best_phases', 'best_p_map', ...
    'focus', 'dx', 'Nx', 'Ny', 'Nz');
disp('Best beam and phase map saved.');

% --- Final plots ---
figure;
plot(1:iterations, FSLR_log, '-o', 'LineWidth', 2);
xlabel('Iteration'); ylabel('FSLR (dB)');
title('Focus-to-Sidelobe Ratio Across Iterations');
grid on;

figure;
imagesc(best_phases); axis image; colormap jet; colorbar;
title(sprintf('Best Phase Map (Iteration %d)', best_iter));

figure;
imagesc(squeeze(best_p_map(:, :, focus(3))));
axis image; colormap hot; colorbar;
title(sprintf('Best Axial Slice at Focus (Iter %d)', best_iter));

% --- Sagittal slice through the focus (Y-Z view) ---
figure;
imagesc(squeeze(best_p_map(focus(1), :, :))');
axis image; colormap hot; colorbar;
xlabel('Y'); ylabel('Z');
title(sprintf('Sagittal Slice at Focus (Iter %d)', best_iter));

z_profile = squeeze(best_p_map(focus(1), focus(2), :));
figure;
plot((1:Nz) * dx * 1e3, z_profile);
xlabel('Depth (mm)'); ylabel('Max Pressure');
title('Axial Pressure Profile');
grid on;


% % --- Sagittal slice through the focus (Y-Z view), log scale ---
% slice = squeeze(best_p_map(focus(1), :, :))';
% log_slice = log10(abs(slice) + 1e-6);  % Add epsilon to avoid log(0)
% 
% figure;
% imagesc(log_slice);
% axis image; colormap hot; colorbar;
% xlabel('Y'); ylabel('Z');
% title(sprintf('Sagittal Slice at Focus (Log Scale, Iter %d)', best_iter));
% caxis([8 15]);  % Adjust as needed to visualise contrast
% 
% 
% 
% % Beam profile along Z
% figure;
% line_profile = squeeze(best_p_map(focus(1), focus(2), :));
% plot((1:Nz)*dx*1e3, line_profile);
% xlabel('Z (mm)'); ylabel('Peak Pressure');
% title('Axial Pressure Profile');



% % Define which Z slice to view (e.g. central sagittal slice)
% x_idx = round(Nx / 2);
% 
% p = reshape(sensor_data, Nx, Ny, Nz, Nt);
% 
% % Loop over time
% for t_idx = 1:10:Nt  % Every 10th frame for speed
%     slice = squeeze(p(x_idx, :, :, t_idx))';  % Sagittal slice (Y-Z plane)
% 
%     imagesc(slice);
%     axis image;
%     caxis([-1 1] * max(abs(p(:))));  % consistent colour scale
%     colormap hot;
%     colorbar;
%     title(sprintf('Sagittal Pressure Field at t = %.2f µs', t(t_idx) * 1e6));
%     xlabel('Y'); ylabel('Z');
% 
%     drawnow;
%     pause(0.05);
% end



% --- Helper function ---
function w = myTukeywin(N, r)
    if r <= 0, w = ones(N,1); return; end
    if r >= 1, w = hann(N); return; end
    w = zeros(N,1); L = floor(r*(N-1)/2);
    for n = 1:N
        if n <= L
            w(n) = 0.5 * (1 + cos(pi*((2*n)/(r*(N-1)) - 1)));
        elseif n <= N - L
            w(n) = 1;
        else
            w(n) = 0.5 * (1 + cos(pi*((2*n)/(r*(N-1)) - (2/r) + 1)));
        end
    end
end
