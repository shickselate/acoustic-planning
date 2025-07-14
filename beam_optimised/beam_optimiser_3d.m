clear; clc;

% --- Parameters ---
Nx = 44; Ny = 44; Nz = 96; dx = 1e-3;
f0 = 500e3;
c = 1500;
focus = [Nx/2, Ny/2, 70];
iterations = 10;
amplitude = 1e5;
window_r = 0.5;
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

% --- Initial phase map ---
phases = zeros(Nx, Ny);

% --- Target mask and FSLR log ---
target_mask = false(Nx, Ny, Nz);
target_mask(focus(1)-1:focus(1)+1, focus(2)-1:focus(2)+1, focus(3)-1:focus(3)+1) = true;
target_idx = find(target_mask);
non_target_idx = find(~target_mask);
FSLR_log = zeros(iterations, 1);

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
        total_phase = 2 * pi * f0 * delay + phases(ix, iy);
        amp = apod(ix, iy) * amplitude;
        source.p(n, :) = amp * sin(2 * pi * f0 * t - total_phase);
    end

    % Sensor everywhere
    sensor.mask = ones(Nx, Ny, Nz);

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
save('best_beam_iterationX.mat', 'best_phases', 'best_p_map', ...
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
