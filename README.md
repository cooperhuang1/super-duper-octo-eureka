clear; clc; close all;

L = 50; J = 1; kB = 1; Tmin = 1.0; Tmax = 4.0; Tn = 50; sweeps = 100;
lat = ones(L, L);
T = linspace(Tmin, Tmax, Tn);
M = zeros(1, Tn);

E = @(lat, i, j) -J * lat(i, j) * (...
    lat(mod(i-2, L)+1, j) + ...
    lat(mod(i, L)+1, j) + ...
    lat(i, mod(j-2, L)+1) + ...
    lat(i, mod(j, L)+1));

for t = 1:length(T)
    beta = 1 / (kB * T(t));
    for sweep = 1:sweeps
        for step = 1:L^2
            i = randi(L);
            j = randi(L);
            dE = -2 * E(lat, i, j);
            if dE <= 0 || rand < exp(-beta * dE)
                lat(i, j) = -lat(i, j);
            end
        end
    end
    M(t) = mean(lat(:));
    figure(1);
    imagesc(lat);
    axis equal tight;
    title(sprintf('T = %.2f', T(t)));
    colorbar;
    pause(0.1);
end

figure(2);
plot(T, M, 'r-o', 'LineWidth', 2);
xlabel('T');
ylabel('\langle M \rangle');
title('T vs \langle M \rangle');
grid on;
ylim([-1, 1]);
