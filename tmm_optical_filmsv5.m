function [R, T, A, r, t, M_global, theta] = tmm_optical_filmsv5(lambda, theta0, pol, n, d)
% tmm_optical_filmsv5 calculates the optical properties of a multilayer thin film
% using the Transfer Matrix Method and plots the spectra.
%
% It supports:
%   (i)   A wavelength sweep mode (if lambda is a vector) to compute Reflection,
%         Transmission, and Absorption versus wavelength.
%   (ii)  A single wavelength mode (if lambda is scalar) to compute and display
%         the optical properties at that wavelength.
%
% Additionally, if run interactively (i.e. no output arguments), after plotting
% the spectra it will prompt you to compute the electric field intensity profile
% (calculated at a single user-specified wavelength).
%
% Usage examples:
%
%   1. Interactive mode (prompts for inputs and then for electric field profile):
%         >> tmm_optical_filmsv4
%
%   2. Single wavelength mode:
%         >> [R, T, A, r, t, M_global, theta] = tmm_optical_filmsv4(500, 0, 's', [1, 1.45, 2.0, 1], [100, 50]);
%
%   3. Wavelength sweep mode:
%         >> lambda = linspace(400,800,200);
%         >> [R, T, A] = tmm_optical_filmsv4(lambda, 0, 's', [1, 1.45, 2.0, 1], [100, 50]);
%         This will plot T (blue) and A (red) versus wavelength.
%
% Inputs:
%   lambda - scalar (single wavelength in nm) or vector (wavelength sweep in nm)
%   theta0 - incidence angle (in degrees)
%   pol    - polarization ('s' for s-polarization or 'p' for p-polarization)
%   n      - vector of refractive indices [n0, n1, ..., n_{M+1}], where n0 is the incident medium
%            and n_{M+1} is the substrate. (Length must equal number of layers + 2.)
%   d      - vector of layer thicknesses (in nm) for the thin-film layers (length = number of layers)
%
% Outputs (for single wavelength mode):
%   R        - Reflectance (|r|^2)
%   T        - Transmittance (with admittance ratio factor)
%   A        - Absorption (A = 1 - R - T)
%   r        - Complex reflection amplitude
%   t        - Complex transmission amplitude
%   M_global - Global transfer matrix for the multilayer structure
%   theta    - Vector of propagation angles in each medium (in radians)
%
% If fewer than 5 input arguments are provided, prompt for them interactively.
if nargin < 5
    lambda = input('Enter wavelength (in nm) or a vector (e.g., linspace(400,800,200)): ');
    theta0 = input('Enter incidence angle (in degrees): ');
    pol = input('Enter polarization (''s'' or ''p''): ', 's');
    n = input('Enter refractive indices as a vector (e.g., [1, 1.45, 2.0, 1]): ');
    d = input('Enter layer thicknesses as a vector (e.g., [100, 50]): ');
end

% Validate that the refractive index vector length equals number of layers + 2.
if length(n) ~= length(d) + 2
    error('Refractive index vector must have length equal to (number of layers + 2).');
end

%---------------------------------------------------------------------------
% Check if lambda is a scalar (single wavelength) or a vector (wavelength sweep).
if numel(lambda) == 1
    % Single wavelength computation.
    [R, T, A, r, t, M_global, theta] = computeOpticalProps(lambda, theta0, pol, n, d);
    
    % Plot a bar graph for the optical properties.
    figure;
    bar([R, T, A]);
    set(gca, 'XTickLabel', {'Reflectance (R)', 'Transmittance (T)', 'Absorption (A)'});
    title(sprintf('Optical properties at \\lambda = %.2f nm', lambda));
    ylabel('Value');
else
    % Wavelength sweep: initialize arrays.
    numWave = length(lambda);
    R_array = zeros(1, numWave);
    T_array = zeros(1, numWave);
    A_array = zeros(1, numWave);
    
    % Loop over all wavelengths in the sweep.
    for k = 1:numWave
        [R_temp, T_temp, A_temp] = computeOpticalProps(lambda(k), theta0, pol, n, d);
        R_array(k) = R_temp;
        T_array(k) = T_temp;
        A_array(k) = A_temp;
    end
    
    % Plot the Transmittance and Absorption spectra.
    figure;
    plot(lambda, T_array, 'b-', 'LineWidth', 2);
    hold on;
    plot(lambda, A_array, 'r-', 'LineWidth', 2);
    xlabel('Wavelength (nm)');
    ylabel('Value');
    legend('Transmittance (T)', 'Absorption (A)');
    title('Optical Spectra: Transmittance and Absorption');
    grid on;
    
    % Return computed arrays.
    R = R_array;
    T = T_array;
    A = A_array;
    r = [];
    t = [];
    M_global = [];
    theta = [];
end

%---------------------------------------------------------------------------
% If the function is run interactively (no output arguments), ask if the user
% wants to compute the electric field intensity profile at a single wavelength.
if nargout == 0
    prompt = 'Do you want to compute the electric field intensity profile? (y/n): ';
    computeField = input(prompt, 's');
    if lower(computeField) == 'y'
        field_lambda = input('Enter wavelength (in nm) for electric field calculation: ');
        computeElectricFieldProfile(field_lambda, theta0, pol, n, d);
    end
end

end

%==========================================================================
function [R, T, A, r, t, M_global, theta] = computeOpticalProps(lambda, theta0, pol, n, d)
% computeOpticalProps computes the optical properties at a single wavelength.
%
% Inputs:
%   lambda - scalar wavelength (in nm)
%   theta0, pol, n, d - same as in the main function.
%
% Outputs:
%   R, T, A - Reflectance, Transmittance, and Absorption.
%   r, t, M_global, theta - as defined in the main function.
%
% Step 1. Compute propagation angles using Snell's law.
theta0_rad = theta0 * pi/180;   % Convert incidence angle to radians
M_layers = length(d);           % Number of thin-film layers
N_total = length(n);            % Total media count = M_layers + 2

theta = zeros(1, N_total);
theta(1) = theta0_rad;
for j = 2:N_total
    sin_theta_j = n(1) * sin(theta0_rad) / n(j);
    theta(j) = asin(sin_theta_j);
end

% Step 2. Compute admittance for each medium.
% For s-polarization: s = n*cos(theta)
% For p-polarization: s = cos(theta)/n
s = zeros(1, N_total);
switch lower(pol)
    case 's'
        for j = 1:N_total
            s(j) = n(j) * cos(theta(j));
        end
    case 'p'
        for j = 1:N_total
            s(j) = cos(theta(j)) / n(j);
        end
    otherwise
        error('Polarization must be either ''s'' or ''p''.');
end

% Step 3. Build the global transfer matrix.
M_global = eye(2);  % Start with the identity matrix.
% Interface from incident medium (n0) to first layer (n1)
M_global = M_global * InterfaceMatrix(s(1), s(2));

for j = 1:M_layers
    % Compute phase shift in layer j (medium index = j+1)
    phi = (2*pi/lambda) * n(j+1) * d(j) * cos(theta(j+1));
    % Propagation (phase) matrix for layer j:
    P = [exp(-1i*phi), 0; 0, exp(1i*phi)];
    M_global = M_global * P;
    % Interface from layer j (medium j+1) to next medium (j+2)
    M_global = M_global * InterfaceMatrix(s(j+1), s(j+2));
end

% Step 4. Apply boundary conditions.
% In the substrate (last medium), assume no backward (reflected) wave.
t = 1 / M_global(1,1);
r = M_global(2,1) / M_global(1,1);

% Step 5. Compute Reflectance, Transmittance, and Absorption.
R = abs(r)^2;
T = (real(s(end)) / real(s(1))) * abs(t)^2;
A = 1 - R - T;

end

%==========================================================================
function computeElectricFieldProfile(lambda_field, theta0, pol, n, d)
% computeElectricFieldProfile computes and plots the electric field intensity
% profile (|E|^2) inside the multilayer stack for a single wavelength.
%
% Inputs:
%   lambda_field - scalar wavelength (in nm) for field calculation.
%   theta0, pol, n, d - same as in the main function.
%
% The routine computes the field inside each layer using the local
% forward and backward wave amplitudes.
%
% Step 1. Compute propagation angles and admittances.
theta0_rad = theta0 * pi/180;
M_layers = length(d);
N_total = length(n);

theta = zeros(1, N_total);
theta(1) = theta0_rad;
for j = 2:N_total
    sin_theta_j = n(1) * sin(theta0_rad) / n(j);
    theta(j) = asin(sin_theta_j);
end

s = zeros(1, N_total);
switch lower(pol)
    case 's'
        for j = 1:N_total
            s(j) = n(j) * cos(theta(j));
        end
    case 'p'
        for j = 1:N_total
            s(j) = cos(theta(j)) / n(j);
        end
    otherwise
        error('Polarization must be either ''s'' or ''p''.');
end

% Step 2. Compute the global reflection coefficient r at lambda_field.
[~, ~, ~, r, ~, ~, ~] = computeOpticalProps(lambda_field, theta0, pol, n, d);

% Step 3. Compute the field inside each layer.
% We assume the incident field amplitude is 1; the field vector in the incident medium is [1; r].
% We compute the cumulative transfer matrix to the beginning of each layer.
M_current = InterfaceMatrix(s(1), s(2)); % Field at beginning of layer 1.
x_total = [];   % To store the spatial coordinate (in nm)
E_total = [];   % To store the complex electric field
x_offset = 0;   % Cumulative thickness offset

for j = 1:M_layers
    % The current layer corresponds to medium index = j+1.
    % Compute the amplitude vector in this layer:
    F_layer = M_current * [1; r];   % [A_j; B_j]
    
    % Compute the x-component of the wavevector in layer j:
    k_j = (2*pi/lambda_field) * n(j+1) * cos(theta(j+1));
    
    % Create a spatial grid within layer j:
    x_j = linspace(0, d(j), 200); % 200 points within this layer
    % Electric field in layer j:
    E_j = F_layer(1) * exp(1i * k_j * x_j) + F_layer(2) * exp(-1i * k_j * x_j);
    
    % Append the positions (shifted by the cumulative offset) and fields.
    x_total = [x_total, x_offset + x_j];
    E_total = [E_total, E_j];
    
    % Update the cumulative offset.
    x_offset = x_offset + d(j);
    
    % Update the transfer matrix to the beginning of the next layer:
    % Full propagation through layer j:
    P_full = [exp(-1i * k_j * d(j)), 0; 0, exp(1i * k_j * d(j))];
    M_current = M_current * P_full * InterfaceMatrix(s(j+1), s(j+2));
end

% Normalize intensity relative to the incident intensity |1|^2.
I_total = abs(E_total).^2;

% Plot the electric field intensity profile.
figure;
plot(x_total, I_total, 'm-', 'LineWidth', 2);
xlabel('Position inside stack (nm)');
ylabel('Normalized Electric Field Intensity |E|^2');
title(sprintf('Electric Field Intensity Profile at \\lambda = %.2f nm', lambda_field));
grid on;

end

%==========================================================================
function T_mat = InterfaceMatrix(s1, s2)
% InterfaceMatrix computes the interface matrix between two media.
%
% The interface matrix is given by:
%    T_mat = 1/(2*s1) * [ s1+s2,  s1-s2;
%                         s1-s2,  s1+s2 ]
T_mat = 1/(2*s1) * [ s1+s2, s1-s2; s1-s2, s1+s2 ];
end
