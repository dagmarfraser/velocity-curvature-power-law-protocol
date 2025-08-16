function [dx, dy] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs)
%% Kinematic differentiation with multiple filter strategies
% Implementation for Fraser et al. (2024) Experimental Brain Research review
% 
% Correspondence: Dagmar Scott Fraser, d.s.fraser@bham.ac.uk
% Created: May 2024, Updated: June 2024
%
%% INPUTS
% x, y          - Trajectory coordinates (column vectors)
% filterType    - Differentiation method:
%   1: MATLAB diff (finite differences)
%   2: Butterworth low-pass → finite differences  
%   3: Finite differences → Butterworth low-pass
%   4: Savitzky-Golay smoothing differentiator as per Crenna et al. (2021)
%   5: Lacquaniti et al. (1983): exponential filter → Lagrange 5-point
%   6: Savitzky-Golay with temporal scaling
%   7: Savitzky-Golay with bandwidth equivalence as per Schafer 2011
% filterParams  - Filter parameters [order, frequency/length, zeroLag]:
%   Cases 1-3: [order, cutoff_Hz, zeroLag] where zeroLag: 0=filter, 1=filtfilt
%   Case 4:    [polynomial_order, window_length, ~]
%   Case 5:    [~, cutoff_Hz, ~] 
%   Case 6:    [polynomial_order, reference_window] temporal scaling from 100Hz reference of Case 4
%   Case 7:    [polynomial_order, target_cutoff_Hz] bandwidth matching Case 4 for =/= 100Hz
% fs            - Sampling frequency (Hz)
%
%% OUTPUTS  
% dx, dy        - N×4 matrices containing [position, velocity, acceleration, jerk]
%                 Position column contains smoothed coordinates for SG methods

dx = zeros(length(x), 4);
dy = zeros(length(y), 4);
dt = 1/fs;

switch filterType

    case 1 % Finite differences
        dx(:,1) = x;
        dy(:,1) = y;
        
        dx(2:end,2) = diff(x,1) * fs;           % velocity
        dy(2:end,2) = diff(y,1) * fs;
        
        dx(2:end-1,3) = diff(x,2) * fs^2;       % acceleration  
        dy(2:end-1,3) = diff(y,2) * fs^2;
        
        dx(2:end-2,4) = diff(x,3) * fs^3;       % jerk
        dy(2:end-2,4) = diff(y,3) * fs^3;

    case 2 % Butterworth filter → finite differences
        dx(:,1) = x;
        dy(:,1) = y;
        
        Fp = filterParams(2);
        N = filterParams(1);
        zeroLag = filterParams(3);
        
        [b,a] = butter(N, Fp/(fs/2));
        
        if zeroLag
            x_filt = filtfilt(b, a, x);         % zero-phase filtering
            y_filt = filtfilt(b, a, y);
        else
            x_filt = filter(b, a, x);           % causal filtering
            y_filt = filter(b, a, y);
        end
        
        dx(2:end,2) = diff(x_filt) * fs;        % velocity
        dy(2:end,2) = diff(y_filt) * fs;
        
        dx(2:end-1,3) = diff(dx(2:end,2)) * fs; % acceleration
        dy(2:end-1,3) = diff(dy(2:end,2)) * fs;
        
        dx(2:end-2,4) = diff(dx(2:end,2),2) * fs^2; % jerk
        dy(2:end-2,4) = diff(dy(2:end,2),2) * fs^2;

    case 3 % Finite differences → Butterworth filter
        dx(:,1) = x;
        dy(:,1) = y;
        
        Fp = filterParams(2);
        N = filterParams(1);
        zeroLag = filterParams(3);
        
        [b,a] = butter(N, Fp/(fs/2));
        
        if zeroLag
            dx(2:end,2) = filtfilt(b, a, diff(x) * fs);
            dy(2:end,2) = filtfilt(b, a, diff(y) * fs);
        else
            dx(2:end,2) = filter(b, a, diff(x) * fs);
            dy(2:end,2) = filter(b, a, diff(y) * fs);
        end
        
        dx(2:end-1,3) = diff(dx(2:end,2)) * fs;
        dy(2:end-1,3) = diff(dy(2:end,2)) * fs;
        
        dx(2:end-2,4) = diff(dx(2:end,2),2) * fs^2;
        dy(2:end-2,4) = diff(dy(2:end,2),2) * fs^2;

    case 4 % Savitzky-Golay smoothing differentiator
        order = filterParams(1);
        framelen = filterParams(2);
        
        [~,g] = sgolay(order, framelen);
        
        % Compute derivatives using SG coefficients scaled for sampling interval
        for p = 0:3
            scale_factor = factorial(p) / (-dt)^p;
            
            dx(:,p+1) = conv(x, scale_factor * g(:,p+1), 'same');
            dy(:,p+1) = conv(y, scale_factor * g(:,p+1), 'same');
        end

    case 5 % Lacquaniti et al. (1983) method
        % Double-sided exponential filter followed by Lagrange 5-point differentiation
        
        Fp = filterParams(2);
        alpha = 1 - exp(-2 * pi * Fp / fs);
        
        % Apply double-sided exponential filtering
        x_filt_fwd = exponentialFilter(x, alpha);
        y_filt_fwd = exponentialFilter(y, alpha);
        
        x_filt = flip(exponentialFilter(flip(x_filt_fwd), alpha));
        y_filt = flip(exponentialFilter(flip(y_filt_fwd), alpha));
        
        dx(:,1) = x_filt;                       % smoothed position
        dy(:,1) = y_filt;
        
        % Apply Lagrange 5-point differentiation
        dx(:,2) = lagrange5PointDiff(x_filt, dt); % velocity
        dy(:,2) = lagrange5PointDiff(y_filt, dt);
        
        dx(:,3) = lagrange5PointDiff(dx(:,2), dt); % acceleration
        dy(:,3) = lagrange5PointDiff(dy(:,2), dt);
        
        dx(:,4) = lagrange5PointDiff(dx(:,3), dt); % jerk
        dy(:,4) = lagrange5PointDiff(dy(:,3), dt);

    case 6 % Temporal-scaled Savitzky-Golay (Fraser et al.)
        % Maintains equivalent temporal window across sampling frequencies
        
        order = filterParams(1);
        reference_framelen = filterParams(2);
        reference_fs = 100;                     % Crenna et al. reference
        
        % Scale window to preserve temporal width
        temporal_width = reference_framelen / reference_fs;
        scaled_framelen = round(temporal_width * fs);
        
        % Ensure valid SG constraints (odd length, minimum order+1)
        min_window = order + 1;
        if scaled_framelen < min_window
            scaled_framelen = min_window;
        end
        
        % Enforce form 4k+1 for compatibility with legacy padding logic
        remainder = mod(scaled_framelen - 1, 4);
        scaled_framelen = scaled_framelen + (4 - remainder) * (remainder ~= 0);
        
        if mod(scaled_framelen, 2) == 0
            scaled_framelen = scaled_framelen + 1;
        end
        
        framelen = scaled_framelen;
        [~,g] = sgolay(order, framelen);
        
        for p = 0:3
            scale_factor = factorial(p) / (-dt)^p;
            
            dx(:,p+1) = conv(x, scale_factor * g(:,p+1), 'same');
            dy(:,p+1) = conv(y, scale_factor * g(:,p+1), 'same');
        end

    case 7 % Bandwidth-equivalent Savitzky-Golay (Schafer 2011)
        % Matches frequency response characteristics across sampling rates
        
        order = filterParams(1);
        target_cutoff_hz = filterParams(2);
        bandwidth_factor = 2.4;                 % Empirical constant for 4th order SG
        
        % Calculate window size for target bandwidth
        window_half_width = fs / (bandwidth_factor * target_cutoff_hz);
        framelen = 2 * round(window_half_width) + 1;
        
        % Apply SG constraints
        min_window = order + 1;
        if framelen < min_window
            framelen = min_window;
            if mod(framelen, 2) == 0
                framelen = framelen + 1;
            end
        end
        
        [~,g] = sgolay(order, framelen);
        
        for p = 0:3
            scale_factor = factorial(p) / (-dt)^p;
            
            dx(:,p+1) = conv(x, scale_factor * g(:,p+1), 'same');
            dy(:,p+1) = conv(y, scale_factor * g(:,p+1), 'same');
        end
end

end

function y_filtered = exponentialFilter(x, alpha)
% Single-sided exponential filter: y[n] = α·x[n] + (1-α)·y[n-1]

y_filtered = zeros(size(x));
y_filtered(1) = x(1);

for i = 2:length(x)
    y_filtered(i) = alpha * x(i) + (1 - alpha) * y_filtered(i-1);
end
end

function dx_dt = lagrange5PointDiff(x, dt)
% Lagrange 5-point differentiation: f'(x₀) = [-f(x₋₂) + 8f(x₋₁) - 8f(x₁) + f(x₂)] / (12h)

dx_dt = zeros(size(x));

% Edge cases using lower-order approximations
dx_dt(1) = (x(2) - x(1)) / dt;                          % forward difference
dx_dt(2) = (x(3) - x(1)) / (2*dt);                      % 3-point central

% Central 5-point formula for interior points
for i = 3:length(x)-2
    dx_dt(i) = (-x(i-2) + 8*x(i-1) - 8*x(i+1) + x(i+2)) / (12*dt);
end

% Edge cases
dx_dt(end-1) = (x(end) - x(end-2)) / (2*dt);           % 3-point central
dx_dt(end) = (x(end) - x(end-1)) / dt;                  % backward difference
end
