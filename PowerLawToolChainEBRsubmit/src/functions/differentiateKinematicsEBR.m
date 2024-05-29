function [dx, dy] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs)
%% demo code for Exp Brain Res review paper of Fraser et al., 2024
% Created May 2024
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk
%
%% inputs 
% x and y coordinates
% filterType - chooses filter and differentiation
%   1 MATLAB diff, scaled by the sample rate to give an approximation of the 
%   2 Nth Order Fp Hz Low pass filter, filtfilt for zero lag followed by Finite Differences 
%   3 Finite Differences followed by Nth Order Fp Hz Low pass filter filtfilt for zero lag 
%   4 Savitzky-Golay smoothing differential filter.
% filterParams - [filter order, Fc Low Pass Cutt off for Butterworth OR
% width for S-G filter zeroLag] 
% - where zeroLag = 1 filtfilt, zeroLag = 0 just employ filter
% fs - sampling frequency of the data
%% outputs
% dx and dy - N x 4, where N = length(input data)
% rows with original data (smoothed only in the case of S-G)
% velocity
% acceleration
% jerk
%   
dx = zeros(length(x),4);
dy = zeros(length(y),4);
dt = 1/fs;

switch filterType

    case 1 % Finite Differences Differentation

        dx = x; % not smoothed!
        dy = y;

        dx(2:end,2) = diff(x,1) * fs; % raw diff -  scaled by the fs
        dy(2:end,2) = diff(y,1) * fs; % velocity

        dx(2:end-1,3) = diff(x,2) * fs * fs; % raw diff -  scaled by the fs^2
        dy(2:end-1,3) = diff(y,2) * fs * fs; % acceleration

        dx(2:end-2,4) = diff(x,3) * fs * fs * fs; % raw diff -  scaled by the fs^3
        dy(2:end-2,4) = diff(y,3) * fs * fs * fs; % jerk

    case 2 % Nth Order Fp Hz Low pass filter, filtfilt for zero lag followed by Finite Differences 

        dx = x; % not smoothed!
        dy = y;
        
        Fp = filterParams(2);
        N = filterParams(1);
        zeroLag = filterParams(3);
        %% make an Nth order Butterworth zero lag low pass filter with corner frequency Fp
        fc = Fp;
        bOrder = N; % this is filter order, will be 2*bOrder when used with filtfilt below
        [b,a] = butter(bOrder,fc/(fs/2));

        % use zerophase digital filtering i.e. filtfilt
        % diff the output of the low pass filter
        if ~zeroLag
            dx(2:end,2) = diff(filter(b,a, (x))) * fs; % raw diff -  scaled by the fs, and now diminshed by the filter
            dy(2:end,2) = diff(filter(b,a, (y))) * fs; % velocity
        else
            dx(2:end,2) = diff(filtfilt(b,a, (x))) * fs; % raw diff -  scaled by the fs, and now diminshed by the filter
            dy(2:end,2) = diff(filtfilt(b,a, (y))) * fs; % velocity
        end

        dx(2:end-1,3) = diff(dx(2:end,2),1) * fs; % raw diff -  scaled by the fs
        dy(2:end-1,3) = diff(dy(2:end,2),1) * fs; % acceleration

        dx(2:end-2,4) = diff(dx(2:end,2),2) * fs * fs; % raw diff - scaled by the fs
        dy(2:end-2,4) = diff(dy(2:end,2),2) * fs * fs; % jerk

    case 3 % Finite Differences with Nth Order Fp Hz Low pass filter

        dx = x; % not smoothed!
        dy = y;

        Fp = filterParams(2);
        N = filterParams(1);
        zeroLag = filterParams(3);

        %% make an Nth order Butterworth zero lag low pass filter with corner frequency Fp
        fc = Fp;
        bOrder = N; % this is filter order, will be 2*bOrder when used with filtfilt below
        [b,a] = butter(bOrder,fc/(fs/2)); % we pass the order and the cut off freq / Nyquist Frequency (i.e. half sampling rate)

        % use zerophase digital filtering i.e. filtfilt
        % low pass filter the output of the diff
        if ~zeroLag
            dx(2:end,2) = filter(b,a, (diff(x,1)*fs)); % raw diff -  scaled by the fs, and now diminshed by the filter
            dy(2:end,2) = filter(b,a, (diff(y,1)*fs));
        else
            dx(2:end,2) = filtfilt(b,a, (diff(x,1)*fs)); % raw diff -  scaled by the fs, and now diminshed by the filter
            dy(2:end,2) = filtfilt(b,a, (diff(y,1)*fs));
        end
        dx(2:end-1,3) = diff(dx(2:end,2),1) * fs; % raw diff -  scaled by the fs
        dy(2:end-1,3) = diff(dy(2:end,2),1) * fs;

        dx(2:end-2,4) = diff(dx(2:end,2),2) * fs * fs; % raw diff - scaled by the fs
        dy(2:end-2,4) = diff(dy(2:end,2),2) * fs * fs;

    case 4 % Savitzky-Golay smooth differentiation
        % adapted from https://uk.mathworks.com/help/signal/ref/sgolay.html
        % example

        order = filterParams(1);
        framelen = filterParams(2);
        padding = (framelen-1)/2;

        [b,g] = sgolay(order, framelen);

        for p = 0:3 % smoothed displacement, velocity, acceleration, jerk.

            dxElement = [ conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same')];
            dyElement = [ conv(y, factorial(p)/(-dt)^p * g(:,p+1), 'same')];
            % Divide the columns by powers of dt to scale the derivatives correctly.

            if p > 0 
                % deal with edge effects / returned smaller derivative
                % arrays
                dxPadded = padarray(dxElement * p, padding/2, 0, 'pre');
                dyPadded = padarray(dyElement * p, padding/2, 0, 'pre');

                dx(:,p+1) = dxPadded(2:length(x)+1);
                dy(:,p+1) = dyPadded(2:length(x)+1);

            end
        end

    case 5

        disp('NOT IMPLEMNTED')

end