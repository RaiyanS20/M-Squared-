% --- Function Definition Section ---

function [wx_physical, wy_physical, centroid_x_physical, centroid_y_physical] = calculateSecondMomentWidth(imageData, pixelSizeX, pixelSizeY, backgroundLevel)
% Calculates the second moment beam widths (wx, wy) and centroid (x, y)
% from a 2D image matrix according to ISO 11146 standard.
%
% Inputs:
%   imageData     : 2D matrix representing the beam intensity profile.
%                   Should be double precision (use im2double if needed).
%   pixelSizeX    : Physical size of one pixel in the horizontal direction
%                   (e.g., in mm or um).
%   pixelSizeY    : Physical size of one pixel in the vertical direction
%                   (e.g., in mm or um).
%   backgroundLevel : Estimated background level to subtract. Can be:
%                     - A scalar value (e.g., measured from dark frame).
%                     - 'auto' to estimate from image corners.
%                     - 0 or 'none' for no subtraction.
%                     Defaults to 'none'.
%
% Outputs:
%   wx_physical   : Second moment beam radius along x-axis (in same units as pixelSizeX).
%   wy_physical   : Second moment beam radius along y-axis (in same units as pixelSizeY).
%   centroid_x_physical: x-coordinate of the beam centroid (in same units as pixelSizeX).
%   centroid_y_physical: y-coordinate of the beam centroid (in same units as pixelSizeY).
%
% Note: The origin (0,0) of the physical coordinate system corresponds to
%       the center of the top-left pixel (index 1,1).

if nargin < 4 || isempty(backgroundLevel) || strcmpi(backgroundLevel, 'none')
    backgroundLevel = 0;
    performSubtraction = false;
else
    performSubtraction = true;
end

if ~isa(imageData, 'double')
    warning('Image data is not double precision. Converting using im2double.');
    imageData = im2double(imageData);
end

[numRows, numCols] = size(imageData);

% --- Background Subtraction ---
if performSubtraction
    if ischar(backgroundLevel) && strcmpi(backgroundLevel, 'auto')
        % Estimate background from corners (e.g., 5% of pixels)
        cornerSize = ceil(min(numRows, numCols) * 0.05);
        corners = [imageData(1:cornerSize, 1:cornerSize); ...
                   imageData(1:cornerSize, end-cornerSize+1:end); ...
                   imageData(end-cornerSize+1:end, 1:cornerSize); ...
                   imageData(end-cornerSize+1:end, end-cornerSize+1:end)];
        bg = mean(corners(:));
        fprintf('Auto-estimated background level: %.4f\n', bg);
    elseif isnumeric(backgroundLevel) && isscalar(backgroundLevel)
        bg = backgroundLevel;
        fprintf('Using specified background level: %.4f\n', bg);
    else
        warning('Invalid backgroundLevel specified. No subtraction performed.');
        bg = 0;
    end

    imageData = imageData - bg;
    % Clamp negative values to zero after subtraction
    imageData(imageData < 0) = 0;
else
     fprintf('No background subtraction performed.\n');
end

% --- Calculate Total Power (Zeroth Moment) ---
totalPower = sum(imageData(:));

if totalPower <= 0
    warning('Total power is zero or negative after background subtraction. Cannot calculate width.');
    wx_physical = NaN;
    wy_physical = NaN;
    centroid_x_physical = NaN;
    centroid_y_physical = NaN;
    return;
end

% --- Create Coordinate Grids (using pixel indices first) ---
% Note: MATLAB indices start from 1.
% meshgrid provides matrix-based coordinates
[X_pix, Y_pix] = meshgrid(1:numCols, 1:numRows);

% --- Calculate Centroid (First Moments) in Pixels ---
x_centroid_pix = sum(sum(X_pix .* imageData)) / totalPower;
y_centroid_pix = sum(sum(Y_pix .* imageData)) / totalPower;

% --- Calculate Variances (Second Central Moments) in Pixels ---
x_dev_sq = (X_pix - x_centroid_pix).^2; % Squared deviation from centroid x
y_dev_sq = (Y_pix - y_centroid_pix).^2; % Squared deviation from centroid y

sigma_x_sq_pix = sum(sum(x_dev_sq .* imageData)) / totalPower;
sigma_y_sq_pix = sum(sum(y_dev_sq .* imageData)) / totalPower;

% --- Calculate Second Moment Widths in Pixels ---
% According to ISO 11146, radius w = 2*sigma
wx_pix = 2 * sqrt(sigma_x_sq_pix);
wy_pix = 2 * sqrt(sigma_y_sq_pix);

% --- Convert to Physical Units ---
% Adjust centroid origin: subtract 0.5 because pixel (1,1) center is at (0.5, 0.5)
% in a system where the top-left corner is (0,0). Then multiply by pixel size.
centroid_x_physical = (x_centroid_pix - 0.5) * pixelSizeX;
centroid_y_physical = (y_centroid_pix - 0.5) * pixelSizeY;

wx_physical = wx_pix * pixelSizeX;
wy_physical = wy_pix * pixelSizeY;

fprintf('Centroid (X, Y): (%.2f, %.2f) pixels -> (%.4f, %.4f) physical units\n', ...
        x_centroid_pix, y_centroid_pix, centroid_x_physical, centroid_y_physical);
fprintf('Second Moment Width (Wx, Wy): (%.2f, %.2f) pixels -> (%.4f, %.4f) physical units\n', ...
        wx_pix, wy_pix, wx_physical, wy_physical);

end % function calculateSecondMomentWidth

% --- End of Function Definition Section ---