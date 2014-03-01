function gaussian = make2DGaussianCentered(w, h, x0, y0, sigma, gaussAmp)
%
% function gaussian = make2DGaussian(w, h, x0, y0, sigma, gaussAmp)

% make sure width and height are even
if mod(w, 2)
    w = w - 1 ; 
    fprintf('make2DGaussianCentered: adjusting width, w, to be even ...\n')
end
if mod(h, 2)
    h = h - 1 ; 
    fprintf('make2DGaussianCentered: adjusting height, h, to be even ...\n')
end

[x y] = meshgrid(-w/2+1:w/2, -h/2+1:h/2);

gaussian = gaussAmp*exp(-(((x-x0)/(2*sigma)).^2)-(((y-y0)/(2*sigma)).^2));
