% Octave implementation of smoothing by a modified sinc
% kernel (MS), as described in M. Schmid, D. Rath and U. Diebold,
% 'Why and how Savitzky-Golay filters should be replaced',
% ACS Measurement Science Au, 2022 
%
% The term 'degree' (variable 'deg' to avoid conflict with the
% Matlab 'degree' function) is defined in analogy to Savitzky-Golay
% (SG) filters; the MS filters have a similar frequency response as SG filters
% of the same deg (2, 4, ... 10).
%
% Note: This file contains more than one function, therefore it is not a
% function file but rather a script file.
% For actual usage, delete the test code at the very end and
% include it in your code with 'source("modifiedSincSmoother.m")'
% When using only one smooth type (smoothMS _or_ smoothMS1),
% this file may be converted to a single-function file by renaming
% according to the name of the main function and placing the other
% functions as inner functions. The, delete the following line.
%
% Copyright notice
% This code is licensed under GNU General Public License (GPLv3)
% When using and/or modifying this program for scientific work
% and the paper on it has been published, please cite the paper.
%
%  Author: Michael Schmid, IAP/TU Wien, 2021.
%          https://www.iap.tuwien.ac.at/www/surface/group/schmid


% The same with the shorter MS1 kernel
function smoothedData = smoothMS1(data, deg, m)
  if (nargin != 3)
    error("Usage: smoothMS1(dataRowVector, degree, m)");
  end
  if (columns(data) < 2 || rows(data) != 1)
    error("Less than two data points or not a row vector");
  end
  kernel = kernelMS1(deg, m);
  fitWeights = edgeWeights1(deg, m);
  extData = extendData(data, m, fitWeights);
  smoothedExtData = conv(extData, kernel, "same");
  smoothedData = smoothedExtData(m+1 : end-m);
end
