% Gaussian-like window function for the MS and MS1 kernels.
% The function reaches 0 at x=+1 and x=-1 (where the kernel ends);
% at these points also the 1st derivative is very close to zero.
function w = windowMS(x, alpha)
  w = exp(-alpha.*x.*x) + exp(-alpha.*(x+2).*(x+2)) + exp(-alpha.*(x-2).*(x-2)) - (2*exp(-alpha)+exp(-9*alpha));
end
