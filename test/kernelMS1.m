% Calculates the MS1 convolution kernel
function kernel = kernelMS1(deg, m)
  coeffs = corrCoeffsMS1(deg);
  kappa = [];
  for abcd = coeffs'
    kappa(end+1) = abcd(1) + abcd(2)/cube(abcd(3)-m);
  end
  kernel(m+1) = windowMS(0, 2); #center element
  for i = 1 : m
    x = i/(m+1);
    w = windowMS(x, 2);
    a = sin((0.5*deg+1)*pi*x)/((0.5*deg+1)*pi*x);
    for j = 1 : length(kappa)
      a = a + kappa(j)*x*sin(j*pi*x);
    end
    a = a*w;
    kernel(m+1-i) = a;
    kernel(m+1+i) = a;
  end
  norm = sum(kernel);
  kernel = kernel./norm;
end
