% Calculates the MS convolution kernel
function kernel = kernelMS(deg, m)
  coeffs = corrCoeffsMS(deg);
  kappa = []; #correction coefficients for the kernel
  for abcd = coeffs'
    kappa(end+1) = abcd(1) + abcd(2)/cube(abcd(3)-m);
  end
  if (rem(deg/2, 2) == 1) #degree 6, 10
    nuMinus2 = -1;
  else
    nuMinus2 = 0;
  end
  kernel(m+1) = windowMS(0, 4); #center element
  for i = 1 : m
    x = i/(m+1);
    w = windowMS(x, 4);
    a = sin((0.5*deg+2)*pi*x)/((0.5*deg+2)*pi*x);
    for j = 1 : length(kappa)
      a = a + kappa(j)*x*sin((2*j+nuMinus2)*pi*x);
    end
    a = a*w;
    kernel(m+1-i) = a;
    kernel(m+1+i) = a;
  end
  norm = sum(kernel);
  kernel = kernel./norm;
end
