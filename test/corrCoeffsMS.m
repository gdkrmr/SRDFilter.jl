% Correction coeffficients for a flat passband of the MS kernel
function coeffs = corrCoeffsMS(deg)
  switch (deg)
    case (2)
      coeffs = [];
    case (4)
      coeffs = [];
    case (6)
      coeffs = [0.001717576, 0.02437382, 1.64375];
    case(8)
      coeffs = [0.0043993373, 0.088211164, 2.359375;
                0.006146815, 0.024715371, 3.6359375];
    case (10)
      coeffs = [0.0011840032, 0.04219344, 2.746875;
                0.0036718843, 0.12780383, 2.7703125];
    otherwise error("Invalid deg");
  end
end
