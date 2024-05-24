% Correction coeffficients for a flat passband of the MS1 kernel
function coeffs = corrCoeffsMS1(deg)
  switch (deg)
    case (2)
      coeffs = [];
    case (4)
      coeffs = [0.021944195, 0.050284006, 0.765625];
    case (6)
      coeffs = [0.001897730, 0.00847681, 1.2625;
                0.023064667, 0.13047926, 1.2265625];
    case(8)
      coeffs = [0.006590300, 0.05792946, 1.915625;
                0.002323448, 0.01029885, 2.2726562;
                0.021046653, 0.16646601, 1.98125];
    case (10)
      coeffs = [9.749618E-4, 0.00207429, 3.74375;
                0.008975366, 0.09902466, 2.707812;
                0.002419541, 0.01006486, 3.296875;
                0.019185117, 0.18953617, 2.784961];
    otherwise error("Invalid deg");
  end
end
