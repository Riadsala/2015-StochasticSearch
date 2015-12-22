
function p = logisticReg_model(beta, x, y)

 
switch beta
     case 1.6
          b1 = -0.63;
          bx = -0.06;
          by = -0.05;
      case 1.65   
          b1 = -0.63 + 1.32;
          bx = -0.06 + 0.006;
          by = -0.05 - 0.05;
     case 1.7
           b1 = -0.63 + 2.20;
          bx = -0.06 - 0.001;
          by = -0.05 - 0.035;
end
 lm = b1 + bx.*x.^2 + by.*y.^2;

p = 1./(1+exp(-lm));

end