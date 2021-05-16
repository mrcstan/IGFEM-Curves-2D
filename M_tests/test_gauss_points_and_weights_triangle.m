path(path,'../M_FEM')
npts = 33;
% perform the integration int_0^1 int_0^{1-x} x^p y^q dy dx, 
% where p and q are non-negative integers on a triangular domain
% the analytical answer is p!(q+1)!/((1+q)(2+p+q)!)
maxDeg = 13;
for p = 1:maxDeg
    q = maxDeg - p;
    analyticalInt = factorial(p)*factorial(q+1)/((q+1)*factorial(2+p+q));


    func = @(x,y) x^p*y^q;

    gauss = gauss_points_and_weights(true,npts,2,'combined');
    gaussInt = 0;
    for i = 1:size(gauss,2)
        gaussInt = gaussInt + func(gauss(1,i),gauss(2,i))*gauss(3,i);
    end

    fprintf('absolute error = %g \n',abs(analyticalInt-gaussInt))
    fprintf('relative error = %g \n',abs(1-gaussInt/analyticalInt))
end