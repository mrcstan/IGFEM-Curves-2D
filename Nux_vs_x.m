close all
clear all
x = linspace(0,1,51);
gamsq = [25.68,83.86,174.2,296.5,450.9];
A = [7.630,2.058,0.901,0.487,0.297]*1e-3;
Nuinf = 4.364;
Nux = nan(size(x));
for i = 1:numel(Nux)
    Nux(i) = 1./(1/Nuinf - 0.5*sum(exp(-gamsq.*x(i))./(A.*gamsq.^2)));
end
plot(x,Nux)