% f(x,d) = d^2 x - a
clear all
close all
path(path,'../../MatlabUsefulFunctions/presentation')
a = 1;
p = 2:1:8;
d = linspace(-1,1,200);

dsq = d.^2;

figure
hold on
np = numel(p);
colors = distinguishable_colors(np);
h = zeros(np,1);
lgd = cell(np,1);
for i = 1:np
    theta = (((dsq-a).^(p(i)+1)-(-a)^(p(i)+1))./(dsq*(p(i)+1))).^(1/p(i));
    h(i) = plot(d,theta,'linestyle','-','color',colors(i,:),'linewidth',2);
    lgd{i} = num2str(p(i));
end

legend(h,lgd)
set(gca,'fontsize',30)
xlabel('p','fontsize',30)
ylabel('\theta_p','fontsize',30)



