function [out1,out2] = output_tilde()
if (nargout == 2)
    out1 = 10;
    out2 = 20;
    disp('two output')
elseif (nargout == 1)
    out1 = 1;
    disp('one output')
end
end