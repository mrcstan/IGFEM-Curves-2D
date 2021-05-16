A = rand(3);
B = ones(5);

for i = 1:5
    try
       C = A*B;
    catch
       %{
       % Give more information for mismatch.
       if (strcmp(err.identifier,'MATLAB:catenate:dimensionMismatch'))

          msg = ['Dimension mismatch occurred: First argument has ', ...
                num2str(size(A,2)), ' columns while second has ', ...
                num2str(size(B,2)), ' columns.'];
          error('MATLAB:myCode:dimensions', msg);

       % Display any other errors as usual.
       else
          disp('hahaha')
       end
        %}
        disp('hahaha')
    end  % end try/catch
end