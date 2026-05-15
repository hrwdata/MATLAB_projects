function value = evaluate_symbolic_at_point(expr, vars, x)
%EVALUATE_SYMBOLIC_AT_POINT Evaluate a symbolic expression at a numeric point.

vars = vars(:).';
x = x(:).';

if numel(vars) ~= numel(x)
    error('evaluate_symbolic_at_point:DimensionMismatch', ...
          'Expected %d variables but received %d values.', ...
          numel(vars), numel(x));
end

value = double(subs(expr, vars, x));
end
