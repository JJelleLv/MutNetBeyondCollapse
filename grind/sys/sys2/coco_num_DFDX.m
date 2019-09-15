function J = coco_num_DFDX(F, x, p, varargin)
[m, n] = size(x);
fr = F(x(:, 1), p(:, 1), varargin{:});
l  = size(fr, 1);
J  = zeros(l, m, n);
for j = 1:n
    x0 = x(:, j);
    h  = 1.0e-8 * ( 1.0 + abs(x0) );
    hi = 0.5 ./ h;
    for i = 1:m
        xx       = x0;
        xx(i)    = x0(i) + h(i);
        fr       = F(xx, p(:, j), varargin{:});
        xx(i)    = x0(i) - h(i);
        fl       = F(xx, p(:, j), varargin{:});
        J(:, i, j) = hi(i) * (fr - fl);
    end
end
