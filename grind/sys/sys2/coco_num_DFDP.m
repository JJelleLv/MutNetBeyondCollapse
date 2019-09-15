function J = coco_num_DFDP(F, x, p, varargin)
x = x(:, :);
p = p(:, :);
m = size(p, 1);
n = size(x, 2);
fr = F(x(:, 1), p(:, 1), varargin{:});
l  = size(fr, 1);
J  = zeros(l, m, n);
for j = 1:n
    p0 = p(:, j);
    h  = 1.0e-8 * ( 1.0 + abs(p0) );
    hi = 0.5 ./ h;
    for i = 1:m
        p1    = p0;
        p1(i) = p0(i) + h(i);
        fr    = F(x(:, j), p1, varargin{:});
        p1(i) = p0(i) - h(i);
        fl    = F(x(:, j), p1, varargin{:});
        J(:, i, j) = hi(i) .* (fr - fl);
    end
end
