function X = vm_1x2_to_1x10(x)

X = zeros(1,9);
X(1:2) = -x;
X(3:5) = 1/2*[x.^2 sqrt(2)*x(1)*x(2)];
X(6:9) = -1/6*[x.^3 sqrt(3)*x(1).^2*x(2) sqrt(3)*x(1).*x(2).^2];
X=[1 X];
end
