function out = fibonacci(n)

%f = zeros(size(1,n+1));
f = zeros(size(1,n));

f(1) = 1;
f(2) = 1;

%for i = 3 : n+1
for i = 3 : n
    f(i) = f(i-1) + f(i-2);
    %str = [num2str(f(i))];
    %disp(str)
end

out = f(end);

end