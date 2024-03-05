M=[[1:24].' randi(100,24,5)]       % a sample array
b=[1,4,7,8,10,12];                  % the indexing vector
N=M(b,:) 


% Example data
data = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];

% List of indices to select
indexList = [2, 5, 8];

% Select elements at the specified indices
selectedData = data(indexList);

% Display the result
disp(selectedData);




t = 0:0.01:5*pi;
y = sin(t);
z = 0:pi/12:2*pi;
disp(z)
clf()
axes()
hold on
for i = 1:numel(z)
    %disp(z(i)*ones(size(t)))

    disp(size(t))
    disp(size((z(i)*ones(size(t)))))
    disp(size(y))
    plot3(t,z(i)*ones(size(t)),y);
end
grid on
xlabel('t')
ylabel('o.5 pi')  % Replace with ylabel(['0.5',char(960)]) for pi symbol
zlabel('sin(t)')
view(12.6, 27.6)



value = 42;  % Replace with your actual value
valueSign = sign(value);

disp(['The sign of ', num2str(value), ' is ', num2str(valueSign)]);