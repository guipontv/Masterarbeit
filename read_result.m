
data = fopen('result.txt','r');
s = fgets(data);
B = [str2double(s)];

while ischar(s)
    s = fgets(data);
    B = cat(1,B,[str2double(s)]);
end
B = B(1:end-1,:);


figure
subplot(2,2,1)
plot(imag(B))
xlim([0 128])
subplot(2,2,2)
plot(imag(B))
xlim([500 628])
subplot(2,2,3)
plot(imag(B))
xlim([1000 1128])
subplot(2,2,4)
plot(imag(B))
xlim([2000 2128])