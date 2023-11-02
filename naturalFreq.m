data2 = [-2.67	-2.73	1.43	1.55	1.71	1.9	2.05	2.14	-1.92	-1.99	-2.15	-2.32	-2.5	-2.63	-2.69	1.42	1.53	1.67	1.84	2	2.09	-1.92	-1.98	-2.11	-2.28	-2.46	-2.59	-2.65	1.41	1.49	1.64	1.8];
clf;
hold on 
plot(data2)
datamin = islocalmin(data2);%generates array of bools with ones at mins
time = zeros(size(data1));
for i = 1:size(data2, 2)
    time(1,i) = i * .1;
end
datamin = time.*datamin;
datamin = nonzeros(datamin);
period = (datamin(end)-datamin(1))/(size(datamin,1)-1);
freq = 1/period*2*pi;
g = 9.8;
l = g/((freq)^2);

