file = fopen('freq.txt','r');
fid = fscanf(file,'%f',[8,8]);
fclose(file);

exact = fid(1,:);
elem1 = fid(2,:);
elem2 = fid(3,:);
elem3 = fid(4,:);
elem4 = fid(5,:);
elem5 = fid(6,:);
elem10 = fid(7,:);
elem15 = fid(8,:);

i = 8;

e = exact(i);
mode1(1) = abs(e-elem1(i))/e;
mode1(2) = abs(e-elem2(i))/e;
mode1(3) = abs(e-elem3(i))/e;
mode1(4) = abs(e-elem4(i))/e;
mode1(5) = abs(e-elem5(i))/e;
mode1(6) = abs(e-elem10(i))/e;
mode1(7) = abs(e-elem15(i))/e;

plotArray = [1, mode1(1); 2, mode1(2); 3, mode1(3); 4, mode1(4); 5, mode1(5); 10, mode1(6);15, mode1(7)];
plot(plotArray(:,1), plotArray(:,2),'*'); grid on;

xlabel('Number of elements'); ylabel('Relative error on eigenfrequency');