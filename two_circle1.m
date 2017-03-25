


rad1 = rand(1,1000)*(2*pi);
rr1 = rand(1,1000)*0.15+0.25;

x1 = zeros(1,1000);
y1 = zeros(1,1000);
for n = 1:1000
    x1(n) = rr1(n)*cos(rad1(n))+0.0;
    y1(n) = rr1(n)*sin(rad1(n))+0.0;
end
figure;
plot(x1,y1,'+');
axis([-1 1 -1 1])
hold on

%rad2 = (1+rand(1,1000))*(pi);
rad2 = rand(1,1000)*(2*pi);
rr2 = rand(1,1000)*0.2+0.0;

x2 = zeros(1,1000);
y2 = zeros(1,1000);
for n = 1:1000
    x2(n) = rr2(n)*cos(rad2(n))+0.0;
    y2(n) = rr2(n)*sin(rad2(n))+0.0;
end
plot(x2,y2,'x');
axis([-1 1 -1 1])

fid1=fopen('circle1.txt','w+');
n = 0;
for i = 1:1000
    %if(r1(i,1)<1 && r1(i,1)>0 && r1(i,2)<1 && r1(i,2)>0)
        n = n+1;
        fprintf('%g\n',n);
        %fprintf(fid1,'%g %g %g %g %g %g %g\n',cos(x1(i)*pi/2),sin(x1(i)*pi/2),cos(y1(i)*pi/2),sin(y1(i)*pi/2),1,x1(i),y1(i));
        fprintf(fid1,'%g %g %g\n',x1(i),y1(i),1);
        %fprintf('%g %g %g %g %g %g %g\n',cos(x1(i)*pi/2),sin(x1(i)*pi/2),cos(y1(i)*pi/2),sin(y1(i)*pi/2),1,x1(i),y1(i));
        %hold on
        %plot(x1(i),y   1(i),'x');
    %end
    %fprintf(fid1,'%g %g %g\n',cos(r1(i,2)*pi/2),sin(r1(i,2)*pi/2),1);
    %fprintf('%g\r\n',r1(:,2));
end
fclose(fid1);

fid2=fopen('circle2.txt','w+');
n = 0;
for i = 1:1000
    %if(r2(i,1)<1 && r2(i,1)>0 && r2(i,2)<1 && r2(i,2)>0)
        n = n+1;
        %fprintf('%g\n',n);
        %fprintf(fid2,'%g %g %g %g %g %g %g\n',cos(x2(i)*pi/2),sin(x2(i)*pi/2),cos(y2(i)*pi/2),sin(y2(i)*pi/2),2,x2(i),y2(i));
        fprintf(fid2,'%g %g %g\n',x2(i),y2(i),2);
        %hold on
        %plot(x2(i),y2(i),'x');
    %end
    %fprintf('%g\r\n',r1(:,2));
end
fclose(fid2);





