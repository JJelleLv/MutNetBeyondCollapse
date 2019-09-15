function cobweb=i_getcobweb(yy)
cobweb = zeros(size(yy, 1) * 2 - 1, 2);
cobweb(1, 1) = yy(1);
cobweb(1, 2) = 0;
for i = 2:size(yy, 1)
   cobweb((i - 1) * 2, 1) = yy(i - 1);
   cobweb((i - 1) * 2, 2) =yy(i);
   cobweb((i - 1) * 2 + 1, 1) = yy(i);
   cobweb((i - 1) * 2 + 1, 2) = yy(i);
end


