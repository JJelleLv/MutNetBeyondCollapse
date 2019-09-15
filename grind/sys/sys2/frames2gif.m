function frames2gif(filename,mov)
for i = 1:length(mov)
im = frame2im(mov(i));
[imind,cm] = rgb2ind(im,256);
if i == 1
  imwrite(imind,cm,filename,'gif','DelayTime',0.06, 'Loopcount',inf);
else
  imwrite(imind,cm,filename,'gif','DelayTime',0.06,'WriteMode','append');
end
end