%% The following commented section reads the image from the input, converts it to RGB images and then save those RGB images
a=dir('E:\Sem 1\ENPM 673\project 2\Oxford_dataset\Oxford_dataset\stereo\centre /*.png');
sd=size(a,1);

for i = 1:sd
   filename = strcat('E:\Sem 1\ENPM 673\project 2\Oxford_dataset\Oxford_dataset\stereo\centre\',a(i).name);

  BW = imread(filename);
  RGB = demosaic (BW, 'gbrg');
  imwrite(RGB, sprintf('%d.png', RGB(i)));

end