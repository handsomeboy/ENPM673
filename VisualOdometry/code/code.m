clear all;
close all;
warning off;





Colour=dir('E:\Sem 1\ENPM 673\project 2\code\*.png');


%% Finding camera parameters 
    
    [fx, fy, cx, cy, G_camera_image, LUT] = ReadCameraModel('E:\Sem 1\ENPM 673\project 2\Oxford_dataset\Oxford_dataset\stereo\centre','E:\Sem 1\ENPM 673\project 2\Oxford_dataset\Oxford_dataset\model');

    K = [fx 0 cx; 0 fy cy; 0 0 1];
    cameraParams = cameraParameters('IntrinsicMatrix',K');
 %% pos and relative pose are defined 
pos = [0 0 0];
Rpos = [1 0 0;
        0 1 0;
        0 0 1];
post = [0 0 0];
Rpost = [1 0 0;
        0 1 0;
        0 0 1];
   


for i=37:39
    i
    I1 = imread(sprintf('%d.png',i));
    I2 = imread(sprintf('%d.png',i+1));

    
%% Undistorting images and then appling Gaussian filter to those images  
    [UI1] = UndistortImage(I1, LUT);
    [UI2] =  UndistortImage(I2,LUT);
    
 
    Ig1 = rgb2gray(UI1);
    Ig2 = rgb2gray(UI2);
    Ig1 = imgaussfilt(Ig1, .8);
    Ig2 = imgaussfilt(Ig2, .8);
%% Detecting corresponding points   
    points1 = detectSURFFeatures(Ig1);
    points2 = detectSURFFeatures(Ig2);

    [features1,validpoints1] = extractFeatures(Ig1,points1);
    [features2,validpoints2] = extractFeatures(Ig2,points2);

    indexPairs = matchFeatures(features1,features2);

    matchedPoints1 = validpoints1(indexPairs(:,1),:);
    matchedPoints2 = validpoints2(indexPairs(:,2),:);

   

    matches1 = matchedPoints1.Location;
    matches2 = matchedPoints2.Location;
   
    
    matches_1_x = matches1(:,1);
    matches_1_y = matches1(:,2);
    matches_2_x = matches2(:,1);
    matches_2_y = matches2(:,2);

   
    s = size(matches_1_x);
    
     
    X1 = [matches_1_x'; matches_1_y']';
    X2 = [matches_2_x'; matches_2_y']';
     

     [Tl,X1n] = normalize2(X1);
     [Tr,X2n] = normalize2(X2);
     
     
 %% Finding inlier points by applying RANSAC and then finding Fundamental matrix

num_samples = 1200;   
thres = .001;
bestinliers = 0;
for i = 1: num_samples
    
    testindices = randperm(size(X1,1), 8);
    testsample1 = X1(testindices, :);
    testsample2 = X2(testindices, :);
     [T1,testsample1n] = normalize2(testsample1);
     [T2,testsample2n] = normalize2(testsample2);
    

     testF = FindFundamentalMatrix(testsample1n, testsample2n); 
    
    
  
    testcountinliers = 0;
    testinlierindices = [];

    difftemp = zeros(8, 1) + 1; 
   for j = 1: 8
       eval = [testsample2n(j,:) 1] * testF * [testsample1n(j,:) 1]';
       if (abs(eval) < thres)
            testcountinliers = testcountinliers + 1;
            testinlierindices = [testinlierindices; testindices(1,j)];
            difftemp(j) = abs(eval); 
       end
   end    
   if (testcountinliers > bestinliers)

      bestinliers = testcountinliers;
       B_F_temp = T2'*testF*T1;
       B_F = B_F_temp/norm(B_F_temp);
       if B_F(end) < 0
       B_F = -B_F; 
    end
       diffbest = difftemp;
       best_inliers = testinlierindices;
       if bestinliers == 8
           break
       end
   end
end
    inliers1 = X1(best_inliers,:);
    inliers2 = X2(best_inliers,:);
    p1 = [inliers1(1,:) 1]';
    p2 = [inliers2(1,:) 1]';
 %% Finding Essential matrix    
    E = K' * B_F * K;
    
    
    [U, D, V] = svd(E);
    e = (D(1,1) + D(2,2)) / 2;
    D(1,1) = 1;
    D(2,2) = 1;
    D(3,3) = 0;
    E = U * D * V';

   %% Finding R and T from the Essential matrix 
    [U, ~, V] = svd(E);

    W = [0 -1 0;
         1 0 0; 
         0 0 1];
     
    R1 = U * W * V';
    if det(R1) < 0
        R1 = -R1;
    end
    R2 = U * W' * V';
    if det(R2) < 0
        R2 = -R2;
    end
    
    T1 = U(:,3)';
    T2 = -T1;
    
    Rs = cat(3, R1, R1, R2, R2);
    Ts = cat(1, T1, T2, T1, T2);
    
    %% Choose the right solution for rotation and translation with the help of triangulation
    numNeg = zeros(1, 4);
    P1 = cameraMatrix(cameraParams, eye(3), [0,0,0]);
 
    for k = 1:size(Ts, 1)
       P2 = cameraMatrix(cameraParams,Rs(:,:,k)', Ts(k, :));
    
       p3Dp_1 = zeros(size(inliers1, 1), 3, 'like', inliers1);
       P11 = P1';
       P22 = P2';
       
       M11 = P11(1:3, 1:3);
       M22 = P22(1:3, 1:3);
       
       c11 = -M11 \ P11(:,4);
       c22 = -M22 \ P22(:,4);

       for b = 1:size(inliers1,1)
          u11 = [inliers1(b,:), 1]';
          u22 = [inliers2(b,:), 1]';
          a11 = M11 \ u11;
          a22 = M22 \ u22;
          A = [a11, -a22];
          y = c22 - c11;
          
          alpha = (A' * A) \ A' * y;
          p = (c11 + alpha(1) * a11 + c22 + alpha(2) * a22) / 2;
          p3Dp_1(b, :) = p';
       end
      
       p3Dp_2 = bsxfun(@plus, p3Dp_1 * Rs(:,:,k)', Ts(k, :));
       numNeg(k) = sum((p3Dp_1(:,3) < 0) | (p3Dp_2(:,3) < 0));
    end
    
    [val, idx] = min(numNeg);
    R = Rs(:,:,idx)';
    T = Ts(idx, :);
    tNorm = norm(T);
    if tNorm ~= 0
        T = T ./ tNorm;
    end
    
    
    %% R & T
    R = R';
    T = -T * R;
    
    %% Finding the trajectory of the camera
    Rpos = R * Rpos;
    
    pos = pos + T * Rpos;
    
    %% Finding trajectory of the camera using MATLAB Image Processing Toolbox
    
    Fthe = estimateFundamentalMatrix(inliers1,inliers2,'Method','Norm8Point');
    
    [relativeOrientation,relativeLocation] = relativeCameraPose(Fthe,cameraParams,inliers1,inliers2);
    
    Rpost = relativeOrientation * Rpost;
    
    post = post + relativeLocation * Rpost;
   %% Plotting the found trajectories 
    figure(8)
    subplot(2,2,[3,4])
    title('Video From Camera')
    imshow(UI2)
    subplot(2,2,1)
    title('Path')
    
    plot(pos(1),pos(3),'+','MarkerEdgeColor','r')
    
    hold on;
    subplot(2,2,2)
    
    title('Path obtained by MATLAB image processing toolbox')
    plot(post(1),post(3),'+','MarkerEdgeColor','g')
   
    hold on;
    
    
    
    pause(0.005)

end 








    
    
    
    
     
     
     