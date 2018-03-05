% AAE 568 Homework 2, Problem 5
% Reese Robinson
% February 26, 2018

clear all

%Load variables
load('HW2_Problem5.mat', '-mat')
y = [y1; y2]; %concatenated values, 500 x 1 vector
H = ones(1, 500);
H = H'; %500 x 1 matrix, full of 1's
% Create the R matrix (diagonal with variances in middle
diagStDev1 = ones(1, 320)*stdev1;
diagStDev1 = diagStDev1';
diagStDev1 = diag(diagStDev1);

diagStDev2 = ones(1, 180)*stdev2;
diagStDev2 = diagStDev2';
diagStDev2 = diag(diagStDev2);

R = blkdiag(diagStDev1, diagStDev2);

%First, solve for Least Squares Estimate of x
tic
xLS_NoWeight = (H'*H)^-1*H'*y;
toc
fprintf('The Unweighted Least Square estimate is: %f \n', xLS_NoWeight);

%Solve for Weighted Least Squares Estimate of x
tic
xLS_Weighted = (H'*R^-1*H)^-1*H'*R^-1*y;
toc
fprintf('The Weighted Least Square estimate is: %f \n', xLS_Weighted);

%Solve Recursive Least Squares Estimate
figure('Name', 'Unweighted Recursive LS');
k = [1, .1, .01];
disp('Computing Recursive LS')
for j = 1:3
    clear K R H S W y_RLS x_RLS_NW
 
    tic
    K = k(j)^2;
    R = diag(1);
    H = [1]';
    S = H*K*H' + R;
    W = K*H'*S^-1;
    y_RLS = [y(1)]';

    x_RLS_NW(1) = xhat_init + W*(y_RLS - H*xhat_init);
    K = K - W*S*W';
    trackK = [K];
  
    for i = 2:500
       H = [H; 1];
       R = blkdiag(R, [1]);
       S = H*K*H' + R;
       W = K*H'*S^-1;
       y_RLS = [y_RLS; y(i)];
       K = K - W*S*W';
  
       x_RLS_NW(i) = x_RLS_NW(i-1) + W*(y_RLS - H*x_RLS_NW(i-1));
    end
    toc
    hold on
    plot(x_RLS_NW)
end

xLS_NoWeight = xLS_NoWeight*ones(500);
plot(xLS_NoWeight)
xLS_Weighted = xLS_Weighted*ones(500);
plot(xLS_Weighted)
legend('K = 1^2', 'K = 0.1^2', 'K = 0.01^2', 'Batch LS', 'Batch Weighted LS')
ylabel({'$\hat{X}_{LS}$'}, 'Interpreter', 'latex')
xlabel('Iterations')
title('Recursive Least Square estimates')


%Solve Recursive Weighted Least Squares Estimate
k = [1, .1, .01];
disp('Computing Recursive Weighted LS')
figure('Name', 'Weighted Recursive LS');
for j = 1:3
    clear K R H S W y_RLS
    
    tic
    K = k(j)^2;
    R = diag(stdev1);
    H = [1]';
    S = H*K*H' + R;
    W = K*H'*S^-1;
    y_RLS = [y(1)]';

    x_RLS_W(1) = xhat_init + W*(y_RLS - H*xhat_init);
    K = K - W*S*W';
    trackK = [K];
  
    for i = 2:500
       H = [H; 1];
       if i < 321
           R = blkdiag(R, stdev1);
       else
           R = blkdiag(R, stdev2);
       end
       S = H*K*H' + R;
       W = K*H'*S^-1;
       y_RLS = [y_RLS; y(i)];
       K = K - W*S*W';
  
       x_RLS_W(i) = x_RLS_W(i-1) + W*(y_RLS - H*x_RLS_W(i-1));
    end
    toc
    hold on
    plot(x_RLS_W)
end
plot(xLS_NoWeight)
plot(xLS_Weighted)
legend('K = 1^2', 'K = 0.1^2', 'K = 0.01^2', 'Batch LS', 'Batch Weighted LS')
ylabel({'$\hat{X}_{LS}$'}, 'Interpreter', 'latex')
xlabel('Iterations')
title('Recursive Weighted Least Square estimates')