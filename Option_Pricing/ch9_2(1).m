clear;clc;
rng(0);
raw_data = textread('data1.txt');
randIndex = randperm(size(raw_data,1));
data_new=raw_data(randIndex,:);
train_data = data_new(1:size(raw_data,1)/2,:);
test_data = data_new(size(raw_data,1)/2+1:size(raw_data,1),:);
%%%%% task1 %%%%%
% 计算先验
p1 = sum(train_data(:,5)==0);
p2 = sum(train_data(:,5)==1);
p3 = sum(train_data(:,5)==2);
% 均值向量
class1 = train_data(train_data(:,5)==0,:);
class2 = train_data(train_data(:,5)==1,:);
class3 = train_data(train_data(:,5)==2,:);
omega1 = mean(class1(:,1:4));
omega2 = mean(class2(:,1:4));
omega3 = mean(class3(:,1:4));
% 协方差矩阵
sigma = cov(train_data(:,1:4));
% 测试
classify_res = [];
for i=1:size(test_data,1)
    res1 = -1/2*(test_data(i,1:4)-omega1)*inv(sigma)*(test_data(i,1:4)-omega1)'+log(p1);
    res2 = -1/2*(test_data(i,1:4)-omega2)*inv(sigma)*(test_data(i,1:4)-omega2)'+log(p2);
    res3 = -1/2*(test_data(i,1:4)-omega3)*inv(sigma)*(test_data(i,1:4)-omega3)'+log(p3);
    res = [res1 res2 res3];
    loc = find(res == min(res));
    classify_res = [classify_res loc];
end
correct_num = 0;
for i=1:size(test_data,1)
    if(classify_res(i) == test_data(i,5))
        correct_num = correct_num + 1;
    end
end
error_rate = (size(test_data,1)-correct_num)/size(test_data,1);

%%%%% task2 %%%%%
% % 计算先验
% p1 = sum(train_data(:,5)==0);
% p2 = sum(train_data(:,5)==1);
% p3 = sum(train_data(:,5)==2);
% % 均值向量
% class1 = train_data(train_data(:,5)==0,:);
% class2 = train_data(train_data(:,5)==1,:);
% class3 = train_data(train_data(:,5)==2,:);
% omega1 = mean(class1(:,1:4));
% omega2 = mean(class2(:,1:4));
% omega3 = mean(class3(:,1:4));
% % 协方差矩阵
% sigma1 = cov(class1(:,1:4));
% sigma2 = cov(class2(:,1:4));
% sigma3 = cov(class3(:,1:4));
% o1 = 1/2*omega1*inv(sigma1)*omega1' - 1/2*log(det(sigma1)) + log(p1);
% o2 = 1/2*omega2*inv(sigma2)*omega2' - 1/2*log(det(sigma2)) + log(p2);
% o3 = 1/2*omega3*inv(sigma3)*omega3' - 1/2*log(det(sigma3)) + log(p3);
% % 测试
% classify_res = [];
% for i=1:size(test_data,1)
%     res1 = test_data(i,1:4)*(-1/2*inv(sigma1))*test_data(i,1:4)'+(inv(sigma1)*omega1')'*test_data(i,1:4)'+o1;
%     res2 = test_data(i,1:4)*(-1/2*inv(sigma2))*test_data(i,1:4)'+(inv(sigma2)*omega2')'*test_data(i,1:4)'+o2;
%     res3 = test_data(i,1:4)*(-1/2*inv(sigma3))*test_data(i,1:4)'+(inv(sigma3)*omega3')'*test_data(i,1:4)'+o3;
%     res = [res1 res2 res3];
%     loc = find(res == min(res));
%     classify_res = [classify_res loc];
% end
% correct_num = 0;
% for i=1:size(test_data,1)
%     if(classify_res(i) == test_data(i,5))
%         correct_num = correct_num + 1;
%     end
% end
% error_rate = (size(test_data,1)-correct_num)/size(test_data,1);
