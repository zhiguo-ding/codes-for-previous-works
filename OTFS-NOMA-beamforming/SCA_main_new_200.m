clear all
eps=[0.005 0.0075 0.01];
 snrdb = [10: 5 : 25];
% for i = 1 : length(eps)
%     [random_result(i,:) my_result(i,:)] = SCA (eps(i),200,0.5,4);
% end

for i = 1 : length(eps)
    [random_result2(i,:) my_result2(i,:)] = SCA (eps(i),400,0.5,4);
end
for i = 1 : length(eps)
    [random_result1(i,:) my_result1(i,:)] = SCA (eps(i),400,0.5,8);
end

%     
% for i = 1 : length(eps)
%     [random_result1(i,:) my_result1(i,:)] = SCA (eps(i),100,0.5,4);
% end
%   for i = 1 : length(eps)
%     [random_result2(i,:) my_result2(i,:)] = SCA (eps(i),100,0.2,4);
% end
%       