clear all
eps=[0.005 0.0075 0.01];
 snrdb = [10: 5 : 25];
for i = 1 : length(eps)
    [random_result(i,:) my_result(i,:)] = SCA (eps(i),50,0.5);
end
    