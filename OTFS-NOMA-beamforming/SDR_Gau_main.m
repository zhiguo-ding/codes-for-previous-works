clear all
eps=[0.005 0.0075 0.01];
 snrdb = [10: 5 : 25]; 
 
 
[random_result1 SCO1 SDR1] = SRD_Gaus(8,200,0.5,10);
[random_result2 SCO2 SDR2] = SRD_Gaus(4,200,0.5,10); 

     
%[random_result1 SCO1 SDR1] = SRDm(8,10,0.2);

%[random_result2 SCO2 SDR2] = SRDm(4,10,0.2);

  [random_result3 SCO3 SDR3] = SRDm(8,200,0.5);

[random_result4 SCO4 SDR4] = SRDm(4,200,0.5);
