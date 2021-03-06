global daytime  tempav  precav  fdParobj  place 

station = 1;
lambda  = 0.01; 
[H_f1, H_f2, Hc_OK, ...
      Hc_Next, Hc_Last, Hc_This, Hc_Enter, Hc_Quit,  ...
      Hc_Weather, Hc_Velocity, Hc_Accel, Hc_HarmAcc, Hc_Residual, ...
      Hc_Temp, Hc_Prec, Hc_CaseNo, Hc_Lambda] = ...
      dailysetup(station, lambda);
pause;
delete(H_f1);
delete(H_f2);
