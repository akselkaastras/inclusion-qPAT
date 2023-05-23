function results = burnthin(results,burnN,thinrate)
XR = results.XR;
results.XR_burnthin = XR(burnN:thinrate:end,:);