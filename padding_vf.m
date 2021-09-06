function vf = padding_vf(vfdata)
% pad the visual field data array to make it suitable for plotting in the
% same format as Elze et al. 2014. Returns a 8x9 matrix

vf = [NaN NaN NaN vfdata(1:4) NaN NaN;
      NaN NaN vfdata(5:10) NaN;
      NaN vfdata(11:18);
      vfdata(19:25) NaN vfdata(26);
      vfdata(27:33) NaN vfdata(34);
      NaN vfdata(35:42);
      NaN NaN vfdata(43:48) NaN;
      NaN NaN NaN vfdata(49:52) NaN NaN];
  
  return