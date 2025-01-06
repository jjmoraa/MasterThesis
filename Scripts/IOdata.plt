# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/ VERSION 20140514+ 
# frame analysis via Matlab interface 
# Sat Dec 14 15:53:29 2024
# G N U P L O T   S C R I P T   F I L E 
set autoscale
unset border
set pointsize 1.0
set xtics; set ytics; set ztics; 
unset zeroaxis
unset key
unset label
set size ratio -1    # 1:1 2D axis scaling 
# set view equal xyz # 1:1 3D axis scaling 
# NODE NUMBER LABELS
set label ' 1' at  0.0000e+000,  0.0000e+000,  0.0000e+000
set label ' 2' at  7.5642e+000,  0.0000e+000,  0.0000e+000
set label ' 3' at  1.5128e+001,  0.0000e+000,  0.0000e+000
set label ' 4' at  2.2693e+001,  0.0000e+000,  0.0000e+000
set label ' 5' at  3.0257e+001,  0.0000e+000,  0.0000e+000
set label ' 6' at  3.7821e+001,  0.0000e+000,  0.0000e+000
set label ' 7' at  4.5385e+001,  0.0000e+000,  0.0000e+000
set label ' 8' at  5.2949e+001,  0.0000e+000,  0.0000e+000
set label ' 9' at  6.0513e+001,  0.0000e+000,  0.0000e+000
set label ' 10' at  6.8078e+001,  0.0000e+000,  0.0000e+000
set label ' 11' at  7.5642e+001,  0.0000e+000,  0.0000e+000
set label ' 12' at  8.3206e+001,  0.0000e+000,  0.0000e+000
set label ' 13' at  9.0770e+001,  0.0000e+000,  0.0000e+000
set label ' 14' at  9.8334e+001,  0.0000e+000,  0.0000e+000
set label ' 15' at  1.0590e+002,  0.0000e+000,  0.0000e+000
set label ' 16' at  1.1346e+002,  0.0000e+000,  0.0000e+000
set label ' 17' at  1.2103e+002,  0.0000e+000,  0.0000e+000
set label ' 18' at  1.2859e+002,  0.0000e+000,  0.0000e+000
set label ' 19' at  1.3616e+002,  0.0000e+000,  0.0000e+000
set label ' 20' at  1.4372e+002,  0.0000e+000,  0.0000e+000
set label ' 21' at  1.5128e+002,  0.0000e+000,  0.0000e+000
# ELEMENT NUMBER LABELS
set label ' 1' at  3.7821e+000,  0.0000e+000,  0.0000e+000
set label ' 2' at  1.1346e+001,  0.0000e+000,  0.0000e+000
set label ' 3' at  1.8910e+001,  0.0000e+000,  0.0000e+000
set label ' 4' at  2.6475e+001,  0.0000e+000,  0.0000e+000
set label ' 5' at  3.4039e+001,  0.0000e+000,  0.0000e+000
set label ' 6' at  4.1603e+001,  0.0000e+000,  0.0000e+000
set label ' 7' at  4.9167e+001,  0.0000e+000,  0.0000e+000
set label ' 8' at  5.6731e+001,  0.0000e+000,  0.0000e+000
set label ' 9' at  6.4295e+001,  0.0000e+000,  0.0000e+000
set label ' 10' at  7.1860e+001,  0.0000e+000,  0.0000e+000
set label ' 11' at  7.9424e+001,  0.0000e+000,  0.0000e+000
set label ' 12' at  8.6988e+001,  0.0000e+000,  0.0000e+000
set label ' 13' at  9.4552e+001,  0.0000e+000,  0.0000e+000
set label ' 14' at  1.0212e+002,  0.0000e+000,  0.0000e+000
set label ' 15' at  1.0968e+002,  0.0000e+000,  0.0000e+000
set label ' 16' at  1.1724e+002,  0.0000e+000,  0.0000e+000
set label ' 17' at  1.2481e+002,  0.0000e+000,  0.0000e+000
set label ' 18' at  1.3237e+002,  0.0000e+000,  0.0000e+000
set label ' 19' at  1.3994e+002,  0.0000e+000,  0.0000e+000
set label ' 20' at  1.4750e+002,  0.0000e+000,  0.0000e+000
# set parametric
# set view 60, 70,  1.00 
# set view equal xyz # 1:1 3D axis scaling 
# unset key
# set xlabel 'x'
# set ylabel 'y'
# set zlabel 'z'
set title "frame analysis via Matlab interface \nanalysis file: IOdata.FMM   deflection exaggeration: 10.0   load case 1 of 1 "
unset clip; 
set clip one; set clip two
set xyplane 0 
  plot 'C:\Users\josej\Downloads\Frame3DD_20140514+\Frame3DD\temp\IOdata-msh' u 2:3 t 'undeformed mesh' w lp lw 1 lt 5 pt 6, 'C:\Users\josej\Downloads\Frame3DD_20140514+\Frame3DD\temp\IOdata-mshf.001' u 1:2 t 'load case 1 of 1' w l lw 1 lt 3
# splot 'C:\Users\josej\Downloads\Frame3DD_20140514+\Frame3DD\temp\IOdata-msh' u 2:3:4 t 'load case 1 of 1' w lp  lw 1 lt 5 pt 6, 'C:\Users\josej\Downloads\Frame3DD_20140514+\Frame3DD\temp\IOdata-mshf.001' u 1:2:3 t 'load case 1 of 1' w l lw 1 lt 3
