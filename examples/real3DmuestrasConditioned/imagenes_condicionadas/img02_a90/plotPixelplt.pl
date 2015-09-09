
$ex=$ARGV[0];

$nx=40;
$ny=60;
$nz=12;
#$nx=100;
#$ny=100;
#$nz=1;
$xmn=0.0;
$ymn=0.0;
$zmn=0.0;
$xsiz=10.0;
$ysiz=10.0;
$zsiz=10.0;


$cmin=-1.0;
$cmax=3.0;
$inc=0.02;

#open(PIXELPLT,">pixelplt_target_nscore.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#targetimage.nscore               -file with gridded data
#5                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#targetimage-nscore_$ex.ps            -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Target (nscore), slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);
#
#
#open(PIXELPLT,">pixelplt_initial_nscore.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#initialimage.nscore               -file with gridded data
#4                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#initialimage-nscore_$ex.ps            -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Initial (nscore), slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);
#
#
#
#open(PIXELPLT,">pixelplt_current_nscore.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#currentimage.nscore               -file with gridded data
#4                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#currentimage-nscore_$ex.ps            -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Current (nscore), slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);
#
#
#
#open(PIXELPLT,">pixelplt_target.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#targetimage.dat               -file with gridded data
#4                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#targetimage_$ex.ps            -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Target, slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);
#
#
#open(PIXELPLT,">pixelplt_initial.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#initialimage.dat               -file with gridded data
#4                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#initialimage_$ex.ps            -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Initial, slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);


for($slice=1;$slice<=12;$slice++){


open(PIXELPLT,">pixelplt_current_scaled_slice$slice.par");
$text=<<EOF; 
                  Parameters for PIXELPLT
                  ***********************

START OF PARAMETERS:
currentimage_scaled.dat               -file with gridded data
4                                -  column number for variable
-1.0e21  1.0e21                  -  data trimming limits
currentimage_scaled_$slice-slice_$ex.ps  -file with PostScript output
1                                -realization number
$nx  $xmn  $xsiz                 -nx,xmn,xsiz
$ny  $ymn  $ysiz                 -ny,ymn,ysiz
$nz  $zmn  $zsiz                 -nz,zmn,zsiz
1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
$slice                           -slice number
Current, slice=$slice    -Title
East                             -X label
North                            -Y label
0                                -0=arithmetic, 1=log scaling
1                                -0=gray scale, 1=color scale
0                                -0=continuous, 1=categorical
$cmin  $cmax  $inc               -continuous:  min, max, increm.
4                                -categorical: number of categories
1     3    Code_One              -category(), code(), name()
2     1    Code_Two 
3     6    Code_Three
4    10    Code_Four


Color Codes for Categorical Variable Plotting:
      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
      13=pink, 14=intermediate green, 15=gray
EOF
print PIXELPLT $text,"\n";
close(PIXELPLT);

open(PIXELPLT,">pixelplt_conditioned_scaled_slice$slice.par");
$text=<<EOF; 
                  Parameters for PIXELPLT
                  ***********************

START OF PARAMETERS:
conditionedimage_scaled.dat               -file with gridded data
4                                -  column number for variable
-1.0e21  1.0e21                  -  data trimming limits
conditionedimage_scaled_$slice-slice_$ex.ps  -file with PostScript output
1                                -realization number
$nx  $xmn  $xsiz                 -nx,xmn,xsiz
$ny  $ymn  $ysiz                 -ny,ymn,ysiz
$nz  $zmn  $zsiz                 -nz,zmn,zsiz
1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
$slice                           -slice number
Conditional, slice=$slice    -Title
East                             -X label
North                            -Y label
0                                -0=arithmetic, 1=log scaling
1                                -0=gray scale, 1=color scale
0                                -0=continuous, 1=categorical
$cmin  $cmax  $inc               -continuous:  min, max, increm.
4                                -categorical: number of categories
1     3    Code_One              -category(), code(), name()
2     1    Code_Two 
3     6    Code_Three
4    10    Code_Four


Color Codes for Categorical Variable Plotting:
      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
      13=pink, 14=intermediate green, 15=gray
EOF
print PIXELPLT $text,"\n";
close(PIXELPLT);



#$ret=`../../gslib90/pixelplt pixelplt_target.par`;
#$ret=`../../gslib90/pixelplt pixelplt_initial.par`;
#$ret=`../../gslib90/pixelplt pixelplt_current.par`;
$ret=`/home/operedo/Dropbox/amtc/apps/gslib90/pixelplt pixelplt_current_scaled_slice$slice.par`;
$ret=`/home/operedo/Dropbox/amtc/apps/gslib90/pixelplt pixelplt_conditioned_scaled_slice$slice.par`;

#$ret=`../../gslib90/pixelplt pixelplt_target_nscore.par`;
#$ret=`../../gslib90/pixelplt pixelplt_initial_nscore.par`;
#$ret=`../../gslib90/pixelplt pixelplt_current_nscore.par`;


}

#$slice=6;
#open(PIXELPLT,">pixelplt_current_scaled_slice$slice.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#currentimage_scaled.dat               -file with gridded data
#4                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#currentimage_scaled_$slice-slice_$ex.ps  -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Current, slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);
#
#open(PIXELPLT,">pixelplt_conditioned_scaled_slice$slice.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#conditionedimage_scaled.dat               -file with gridded data
#4                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#conditionedimage_scaled_$slice-slice_$ex.ps  -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Conditional, slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);
#
#
#
#
#$ret=`/home/operedo/Dropbox/amtc/apps/gslib90/pixelplt pixelplt_current_scaled_slice$slice.par`;
#$ret=`/home/operedo/Dropbox/amtc/apps/gslib90/pixelplt pixelplt_conditioned_scaled_slice$slice.par`;
#
#
#$slice=12;
#open(PIXELPLT,">pixelplt_current_scaled_slice$slice.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#currentimage_scaled.dat               -file with gridded data
#4                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#currentimage_scaled_$slice-slice_$ex.ps  -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Current, slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);
#
#open(PIXELPLT,">pixelplt_conditioned_scaled_slice$slice.par");
#$text=<<EOF; 
#                  Parameters for PIXELPLT
#                  ***********************
#
#START OF PARAMETERS:
#conditionedimage_scaled.dat               -file with gridded data
#4                                -  column number for variable
#-1.0e21  1.0e21                  -  data trimming limits
#conditionedimage_scaled_$slice-slice_$ex.ps  -file with PostScript output
#1                                -realization number
#$nx  $xmn  $xsiz                 -nx,xmn,xsiz
#$ny  $ymn  $ysiz                 -ny,ymn,ysiz
#$nz  $zmn  $zsiz                 -nz,zmn,zsiz
#1                                -slice orientation: 1=XY, 2=XZ, 3=YZ
#$slice                           -slice number
#Conditional, slice=$slice    -Title
#East                             -X label
#North                            -Y label
#0                                -0=arithmetic, 1=log scaling
#1                                -0=gray scale, 1=color scale
#0                                -0=continuous, 1=categorical
#$cmin  $cmax  $inc               -continuous:  min, max, increm.
#4                                -categorical: number of categories
#1     3    Code_One              -category(), code(), name()
#2     1    Code_Two 
#3     6    Code_Three
#4    10    Code_Four
#
#
#Color Codes for Categorical Variable Plotting:
#      1=red, 2=orange, 3=yellow, 4=light green, 5=green, 6=light blue,
#      7=dark blue, 8=violet, 9=white, 10=black, 11=purple, 12=brown,
#      13=pink, 14=intermediate green, 15=gray
#EOF
#print PIXELPLT $text,"\n";
#close(PIXELPLT);
#
#
#
#
#$ret=`/home/operedo/Dropbox/amtc/apps/gslib90/pixelplt pixelplt_current_scaled_slice$slice.par`;
#$ret=`/home/operedo/Dropbox/amtc/apps/gslib90/pixelplt pixelplt_conditioned_scaled_slice$slice.par`;
#
#
