#!/bin/bash
gmt set PS_MEDIA letter
gmt set PROJ_LENGTH_UNIT inch
gmt set PS_PAGE_ORIENTATION landscape
gmt set FONT_LABEL 09p
gmt set MAP_FRAME_WIDTH 4p
gmt set FONT_ANNOT_PRIMARY 12p
gmt set FONT_LABEL 12p
gmt set MAP_TICK_LENGTH_PRIMARY 0.13i

# variables
infile=./GA_GEBCO.nc
outfile=GA_filt.ps
topocpt=../../gmt/mby.cpt
filtcpt=nogreen.cpt
div_color=black
synth_div_col=35/75/225
line_wght=6.0

# draw upper left, raw topo, unflitered
gmt grdgradient $infile -GGA-grad.nc -A345 -Ne0.6 -V
gmt grdhisteq GA-grad.nc -GGA-hist.nc -N -V
gmt grdmath GA-hist.nc 5.8 DIV = GA-int.nc
gmt grdimage $infile -R-7/0.5/33.5/38.5 -Jm0.5 -B0.5WesNa1 -Xa1 -Ya4.7 -IGA-int.nc -C$topocpt -V -K > $outfile
gmt psxy ./northdivide.xy -Jm0.5 -R -W$line_wght,$div_color -Xa1 -Ya4.7 -V -K -O >> $outfile
gmt psxy ./southdivide.xy -Jm0.5 -R -W$line_wght,$div_color -Xa1 -Ya4.7 -V -K -O >> $outfile
gmt set MAP_TICK_LENGTH_PRIMARY 0
gmt psscale -D4.05i/5.05i/0.6i/0.40i -S -L -C$topocpt -K -O >> $outfile
gmt set MAP_TICK_LENGTH 0.0787402i

# draw upper right, filter 1
gmt grdfft $infile -GGAfftur.nc -F-/-/50000/40000 -fg -V
gmt grdgradient GAfftur.nc -GGAfftur-grad.nc -A345 -Ne0.6 -V
gmt grdhisteq GAfftur-grad.nc -GGAfftur-hist.nc -N -V
gmt grdmath GAfftur-hist.nc 5.8 DIV = GAfftur-int.nc
gmt grdimage GAfftur.nc -R-7/0.5/33.5/38.5 -Jm0.5 -B0.5wEsNa1 -Xa5.5 -Ya4.7 -IGAfftur-int.nc -C$filtcpt -V -K -O >> $outfile
gmt pscoast -Jm0.5 -R-7/0.5/33.5/38.5 -Xa5.5 -Ya4.7 -Di -W2,black -A0/0/1 -K -O >> $outfile
gmt psxy ./northdivide.xy -Jm0.5 -R -W$line_wght,$div_color -Xa5.5 -Ya4.7 -V -K -O >> $outfile
gmt psxy ./southdivide.xy -Jm0.5 -R -W$line_wght,$div_color -Xa5.5 -Ya4.7 -V -K -O >> $outfile
gmt psxy north50k.xy -Jm0.5 -R -O -K -W$line_wght,$synth_div_col,-- -Xa5.5 -Ya4.7 -V >> $outfile
gmt psxy south50k.xy -Jm0.5 -R -O -K -W$line_wght,$synth_div_col,-- -Xa5.5 -Ya4.7 -V >> $outfile
gmt set MAP_TICK_LENGTH 0 
gmt psscale -D8.55i/5.05i/0.6i/0.40i -S -L -C$filtcpt -K -O >> $outfile
gmt set MAP_TICK_LENGTH 0.0787402i

# draw the scalebar
psbasemap -J -R -Lf3/34.6/38/500k+lKilometers+jb -O -P >> $outfile
