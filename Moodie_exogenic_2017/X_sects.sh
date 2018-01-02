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
outfile=xsects.ps

# create path coordinates that the cross sections will be drawn along
project -C-3.75/38 -E-1.5/34 -G0.5 -Q > xsect.xy

# add MOHO, LAB, geoid, and topography sequentially to a text file
gmt grdtrack xsect.xy -GGAMOHO.nc > xsectMOHO.xy
gmt grdtrack xsectMOHO.xy -GGALABinv.nc > xsectMOHOandLAB.xy
gmt grdtrack xsectMOHOandLAB.xy -GGAegm08.nc > xsectMOHOandLABandGEOID.xy
gmt grdtrack xsectMOHOandLABandGEOID.xy -GGAGEBCO.nc > xsectMOHOandLABandGEOIDandTOPO.xy

# extract and plot geoid
awk -F "\t" '{print $3 "\t" $6}' xsectMOHOandLABandGEOIDandTOPO.xy > xsectGEOIDplot.xy
gmt psxy xsectGEOIDplot.xy -JX8i/1.5i -R0/500/-10/60 -W10.0pt,red -Bg50a100:"":/g10a20:"Depth (km)":WSne:.: -K -Xa1 -Ya6.3 > $xsectout
# extract and plot topo
awk -F "\t" '{print $3 "\t" $7}' xsectMOHOandLABandGEOIDandTOPO.xy > xsectTOPOplot.xy
gmt psxy xsectTOPOplot.xy -JX8i/1.5i -R0/500/-10/60 -W10.0pt,brown -K -O -Xa1 -Ya6.3 >> $xsectout
# extract and plot MOHO
awk -F "\t" '{print $3 "\t" $4}' xsectMOHOandLABandGEOIDandTOPO.xy > xsectMOHOplot.xy
gmt psxy xsectMOHOplot.xy -JX8i/4.8i -R0/500/-175/10 -W10.0pt,34/139/34 -Bg50a100:"Distance (km)":/g25a25:"Depth (km)":WSne:.: -K -O -Xa1 -Ya1 >> $xsectout
# draw topo again (sometimes needed for clarity)
gmt psxy xsectTOPOplot.xy -JX8i/4.8i -R0/500/-175/10 -W10.0pt,brown  -K -O -Xa1 -Ya1 >> $xsectout
# extract and plot LAB
awk -F "\t" '{print $3 "\t" $5}' xsectMOHOandLABandGEOIDandTOPO.xy > xsectLABplot.xy
gmt psxy xsectLABplot.xy -JX8i/4.8i -R0/500/-175/10 -W10.0pt,blue -O -Xa1 -Ya1 >> $xsectout
