set term pngcairo size 800,800
#set term pdfcairo enhanced size 4in, 4in
set object 1 rect behind from screen 0,0 to screen 1,1 fc rgb "#000000" fillstyle solid 1.0
set border lc rgb "white"
set key tc rgb "white"
unset key
set xlabel "x" tc rgb "white"
set ylabel "y" tc rgb "white" offset 1, 0
set palette rgbformulae 22, 13, -31
set style fill border 6
set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set xtics 0.1
set ytics 0.1
set isosamples 1000
do for [i=0:10]{
#outfile = sprintf("barrett_mfs2%06d.png", i)
outfile = sprintf("barrett_mfs2.png")	
set output outfile
plot for[j=0:i] sprintf("barrett_mfs%06d.dat", 500 * j) using 1:2 w l lw 1.5 lc "skyblue"
#plot sprintf("yoko_kuro%06d.dat", 5 * i) using 1:2 lw 2 lc "skyblue" w l
#plot sprintf("yoko_kuro%06d.dat", i) using 1:2 with filledcurves lc "white" linewidth 3
}