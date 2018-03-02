print """set terminal png size 1280,960

set tmargin 5
set lmargin 10
set ylabel "Phase (degrees)" font ",30" offset 0, +1
set yrange [-180:180]
set bmargin 4
set xlabel "Frequency (GHz)" font ",30" offset 0, -1
set title "Phase slope after instrumental corrections" font "Times,40" offset 0,0
set xtic font ", 20"
set ytic font ", 20"
set key font ", 20" spacing 1.5
"""

antennas = ["Mc", "On", "Sr", "Sv", "Tr", "Ur", "Wb", "Zc"]
print "plot", ", ".join(['"scan48-all-before-s.txt" index {} using 1:2 title "Ef-{}"'.format(i, s)
                 for (i, s) in enumerate(antennas)])
