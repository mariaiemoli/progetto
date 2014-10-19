#set terminal epslatex color solid
set terminal png
#set output "../Parts/ApplicativeExamples/Images/pressure.tex"

set output "./pressure.png"

data = "pressure.txt"

set size 1,1

set palette rgbformulae 33,13,10

#set bmargin at screen 0

set view equal xy
#set view 0,0,20

set view 67,326

#set border 0
#unset key
#unset tics
set cbtics 0.25

splot data with lines pal

set terminal x11
