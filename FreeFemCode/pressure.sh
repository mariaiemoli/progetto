#!/bin/bash

gnuplot pressure.gnuplot

ps2pdf -dEPSCrop ../Parts/ApplicativeExamples/Images/pressure.eps \
../Parts/ApplicativeExamples/Images/pressure.pdf

rm ../Parts/ApplicativeExamples/Images/pressure.eps

sed -i 's,\.\.\/Parts\/ApplicativeExamples\/Images\/,,' ../Parts/ApplicativeExamples/Images/pressure.tex
