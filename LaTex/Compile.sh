#!/bin/bash

rm -r compile/*
cp -r *.tex compile
cp -r *.bib compile
#cp -r classi/* compile
cp -r img/ compile

cd compile

latex Main
bibtex Main
latex Main
latex Main
dvipdf Main

cp Main.pdf ../
cd ..

open Main.pdf
