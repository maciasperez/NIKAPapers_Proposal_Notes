#! /bin/csh -f
pdflatex desert.tex
pdflatex desert.tex
bibtex desert
pdflatex desert.tex
pdflatex desert.tex
rm -rf desert.aux
rm -rf desert.out
rm -rf desert.bbl
rm -rf desert.blg
rm -rf desert.log

