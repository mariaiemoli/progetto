SUFFIXES=.pdf .eps .tex .dia .fig .svg

FIGFIGURES=${wildcard *.fig}
DIAFIGURES=${wildcard *.dia}
SVGFIGURES=${wildcard *.svg}
PNGFIGURES=${wildcard *.png}
ODGFIGURES=${wildcard *.odg}

FIGEPS=$(DIAFIGURES:.dia=.eps)
FIGEPS+=$(FIGFIGURES:.fig=.eps)
FIGEPS+=$(SVGFIGURES:.svg=.eps)
FIGEPS+=$(PNGFIGURES:.png=.eps)
FIGEPS+=$(ODGFIGURES:.odg=.eps)
FIGPDF=$(FIGFIGURES:.fig=.pdf)
FIGPDF+=$(DIAFIGURES:.dia=.pdf)
FIGPDF+=$(SVGFIGURES:.svg=.pdf)
FIGPDF+=$(PNGFIGURES:.png=.pdf)
FIGPDF+=$(ODGFIGURES:.odg=.pdf)
all: pdf

ps: $(FIGEPS)

pdf: $(FIGPDF)


clean:
	for f in *.fig; do rm -f `basename $$f '.fig'`.eps; done
	for f in *.fig; do rm -f `basename $$f '.fig'`.pdf; done
	for f in *.svg; do rm -f `basename $$f '.svg'`.pdf; done
	for f in *.svg; do rm -f `basename $$f '.svg'`.eps; done
	for f in *.dia; do rm -f `basename $$f '.dia'`.pdf; done
	for f in *.dia; do rm -f `basename $$f '.dia'`.eps; done
	for f in *.png; do rm -f `basename $$f '.png'`.eps; done
	for f in *.png; do rm -f `basename $$f '.png'`.pdf; done
	for f in *.odg; do rm -f `basename $$f '.png'`.eps; done
	for f in *.odg; do rm -f `basename $$f '.png'`.pdf; done
%.eps:%.fig
	fig2dev -L eps $< $@

%.pdf:%.fig
	fig2dev -L pdf $< $@


%.eps:%.dia
	dia -n --export=$@ $<


%.eps:%.png
	convert  $< $@

%.pdf:%.png
	convert $< $@

%.pdf:%.svg
	inkscape -D --file=$<  --export-pdf=$@

%.eps:%.svg
	inkscape -D --file=$<  --export-eps=$@

%.pdf:%.eps
	epstopdf --outfile=$@ $< 

%.pdf:%.odg
	/usr/lib/libreoffice/program/sdraw --nospash --nologo --terminate-after-init --convert-to pdf $< 