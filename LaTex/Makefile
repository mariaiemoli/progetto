.PHONY: main.pdf pdf clean distclean

main.pdf: main.tex bibliografia.bib 
		@pdflatex main
		@bibtex main
		@pdflatex main-frn
		@pdflatex  main
#		texi2pdf main

pdf: main.pdf

clean: 
	rm -f *.aux *.log *.out
        
distclean: 
	rm -f *.aux *.log *.out *.pdf *.gz