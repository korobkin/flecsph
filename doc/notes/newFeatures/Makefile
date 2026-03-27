#TeX Makefile
FILE=gw

all: $(FILE).pdf

view: $(FILE).pdf display

.PHONY: clean

clean:
	rm -rf *~ *.idx *.lot *.toc *.lof *.aux *.blg *.out *.bbl *.log *.pdf

display:
	open $(FILE).pdf

$(FILE).pdf: $(FILE).tex
	pdflatex -shell-escape $(FILE)
	bibtex  $(FILE)
	bibtex  $(FILE)
	pdflatex -shell-escape $(FILE)
	pdflatex -shell-escape $(FILE)
