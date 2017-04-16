name = manual
all:
	pdflatex  --shell-escape $(name).tex
	pdflatex  --shell-escape $(name).tex
	rm *.log *.aux
clean:
	rm -f *.out *.bbl *.blg *.log *.aux *.spl
	find . -name *~ -exec rm {} \;
