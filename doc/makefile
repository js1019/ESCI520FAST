NAME = ESCI520FAST
pdf: $(NAME).tex
	pdflatex $(NAME); bibtex $(NAME); pdflatex $(NAME)
clean:
	rm -f *.aux *.dvi *.gz *.log *.blg *bbl *.aut *.pdf $(NAME).out *.backup
