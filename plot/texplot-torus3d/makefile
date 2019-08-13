all: torus_plot.eps 

torus_plot.eps : torus_plot.tex
	latex torus_plot.tex
	dvips -o torus_plot.eps torus_plot.dvi

clean:
	rm -f *.aux *.log *.dvi *.blg *.bbl *.pdf *.out
