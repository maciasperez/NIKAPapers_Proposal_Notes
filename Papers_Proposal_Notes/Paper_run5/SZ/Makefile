FICHIER_TEX= NIKA_tSZ_observation.tex
FICHIER_BIB= biblio.bib
ARCHIVE=$(FICHIER_TEX:%.tex=%.zip)


all : 
	pdflatex ${FICHIER_TEX}
	bibtex $(FICHIER_TEX:%.tex=%)
	pdflatex ${FICHIER_TEX}
################################################
#lance le make de chaque sous rep
figltex : 
	grep includegraphics *.tex > liste.txt
	awk  '{st=index($$0,"{");rem=substr($$0,st+1); st=index(rem,"}"); rem=substr(rem,1,st-1); printf("pdf2ps  %s.pdf %s.ps\nconvert %s.png %s.ps\nps2eps %s.ps\nrm %s.ps\n",rem,rem,rem,rem,rem,rem);}' liste.txt > conv_liste.sh
	bash conv_liste.sh

ltex : 
	latex ${FICHIER_TEX}
	bibtex $(FICHIER_TEX:%.tex=%)
	dvips $(FICHIER_TEX:%.tex=%.dvi)
	ps2pdf $(FICHIER_TEX:%.tex=%.ps)

clean: $(REP)
	@echo "suppression des fichiers de compilation"
	@# fichiers de compilation latex
	@rm -f *.log *.aux *.dvi *.toc *.lot *.lof *.out
	@# fichiers de bibtex
	@rm -f *.bbl *.blg
	@# fichier zip 
	@rm -f *.zip

zip: 
	#rm ${ARCHIVE}
	grep includegraphics *.tex > liste.txt
	awk  '{st=index($$0,"{");rem=substr($$0,st+1); st=index(rem,"}"); rem=substr(rem,1,st-1); printf("zip -r ${ARCHIVE} %s.p*\nzip -r ${ARCHIVE} %s.j*\n",rem,rem);}' liste.txt > zip_liste.sh
	zip -r ${ARCHIVE} Makefile
	zip -r ${ARCHIVE} *.cls
	zip -r ${ARCHIVE} *.tex
	bash zip_liste.sh

zipltex:
	#rm ${ARCHIVE}
	grep includegraphics *.tex > liste.txt
	awk  '{st=index($$0,"{");rem=substr($$0,st+1); st=index(rem,"}"); rem=substr(rem,1,st-1); printf("zip -r ${ARCHIVE} %s.eps\n",rem);}' liste.txt > zip_liste.sh
	zip -r ${ARCHIVE} ${FICHIER_BIB}
	zip -r ${ARCHIVE} Makefile
	zip -r ${ARCHIVE} ${FICHIER_TEX}
	bash zip_liste.sh
