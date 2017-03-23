FICHIERS_ODP = $(wildcard *.odp)
FICHIERS_ODG = $(wildcard *.odg)
FICHIERS_ODS = $(wildcard *.ods)
FICHIERS_DOT = $(wildcard *.dot)
FICHIERS_SVG = $(wildcard *.svg)
FICHIERS_FIG = $(wildcard *.fig)

#FICHIERS_PDF=$(patsubst *.dot,%pdf, $(FICHIERS_DOT))
# 	fig2mpdf $<

FICHIERS_PDF=$(FICHIERS_DOT:%.dot=%.pdf) $(FICHIERS_FIG:%.fig=%.pdf) $(FICHIERS_SVG:%.svg=%.pdf) $(FICHIERS_ODG:%.odg=%.pdf) $(FICHIERS_ODS:%.ods=%.pdf) $(FICHIERS_ODP:%.odp=%.pdf)



#	@echo "creation figures"
#	@echo $(FICHIERS_PDF)
# $@ -> cible
# $< -> nom premiere dependance
# $^ -> liste des dependances
# $? -> dependance plus recente que la cible
# $* -> nom du fichier sans suffixe
#
all: $(FICHIERS_PDF)

clean :
	@echo "suppression des fichiers cibles"
	@rm -f *.pdf

cleaner: clean

zip:

%.pdf : %.dot
	@echo "generation pdf a partir de *.dot"
	dot -Tpdf $< -o $@

%.pdf : %.svg
	@echo "generation pdf a partir de *.svg"
	inkscape -A ~/PDF/$@ $<
	pdfcrop ~/PDF/$@
	mv ~/PDF/$*-crop.pdf  ./$*.pdf
	rm ~/PDF/$*.pdf

%.pdf : %.odg 
	@echo "generation pdf a partir de *.odg"
	libreoffice -pt "PDF" $<
	@echo "odg converted"
	@if test -f  $@ ; then echo "do crop"; else sleep 3; fi
	@sleep 1
	pdfcrop ~/PDF/$@
	mv ~/PDF/$*-crop.pdf  ./$*.pdf
	rm ~/PDF/$*.pdf

%.pdf : %.ods
	@echo "generation pdf a partir de *.ods"
	libreoffice -pt "PDF" $<
	@echo "ods converted"
	@if test -f  $@ ; then echo "do crop"; else sleep 3; fi
	@sleep 1
	pdfcrop ~/PDF/$@
	mv ~/PDF/$*-crop.pdf  ./$*.pdf
	rm ~/PDF/$*.pdf

%.pdf : %.odp
	@echo "generation pdf a partir de *.odp"
	libreoffice -pt "PDF" $<
	@echo "odp converted"
	@if test -f  $@ ; then echo "do crop"; else sleep 3; fi
	@sleep 1
	pdfcrop ~/PDF/$@
	mv ~/PDF/$*-crop.pdf  ./$*.pdf
	rm ~/PDF/$*.pdf

%.pdf : %.fig
	@echo "generation pdf a partir de *.fig"
	fig2dev -Lpdf $< $@
