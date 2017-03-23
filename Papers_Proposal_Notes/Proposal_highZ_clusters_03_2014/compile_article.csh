#! /bin/tcsh -f
#
# NAME:
#     compile_article.csh
#
# PUPROSE:
#     launch LaTeX and BibTeX 
#     during this time, I can write and think instead of waiting 
#     the end of the compilation
#
# EXAMPLE:
#      chod 755 compile_article.csh
#     ./compile_these.csh
#
# MODIFICATION HISTORY:
#     26-Jun-2000 Herve Dole Written, from compile_these
#
set article = 30m_planck_herschel_spitzer_201309

latex $article

latex $article

#bibtex $article

#bibtex $article

latex $article

latex $article




# Generate PDF of Summary
#------------------------
dvips -P pdf -t a4  $article
#dvips -P pdf  $article

ps2pdf13 $article.ps $article.pdf

rm -f $article.ps










