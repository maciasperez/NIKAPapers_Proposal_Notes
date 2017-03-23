%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      sf2a-conf  class 2016 for LaTeX2e              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%-----------------------------------------------------------------------------
%% --> THE CLASS OPTION MATERIAL
%%-----------------------------------------------------------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{sf2a-conf2016}[]
\LoadClass[twoside]{article}
\PassOptionsToClass{twoside}{sf2a-conf2016}
\relax
%
%%-----------------------------------------------------------------------------
%% --> THE TOP MATTER MATERIAL
%%-----------------------------------------------------------------------------
%
%===  maketitle
\renewcommand\maketitle{\par
  \begingroup
    \let\@makefnmark\relax \let\@thefnmark\relax
    \@spcont
      \vbox{\hbox{\fontsize{10}{9}\selectfont\itshape\@NMJRNL@E\hfill}%
            \hbox{\fontsize{8}{6}\selectfont\itshape\@NMJRNL@F\hfill}}
    \@spcc
    \@maketitle
    \thispagestyle{plain}
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\title\relax
  \global\let\author\relax
  \global\let\and\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@title\@empty
  \normalfont\normalsize
}
\def\@maketitle{%
  \@settitle\@thanks\@setauthors
  }
%-A- Headings 
\def\@NMJRNL@F{L. Cambr\'esy, M. Deleuil, E. P\'econtal, C. Reyl\'e, J. Richard, L. Tresse and I. Vauglin (eds)}%
\def\@NMJRNL@E{SF2A 2016}%
\def\@NMJRNL@J{}%
\def\TitreGlobal#1{\global\def\@MSSG@TtrGlbl{#1}}
\def\FirstPage#1{\setcounter{page}{#1}}
%
%-B- title
\def\title#1{\global\def\@title{\uppercase{#1}}
  \global\def\leftmark{\def\\{\ }#1}\global\def\rightmark{\@MSSG@TtrGlbl}}
\newbox\@b@rnngttl
\def\runningtitle#1{\global\setbox\@b@rnngttl=\hbox{#1}
                    \ifdim\wd\@b@rnngttl>0.9\hsize\relax\else\global\def\@MSSG@rnnghd{#1}\fi}
\def\@settitle{\def\thanks##1{\unskip{}}\begin{center}\fontsize{12}{14}\selectfont\sffamily\bfseries
              \@title\end{center}\@pc}
%
%-C- authors and address
\def\@setauthors{\begingroup\trivlist
  \centering
  \item\relax\fontsize{11}{13}\selectfont\rmfamily\@@th@rs\ignorespaces
  \footnote{\box\@b@ddrss@}\endtrivlist\endgroup}
\def\email#1{{e-mail: \tt#1}}
\newcount\@c@thr@\@c@thr@=0
\def\author#1{\global\advance\@c@thr@ by 1
          \global\expandafter\edef\csname @thr@\romannumeral\@c@thr@\endcsname{#1}
          \global\expandafter\def\csname @ddrss@\romannumeral\@c@thr@\endcsname{}
          \global\expandafter\def\csname @scndddrss@\romannumeral\@c@thr@\endcsname{}
          \global\expandafter\def\csname @smddrss@\romannumeral\@c@thr@\endcsname{}}%\author
%    Elaborating the two author lists (First page and heading)
\newcount\@y\newcount\@x
\def\@cnjctn{\ifnum\the\@c@thr@=1\null\else{{\ and\ }}\fi}
\def\@@th@rs{\@x=0\global\@y=\@c@thr@\global\advance\@y by -1
\loop\advance\@x by 1
\ifnum\the\@x<\the\@y\csname  @thr@\romannumeral\@x\endcsname\ignorespaces
            ${}^{\csname  @ddrss@\romannumeral\@x\endcsname
                 \csname  @smddrss@\romannumeral\@x\endcsname
                 \csname  @scndddrss@\romannumeral\@x\endcsname}$,
\repeat
\csname  @thr@\romannumeral\@y\endcsname\ignorespaces
            ${}^{\csname  @ddrss@\romannumeral\@y\endcsname
                 \csname  @smddrss@\romannumeral\@y\endcsname
                 \csname  @scndddrss@\romannumeral\@y\endcsname}$\@cnjctn
\csname  @thr@\romannumeral\@c@thr@\endcsname\ignorespaces
            ${}^{\csname  @ddrss@\romannumeral\@c@thr@\endcsname
                 \csname  @smddrss@\romannumeral\@c@thr@\endcsname
                 \csname  @scndddrss@\romannumeral\@c@thr@\endcsname}$}%\@@th@rs
\def\@qq#1#2{\vrule height#1 depth#2 width0pt}
\def\@hspc{\kern1pt}
\def\@hspcc{\kern2pt}
\def\@hspccc{\kern3pt}
\def\@hspn{\kern-1.5pt}
\newcount\@c@ddrss@\newbox\@b@ddrss@
%
\def\@dd@ddrss@#1{%
\global\setbox51=\vbox{\advance\hsize by-12pt\unvbox\@b@ddrss@
  \vtop{\footnotesize\noindent{${}^{\the\@c@ddrss@}$\ }\@qq{10pt}{0pt}\textrm{#1}}}
  \global\setbox\@b@ddrss@=\vbox{\unvbox51}}%
%
\def\@dd@scndddrss@#1{%
\global\setbox51=\vbox{\advance\hsize by-12pt\unvbox\@b@ddrss@
  \vtop{\footnotesize\noindent{${}^{\the\@c@ddrss@}$\ }\@qq{10pt}{0pt}\textrm{#1}}}
  \global\setbox\@b@ddrss@=\vbox{\unvbox51}}
%
\def\address#1{\global\advance\@c@ddrss@ by 1\@dd@ddrss@{#1}
  \expandafter\edef\csname @ddrss@\romannumeral\@c@thr@\endcsname{\@hspc\number\@c@ddrss@}}
\def\secondaddress#1{\global\advance\@c@ddrss@ by 1\@dd@ddrss@{#1}
  \expandafter\edef\csname @scndddrss@\romannumeral\@c@thr@\endcsname%
{,\@hspcc\number\@c@ddrss@}}%\secondaddress#1
\def\sameaddress#1{\expandafter\edef\csname @smddrss@\romannumeral\@c@thr@\endcsname{\@hspc{}#1}}
%
%-D- first page of the paper : \@setbegin
\def\ps@plain{
  \def\@evenhead{}%
  \def\@oddhead{}%
  \def\@evenfoot{}%
  \def\@oddfoot{\hbox to\hsize{\hfill\lower-5pt\hbox{\fontsize{8}{10}\selectfont\rmfamily
\copyright{} Soci\'et\'e Francaise d'Astronomie et d'Astrophysique (SF2A) 2016}}}
}
%
%-E- footnote and thanks
\addtolength{\skip\footins}{12pt}
\def\thanks#1{\protected@xdef\@thanks{\@thanks
        \protect\footnotetext[0]{#1}}%
}
%-F- figures
\renewcommand\figurename{Fig.}
\renewcommand\tablename{Table}
\def\fnum@figure{\figurename~\thefigure}
\def\fnum@table{\tablename~\thetable}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{{\small{\bf #1.}~#2}}%
  \ifdim \wd\@tempboxa >\hsize
    {\small{\bf #1.}~#2}\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
%-G- abstract
\renewenvironment{abstract}{%
        \small\@pc\begin{center}%
        \quotation\noindent\ignorespaces{\bf \abstractname.}
        }{\endquotation\@pc\end{center}}
%-H- keywords
\newenvironment{keywords}{%
        \small\@pc\begin{center}%
        \quotation\noindent\ignorespaces{Keywords:~}
        }{\endquotation\@pc\end{center}}
%%%
%%-----------------------------------------------------------------------------
%% --> SECTIONING MATERIAL
%%-----------------------------------------------------------------------------
\newdimen\bibindent
\setlength\bibindent{\parindent}
\def\@biblabel#1{}
\def\@cite#1#2{{#1\if@tempswa , #2\fi}}

\renewenvironment{thebibliography}[1]
     {\section*{\refname}
        \small 
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
           \leftmargin\bibindent
           \setlength{\itemindent}{-\leftmargin}
           \setlength{\parsep}{0pt}
           \@openbib@code
           \usecounter{enumiv}%
           \let\p@enumiv\@empty
           \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\let\JFMthebib\thebibliography
\let\JFMendthebib\endthebibliography

\renewcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
%
\renewcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \sffamily\bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
%
\renewcommand{\section}{\@startsection%
{section}%
{1}%
{\z@ }%
{-3.5ex \@plus -1ex \@minus -.2ex}%
{2.3ex \@plus .2ex}%
{\fontsize{10}{12}\selectfont\sffamily\bfseries}}
%
\renewcommand{\subsection}{\@startsection%
{subsection}
{2}
{\z@ }
{-3.25ex\@plus -1ex \@minus -.2ex}
{1.5ex \@plus .2ex}
{\fontsize{10}{13}\selectfont\sffamily\itshape}}
%
\renewcommand{\subsubsection}{\@startsection%
{subsubsection}
{3}{\z@ }
{-3.25ex\@plus -1ex \@minus -.2ex}
{1.5ex \@plus .2ex}
{\sffamily\normalfont}}
%
\@addtoreset{equation}{section}
\renewcommand{\theequation}{\arabic{section}.\arabic{equation}}
\def\Appendix{\section*{Appendix}
\renewcommand{\theequation}{\Alph{section}.\arabic{equation}}
  \renewcommand{\thesubsection}{\Alph{section}.\arabic{subsection}}\appendix\null}
%%-----------------------------------------------------------------------------
%% --> METRIC DATA
%%-----------------------------------------------------------------------------
\setlength{\paperwidth}{210mm}
\setlength{\paperheight}{297mm}

\topmargin -0.5cm  
\setlength\oddsidemargin {-2mm}\setlength\evensidemargin{-2mm}
\setlength{\textwidth}{170mm}\setlength{\textheight}{235mm}
\headheight=12.45pt
%- espaces
\def\@pc{\vskip2mm plus 0mm minus1mm}
\def\@spc{\vskip8mm plus 0mm minus3mm}
\def\@spcx{\vskip13mm plus 0mm minus3mm}
\def\@spcc{\vskip32mm plus 0mm minus 3mm}
\def\@spcont{\vspace*{-42pt}}
\def\@boat#1#2{\setbox50=\null\ht50=#1\dp50=#2\box50}
%%-----------------------------------------------------------------------------
%% --> MESSAGES
%%-----------------------------------------------------------------------------
\def\@MSSG@TtrGlbl{The title will be set by the publisher.}
\def\@MSSG@rnnghd{Give a shorter title using $\backslash${\tt runningtitle}}
%%-----------------------------------------------------------------------------
%% --> HEAD and FOOT
%%-----------------------------------------------------------------------------
\def\ps@headings{\ps@empty
  \def\@evenhead{\normalfont\normalsize{\normalsize\thepage}\hfil \rightmark{}{}\hfil}%
  \def\@oddhead{\normalfont\hfil\ifvoid\@b@rnngttl
         \setbox60=\hbox{\normalsize\normalfont\@thr@i\ifnum\the\@c@thr@=1\relax
                   \else\ifnum\the\@c@thr@=2,\@hspccc\@thr@ii\else{\@hspccc\etal}\fi\fi:\@hspccc\leftmark{}{}}
         \ifdim\wd60>0.9\hsize{\@MSSG@rnnghd}\else\box60\fi%(\ifdim\wd60>0.9\hsize)
                                \else\@MSSG@rnnghd\fi% (\ifvoid\@b@rnngttl)
  \hfil{\normalsize{\thepage}}}%\@oddhead
  \let\@mkboth\markboth}
%
\def\ps@myheadings{\ps@headings\let\@mkboth\@gobbletwo}
\pagestyle{headings}
%%-----------------------------------------------------------------------------
%% --> FLOAT
%%-----------------------------------------------------------------------------
\setcounter{topnumber}{2}
\renewcommand\topfraction{.95}
\setcounter{bottomnumber}{1}
\renewcommand\bottomfraction{.7}
\setcounter{totalnumber}{3}
\renewcommand\textfraction{.05}
\renewcommand\floatpagefraction{.95}
%%-----------------------------------------------------------------------------
%% --> ACKNOWLEDGEMENTS
%%-----------------------------------------------------------------------------
% Define `acknowledgements' environment
\def\acknowledgements{\par\addvspace{17pt}\footnotesize\rmfamily\trivlist\item[]}
\def\endacknowledgements{\endtrivlist\addvspace{6pt}}
%
%-->  LATIN ABBREVIATIONS
%................ ................ ................ ................  ................
\def\cf{\emph{cf.\/}}\def\ie{\emph{i.e.\/}}\def\etc{etc.}
\def\apriori{\emph{a priori\/}}\def\afortiori{\emph{a fortiori\/}}
\def\loccit{\emph{a loc. cit.\/}}\def\etal{\emph{et al.\/}}
\def\vg{\emph{v.g.\/}}
%
\normalsize
\AtBeginDocument{%
  \@ifpackageloaded{natbib}
    {\bibpunct{(}{)}{;}{a}{}{,}%
     \global\let\thebibliography\JFMthebib
     \global\let\endthebibliography\JFMendthebib
}
    {\relax}%
}
%%---------------------------------------------------------------------
% Bibliography and bibfile
%
\def\aj{AJ}%
          % Astronomical Journal
\def\actaa{Acta Astron.}%
          % Acta Astronomica
\def\araa{ARA\&A}%
          % Annual Review of Astron and Astrophys
\def\apj{ApJ}%
          % Astrophysical Journal
\def\apjl{ApJ}%
          % Astrophysical Journal, Letters
\def\apjs{ApJS}%
          % Astrophysical Journal, Supplement
\def\ao{Appl.~Opt.}%
          % Applied Optics
\def\apss{Ap\&SS}%
          % Astrophysics and Space Science
\def\aap{A\&A}%
          % Astronomy and Astrophysics
\def\aapr{A\&A~Rev.}%
          % Astronomy and Astrophysics Reviews
\def\aaps{A\&AS}%
          % Astronomy and Astrophysics, Supplement
\def\azh{AZh}%
          % Astronomicheskii Zhurnal
\def\baas{BAAS}%
          % Bulletin of the AAS
\def\bac{Bull. astr. Inst. Czechosl.}%
          % Bulletin of the Astronomical Institutes of Czechoslovakia 
\def\caa{Chinese Astron. Astrophys.}%
          % Chinese Astronomy and Astrophysics
\def\cjaa{Chinese J. Astron. Astrophys.}%
          % Chinese Journal of Astronomy and Astrophysics
\def\icarus{Icarus}%
          % Icarus
\def\jcap{J. Cosmology Astropart. Phys.}%
          % Journal of Cosmology and Astroparticle Physics
\def\jrasc{JRASC}%
          % Journal of the RAS of Canada
\def\mnras{MNRAS}%
          % Monthly Notices of the RAS
\def\memras{MmRAS}%
          % Memoirs of the RAS
\def\na{New A}%
          % New Astronomy
\def\nar{New A Rev.}%
          % New Astronomy Review
\def\pasa{PASA}%
          % Publications of the Astron. Soc. of Australia
\def\pra{Phys.~Rev.~A}%
          % Physical Review A: General Physics
\def\prb{Phys.~Rev.~B}%
          % Physical Review B: Solid State
\def\prc{Phys.~Rev.~C}%
          % Physical Review C
\def\prd{Phys.~Rev.~D}%
          % Physical Review D
\def\pre{Phys.~Rev.~E}%
          % Physical Review E
\def\prl{Phys.~Rev.~Lett.}%
          % Physical Review Letters
\def\pasp{PASP}%
          % Publications of the ASP
\def\pasj{PASJ}%
          % Publications of the ASJ
\def\qjras{QJRAS}%
          % Quarterly Journal of the RAS
\def\rmxaa{Rev. Mexicana Astron. Astrofis.}%
          % Revista Mexicana de Astronomia y Astrofisica
\def\skytel{S\&T}%
          % Sky and Telescope
\def\solphys{Sol.~Phys.}%
          % Solar Physics
\def\sovast{Soviet~Ast.}%
          % Soviet Astronomy
\def\ssr{Space~Sci.~Rev.}%
          % Space Science Reviews
\def\zap{ZAp}%
          % Zeitschrift fuer Astrophysik
\def\nat{Nature}%
          % Nature
\def\iaucirc{IAU~Circ.}%
          % IAU Cirulars
\def\aplett{Astrophys.~Lett.}%
          % Astrophysics Letters
\def\apspr{Astrophys.~Space~Phys.~Res.}%
          % Astrophysics Space Physics Research
\def\bain{Bull.~Astron.~Inst.~Netherlands}%
          % Bulletin Astronomical Institute of the Netherlands
\def\fcp{Fund.~Cosmic~Phys.}%
          % Fundamental Cosmic Physics
\def\gca{Geochim.~Cosmochim.~Acta}%
          % Geochimica Cosmochimica Acta
\def\grl{Geophys.~Res.~Lett.}%
          % Geophysics Research Letters
\def\jcp{J.~Chem.~Phys.}%
          % Journal of Chemical Physics
\def\jgr{J.~Geophys.~Res.}%
          % Journal of Geophysics Research
\def\jqsrt{J.~Quant.~Spec.~Radiat.~Transf.}%
          % Journal of Quantitiative Spectroscopy and Radiative Transfer
\def\memsai{Mem.~Soc.~Astron.~Italiana}%
          % Mem. Societa Astronomica Italiana
\def\nphysa{Nucl.~Phys.~A}%
          % Nuclear Physics A
\def\physrep{Phys.~Rep.}%
          % Physics Reports
\def\physscr{Phys.~Scr}%
          % Physica Scripta
\def\planss{Planet.~Space~Sci.}%
          % Planetary Space Science
\def\procspie{Proc.~SPIE}%
          % Proceedings of the SPIE
\let\astap=\aap
\let\apjlett=\apjl
\let\apjsupp=\apjs
\let\applopt=\ao
%%-------------------------------------------------------------------------

\endinput

