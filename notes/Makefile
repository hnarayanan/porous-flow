# Modified from the simple build system I used for managing my
# doctoral dissertation.
# (c) Harish Narayanan, 2007 -- 2010

# Fundamental variables:

# A versioning scheme to keep track of progress
VERSION = 0.0.1
BASENAME = adaptive-porous-flow

# Change these to your heart's content.
LATEX   = latex
BIBTEX  = bibtex
DVIPS   = dvips
PS2PDF  = ps2pdf14 #On Mac OS X, you might need to use pstopdf instead.

SOURCES = Makefile ${BASENAME}.tex ${BASENAME}.bib ${BASENAME}.bst \
	  auxiliary/*.tex sections/*.tex \
          images/*/*.eps

# Subtle changes to the command-line flags below can have significant
# impact on the quality of generated documents. Comment-out the ps2pdf
# flags below if on Mac OS X.
LFLAGS  = 					#latex flags
BFLAGS  = -terse				#bibtex flags
DFLAGS  = -q -z -Z -t letter -Ppdf -G0 -o ${BASENAME}.ps  #dvips flags
PFLAGS  = -dPDFSETTINGS=/prepress		#ps2pdf flags
ADPATH  = "/Library/Application Support/Adobe/Adobe PDF/Settings/"
ADFILE  = "High Quality Print.joboptions"	#distiller settings file
#ADFILE  = "Press Quality.joboptions"

# Human-readable targets:

# I find it easiest to use the following names while running make. For
# e.g., running 'make ${BASENAME}' (or just 'make') cleans the source
# distribution of cruft and generates a pristine PDF.

all:	  ${BASENAME}

compose:  ${BASENAME}.dvi

polish:   ${BASENAME}.pdf

preview:  ${BASENAME}.dvi
	  xdvi ${BASENAME}.dvi&

${BASENAME}:   pristine polish
	  make clean

distill:  pristine distiller
	  make clean

archive:  pristine
	  hg archive -ttbz2 ../older/${BASENAME}-${VERSION}.tar.bz2

# Actual drudgery:

# Shoo, here be dragons.

clean:
	@echo "Cleaning cruft:"
	rm -f *~ *.aux *.bbl *.blg *.brf *.dvi *.loa *.lof *.log *.lot *.out \
	*.ps *.toc *.tmp *.spl
	(cd auxiliary; make clean)
	(cd sections; make clean)
	(cd images; make clean)

pristine: clean
	@echo "Removing older output files:"
	rm -f *.pdf

${BASENAME}.dvi: ${SOURCES}
	@echo "Creating the dvi file:"
	${LATEX}   ${LFLAGS} ${BASENAME} && ${BIBTEX}  ${BFLAGS} ${BASENAME} && \
	${LATEX}   ${LFLAGS} ${BASENAME} && ${LATEX}   ${LFLAGS} ${BASENAME}
	egrep --color=auto -i '(warn|error|full)' ${BASENAME}.log

${BASENAME}.ps: ${BASENAME}.dvi
	@echo "Creating the postscript file:"
	${DVIPS}   ${DFLAGS} ${BASENAME}.dvi
	@echo "Fix capitalisation to allow creation of elegant PDF bookmarks:"
	perl -pi -e 's/Title \(([A-Z])([A-Z].*)\)/Title (\1\L\2)/' ${BASENAME}.ps

${BASENAME}.pdf: ${BASENAME}.ps
	@echo   "Creating the PDF file:"
	${PS2PDF}  ${PFLAGS} ${BASENAME}.ps

distiller: ${BASENAME}.ps
	@echo	"Creating the PDF file using Acrobat Distiller:"
	osascript distiller.applescript ${BASENAME}.ps ${ADPATH}${ADFILE}