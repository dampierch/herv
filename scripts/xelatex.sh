## this script is meant to be run from snakemake pipeline
## it will compile a complete pdf including figures, tables, and refs
## it will also compile a rough docx file for sharing with collaborators
## sample useage: snakemake manuscript

fig1=$1
fig2=$2
fig3=$3
tab1=$4
tab2=$5
spf1=$6
spf2=$7
spf3=$8
spf4=$9
spf5=${10}
spf6=${11}
spf7=${12}
spt1=${13}
spt2=${14}
spt3=${15}

cp $fig1 ../manuscript/figs/fig1.pdf
cp $fig2 ../manuscript/figs/fig2.pdf
cp $fig3 ../manuscript/figs/fig3.pdf
cp $tab1 ../manuscript/tabs/tab1.tex
cp $tab2 ../manuscript/tabs/tab2.tex
cp $spf1 ../manuscript/figs/spf1.pdf
cp $spf2 ../manuscript/figs/spf2.pdf
cp $spf3 ../manuscript/figs/spf3.pdf
cp $spf4 ../manuscript/figs/spf4.pdf
cp $spf5 ../manuscript/figs/spf5.pdf
cp $spf6 ../manuscript/figs/spf6.pdf
cp $spf7 ../manuscript/figs/spf7.pdf
cp $spt1 ../manuscript/tabs/spt1.tex
cp $spt2 ../manuscript/tabs/spt2.tex
cp $spt3 ../manuscript/tabs/spt3.tex

printf '%s' "Adding texlive to \$PATH\n"
export PATH=$PATH:/home/chd5n/apps/texlive/2020/bin/x86_64-linux; xelatex -v
export PATH=$PATH:/home/chd5n/apps/pandoc-2.10/bin; pandoc -v

cd ../manuscript/
rm main.{aux,bbl,blg,log,pdf}
rm supp.{aux,bbl,blg,log,pdf}

xelatex main
bibtex main
xelatex main
xelatex main

xelatex supp
bibtex supp
xelatex supp
xelatex supp

xelatex cover

pandoc \
    --standalone \
    --from=latex \
    --to=docx \
    --reference-doc=docx/custom-reference.docx \
    --bibliography=main.bib \
    --csl=ama.csl \
    --output=docx/main.docx \
    main.tex

pandoc \
    --standalone \
    --from=latex \
    --to=docx \
    --reference-doc=docx/custom-reference.docx \
    --bibliography=main.bib \
    --csl=ama.csl \
    --output=docx/supp.docx \
    supp.tex
