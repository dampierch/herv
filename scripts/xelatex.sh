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
spt1=$8
spt2=$9

cp $fig1 ../manuscript/figs/fig1.pdf
cp $fig2 ../manuscript/figs/fig2.pdf
cp $fig3 ../manuscript/figs/fig3.pdf
cp $tab1 ../manuscript/tabs/tab1.tex
cp $tab2 ../manuscript/tabs/tab2.tex
cp $spf1 ../manuscript/figs/spf1.pdf
cp $spf2 ../manuscript/figs/spf2.pdf
cp $spt1 ../manuscript/tabs/spt1.tex
cp $spt2 ../manuscript/tabs/spt2.tex

printf '%s' "Adding texlive to \$PATH\n"
export PATH=$PATH:/home/chd5n/apps/texlive/2020/bin/x86_64-linux; xelatex -v
export PATH=$PATH:/home/chd5n/apps/pandoc-2.10/bin; pandoc -v

cd ../manuscript/
rm main.{aux,bbl,blg,log,pdf}
rm ._main.pdf
rm supp.{aux,bbl,blg,log,pdf}
rm ._supp.pdf

xelatex main
bibtex main
xelatex main
xelatex main

xelatex supp
bibtex supp
xelatex supp
xelatex supp

pandoc --standalone --output=docx/main.docx --from=latex --to=docx main.tex
