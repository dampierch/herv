## this script is experimental

fig1=$1
fig2=$2
tab1=$3
tab2=$4
spf1=$5
spt1=$6
spt2=$7

cp $fig1 ../manuscript/figs/fig1.pdf
cp $fig2 ../manuscript/figs/fig2.pdf
cp $tab1 ../manuscript/tabs/tab1.tex
cp $tab2 ../manuscript/tabs/tab2.tex
cp $spf1 ../manuscript/figs/spf1.pdf
cp $spt1 ../manuscript/tabs/spt1.tex
cp $spt2 ../manuscript/tabs/spt2.tex

printf '%s' "Adding texlive to \$PATH\n"
export PATH=$PATH:/home/chd5n/apps/texlive/2020/bin/x86_64-linux

cd ../manuscript/
rm main.{aux,bbl,blg,log,pdf}
rm ._main.pdf

xelatex main
bibtex main
xelatex main
xelatex main
