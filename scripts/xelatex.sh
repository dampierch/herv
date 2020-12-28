## this script is experimental

fig1=$1
fig2=$2
tab1=$3

cp $fig1 ../manuscript/figs/fig1.pdf
cp $fig2 ../manuscript/figs/fig2.pdf
cp $tab1 ../manuscript/tabs/tab1.tex

printf '%s' "Adding texlive to \$PATH"
export PATH=$PATH:/home/chd5n/apps/texlive/2020/bin/x86_64-linux

cd ../manuscript/
rm main.{aux,bbl,blg,log,pdf}
rm ._main.pdf

xelatex main
bibtex main
xelatex main
xelatex main
