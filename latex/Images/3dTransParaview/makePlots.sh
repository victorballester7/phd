# get all the .png files in pngs directory

pngFiles=$(ls pngs/*.png)

texSample="sample.tex"

for file in $pngFiles
do
    fileName=$(basename $file)

    echo "Processing $fileName - $file"

    # substitute the line
    # \def\image{pngs/d4w19_3d_l2crit_15.png}
    # in sample.tex with 
    # \def\image{$file}
    # and save it into pdfs directory with name $fileName.tex
    # Then compile with pdflatex and clean up the auxiliary files

    sed "s|\\\\def\\\\image{.*}|\\\\def\\\\image{../$file}|" $texSample > pdfs/${fileName%.png}.tex
    cd pdfs
    pdflatex ${fileName%.png}.tex
    rm ${fileName%.png}.aux ${fileName%.png}.log
    cd ..
done
