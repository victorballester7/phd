set term pdf
set termoption enhanced
set output "TGV64Energy.pdf"
set xlabel "{/Symbol t}"
set ylabel "E_k"
set key off
plot 'TGV64Energy.eny' using 1:2 with lines