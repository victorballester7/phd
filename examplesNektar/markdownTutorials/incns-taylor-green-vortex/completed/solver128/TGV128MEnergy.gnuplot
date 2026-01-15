set term pdf
set output "TGV128MEnergy.pdf"
set log xy
set format x "10^{%L}"
set format y "10^{%L}"
set xlabel "k"
set ylabel "E_k"
plot [1:32] [1e-9:1] for [i=1:10] "TGV128MEnergy.mdl" using 2:3 every 1::(20*i)*32::(20*i+1)*32-1 with lines title "tau = ".(2*i)