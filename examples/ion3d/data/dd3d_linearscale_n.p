      # Gnuplot script file for plotting data in file "*.dat"
      # This file is called   dd3d.p
      # line colore:
      # 1:red; 2: green; 3:blue;4:purple;5:light blue; 6:yellow;7:black;8:orange;9:grey;0:black

      set   autoscale                        # scale axes automatically
      set key bmargin center horizontal      # put the legend in the bottem center; remove this line to put on right top
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set title "nekcem vs reference paper: linear scale"

      set xlabel "x (nm)"
      set ylabel "densities (cm^{-3})"
      set grid
      set xr [-5:8.5]
      set mxtics
      set mytics
      plot      "cn_paper.dat" using ($1-5):(10**($2+21)) title 'n paper' with lines linetype 1 linecolor 0 linewidth 2, \
               		"cn_nekcem_bdf1.dat" using ($1-5):(10**($2+21)) title 'n nek coarse mesh' with lines linetype 1 linecolor 3 linewidth 2,  \
		"cn_1k_nekcem_bdf1_new.dat" using ($1-5):(10**($2+21)) title 'n nek fine mesh' with lines linetype 1 linecolor 5 linewidth 2
		
	set terminal postscript color portrait dashed enhanced 'Times-Roman'
	set output 'plot_linear_n.eps'
	set size 1,0.5
        replot

