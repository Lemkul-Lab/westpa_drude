cp west.h5 old_west.h5
w_pdist -W old_west.h5
plothist average pdist.h5 0 1
mv hist.pdf 2dhist.pdf
plothist evolution pdist.h5 0
mv hist.pdf 1dhist0.pdf
plothist evolution pdist.h5 2
mv hist.pdf 1dhist1.pdf
w_pdist -W old_west.h5 --construct-dataset construct_funcs.construct_auxdata1 
plothist average pdist.h5 0 1
mv hist.pdf 2dhistaux.pdf
plothist average pdist.h5 1 2
mv hist.pdf 2dhistaux2.pdf
plothist evolution pdist.h5 0
mv hist.pdf 1daux0.pdf
plothist evolution pdist.h5 1
mv hist.pdf 1daux1.pdf
plothist evolution pdist.h5 2
mv hist.pdf 1daux2.pdf
