# OrCAS Target Selection script
![alt text](https://github.com/iancrossfield/OrCAS_target_selection/blob/main/graphics/orcas_logo.png?raw=true)

Target selection and prioritization for OrCAS Survey (cf. Crossfield et al. 2024)

Steps to use:
1) Install the KPF Exposure Time Calculator from https://github.com/California-Planet-Search/KPF-etc
2) Download this repository into your favorite working directory
3) Run orcas_selection.py at the python prompt, IPython terminal, or Jupyter Notebook (as you prefer)
4) Get an output printed to the screen that is something like:

<pre>
          toi   pl_rade  pl_orbper    Vmag      k_rv        TSM   kpf_tottime  SUR  
291   1247.01  2.167310  15.923524   9.080  1.189575  33.501752  8.880000e+03 3.843099e+00  
949   1799.01  1.627715   7.093853   8.980  0.936103  35.371515  8.430000e+03 3.297758e+00  
736   1630.01  2.411880  12.055754   9.624  1.938917  42.407983  1.212000e+04 3.113258e+00  
806   1691.01  2.927280  16.736914  10.134  2.287606  57.343045  1.779000e+04 2.260541e+00  
836   1716.01  2.772520   8.082331   9.406  2.258838  50.845183  1.230000e+04 2.203562e+00  
844   1723.01  2.999970  13.726522   9.664  2.427055  47.290155  1.413000e+04 1.991311e+00  
912   1772.01  2.833771   8.054221   9.955  2.739034  45.797036  1.530000e+04 1.503284e+00  
5754  6054.01  2.527773   7.490327   8.020  1.607805  46.386687  2.174784e+04 1.261493e+00  
874   1744.01  2.633540  22.341674   9.543  1.531183  26.285733  1.864862e+04 1.257229e+00  
....
</pre>

Where the column "SUR" is the so-called "Sub-Neptune Underrepresentation Rate."
