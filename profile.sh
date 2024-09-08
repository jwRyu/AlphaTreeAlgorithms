sudo perf record -g ./AlphaTree 
sudo perf report
gprof ./AlphaTree gmon.out | gprof2dot | dot -Tpng -o output.png
gprof ./AlphaTree gmon.out > output.txt
