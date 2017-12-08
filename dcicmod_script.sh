
# Works on my linux machines at work, but the location of various python libraries and include files may differ

#!/bin/bash
gcc -I/usr/include/python2.7/ -fPIC -fopenmp -c cicdensmodule.c

gcc -shared -Xlinker -export-dynamic -I/usr/include/python2.7/ -L/usr/lib/python2.7 -lm cicdensmodule.o -o cicdens.so


