import sys
import time
import numpy as np
import scipy.sparse as sp
import doubletdetection as dd
 
# NOTE
#- drop warning about densifying and progress/iterations in shell
#  this will only print the file and the number of seconds:
#  $ python time_dblDetection.py /Users/dennis/tt_sparse.npz  2>/dev/null | tail -n 1
#  PS: probably put this in a shell script with a trap to clean up temp files.... 
 
#- Load counts (genes x cells) and reshape
cnts_file = str(sys.argv[1])
cnts      = sp.load_npz(cnts_file)
cnts      = cnts.transpose()
 
#- Perform doublet annotation
tme_s    = time.time()
clf      = dd.BoostClassifier()
res      = clf.fit(cnts)
dbl_cls  = res.predict()
dbl_scs  = np.nanmean(res.all_log_p_values_,axis=1)
tme_e    = time.time()
 
#- Report time (wall-clock, in seconds)
print(cnts_file,"\t",round(tme_e - tme_s))


 
 
 
  
 
