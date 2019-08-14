import sys
import time
import scipy.sparse as sp
import scrublet as scr
 
#- Load counts (genes x cells) and reshape
cnts_file = str(sys.argv[1])
cnts      = sp.load_npz(cnts_file)
cnts      = cnts.transpose()
 
#- Perform doublet annotation
tme_s            = time.time()
scrub            = scr.Scrublet(cnts)
dbl_scs, dbl_cls = scrub.scrub_doublets(verbose=False)
tme_e            = time.time()
 
#- Report time (wall-clock, in seconds)
print(cnts_file,"\t",round(tme_e - tme_s))


 
 
 
  
 
