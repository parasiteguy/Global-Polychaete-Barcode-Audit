This repository contains the full reproducible workflow used in the manuscript Davinack - Where Are the Polychaete Barcodes? A Worldwide Audit Reveals Major Gaps in COI Coverage and Metadata Quality in GenBank.
The scripts in this repository: 
  - download and clean polychaete COI metadata from GenBank
  - reconcile names against an extant-only WoRMS backbone
  - compute COI coverage at species, genus and family levels
  - quantify metadata completeness (voucher-linked, identified-by fields)

All analyses were performed in Python 3.10 using pandas, Bipython, Numpy, matplotlib and seaborn packages. 
