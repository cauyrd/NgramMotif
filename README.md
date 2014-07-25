Introduction
------------
NgramMotif is implemented by Python and it is designed to detect motif enrichment and spatial signatures of TF ChIP-seq data.

Pre-installtalation
-------------------
A. Programming environment and NGS tools

* Python version 2.5 or later (http://www.python.org/)
* BEDtools (http://code.google.com/p/bedtools/)

B. Python packages:

* HTSeq (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)
* Biopython (http://biopython.org/wiki/Main_Page)
* Numpy (http://www.numpy.org/)
* Matplotlib (http://matplotlib.org/)

Running NgramMotif
------------------
#### Command-line usage:
    python NgramMotif.py -b <Bedfile> -m <MOTIFname> -r <reference_genome.fasta>[opts]
#### Options:
	-e <int>    :increase the BED entry by the same number base pairs in each direction (default 150)
	-w <int>    :window size for peak regions (default 1000)
 	-s <int>    :step for sliding windows of peak regions (default 100)
 	-f <int>    :unit of shift base pairs from motif start to peak summit (default 3)
	-n          :transcription factor name of input ChIP-seq, the generated plot will highlight this motif (default None)
	-g          :path of reference genome size (default hg19)
 	-h          :produce this menu
#### Example:
    python NgramMotif.py -b sample_data/AR_DHT_ChIPseq.PeakSummit.bed -m sample_data/motif.txt -n AR
Example data are provide in the directory **sample_data**

Output
------
1. *.enrich.png is the motif enrichment plot
2. *.space.png is motif spatial distribution plot

Note
-------------------
1. Instead of peak summit bed file, user can also provide the peak region file as input bed file and use -e 0 as extension distance.**If this case, the *.space.png file will not be produced.
2. Folder **genome** contains genome size file. Currenely only support **human**. User can provide their own genome file using **-g** to specify the path.
3. Folder **transfact** contains all human TFs PFM matrix files from TRANSFACT database. Users can provide their own PFM matrix file and update the **namelist.txt** file in this folder to add the new term.
    

Contact us
----------
Questions, suggestions, comments, etc?

Author: Rendong Yang

Send email to cauyrd@gmail.com

Referencing NgramMotif
----------------------
updating soon!
