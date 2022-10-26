# A bunch of stuff for downloading and checking stuff on a bunch of ecoli genomes from enterobase
The scripts work together based on the assumption that stuff is being done in a particular project directory. The R scripts take this directory as a command line argument.
# Dependencies
We used [quast](http://bioinf.spbau.ru/quast) for checking basic statistics on the assemblies, [kma](https://bitbucket.org/genomicepidemiology/kma/src/master/) for checking matches to specific genes.
# Downloads
First, metadata containing assemblystats, cgmlst, wgmlst phylotypes and achtman mlst was downloaded off enterobase, in our case such that the the release date was at most 22-09-27.  
download_links.py was used for downloading assembly links.  
wrangle_ecoli_meta.R was used for combining the metadata files and filtering them to get a list of desired samples.  
download_fastas.py and generate_download.R  was used for downloading the desired subset of assemblies fitting filtering criteria in increments of 1000. This might take a few days.