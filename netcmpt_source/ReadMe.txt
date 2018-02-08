NetCmpt ReadMe file

This program is supposed to automate the steps from input (species-specific lists of EC numbers) to output (Competition_matrix.txt).

Content
-----------------------
The NetCmpt package contains:
- example_NetCmpt_input.txt: an example input file 
- reaction_parsed.txt: a generic, cross-species metabolic network
- enzyme_lables.txt,reaction_labeles.txt,ec_reac_mapping.txt,reac_ec_mapping.txt,compounds_labels.txt, compund_labels.txt: index files 
- NetSeed package (NetSeed.pm): NetSeed package (http://elbo.gs.washington.edu/software.html) is implanted in the NetCmpt package.
				   It is used for calculating the species-specific metabolic environment.
- biomass_vector.txt: a list of generic target metabolites
- NetCmpt.exe: the NetCpmt code computing inter-species competition
- runNetcmpt.pl: the wrapper code for the running NetCmpt package



Input file:

In the input file, each row represents a species. The first entry represents an identifier, followed by the species-specific 
EC content where entries are separated by space. 

Output file:

The Output file is a matrix describing all pairwise computed EMO values (Effective Metabolic Overlap). EMO score is not symmetrical and [i,j] and [j,i] 
donate different values. The score in each cell describes the effect of the 

Output file:

The Output file is a matrix describing all pairwise computed EMO values. EMO score is not symmetrical 
and [i,j] and [j,i] donate different values. The score in each  cell describes the effect of the 
column species on the row species. That is, in the example here species 1 is not affected by species 
2 (EMO score 0 at [1,2]) but species 2 is affected by species 1 (EMO score 1 at [2,1]).

********************************************************
How do I get the EC content list?
For non-sequenced species, the EC content can be constructed according to BLAST results
and annotation tools such as BLAST2GO
For sequences species, the EC content can be retrieved from sites such as the DOE JGI site
(genome.jgi.doe.gov). Following choosing a specific species (e.g., 
http://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=642555120)
look at enzymes (under genome statistics). This leads to species-specific lists of EC numbers
(http://img.jgi.doe.gov/cgi-bin/w/main.cgi?section=TaxonDetail&page=enzymes&taxon_oid=642555120)
********************************************************

Running NetCmpt
--------------------------
###To run the program use:"perl RunNetcmpt.pl [input_file]" 
###A typical command might look like:" perl RunNetcmpt.pl example_NetCmpt_input.txt"".

this would:
1. For each input sequence (a row == species) constructs the species-specific metabolic network, species-specific metabolic environment and species-specific
list of target metabolites. 
2. For each pairwise combination computes pair-specific environment and EMO score 

To run this programs you will need:
1. A working installation of PERL (freely available at http://www.perl.org/get.html).
2. read/write permissions for the folder containing the input file.

Contact: shiri.freilich@gmail.com

