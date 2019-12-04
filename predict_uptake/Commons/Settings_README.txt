Uptake Model Settings READ ME:

USS list:  Contains the genome position and USS score for each USS that scores above a designated score threshold, in the genome whose uptake is being modeled.  The file name should include the genome (‘Fake’, ‘NP’, ‘GG’ or ‘Rd’) and the score threshold used (e.g. ‘Fake10kb_USS10’, ‘NP_USS9.5’).  Note that genomes bigger than 100 kb take a long time, especially with long fragment sizes (on a new Mac, ~30 min for 100 kb).

folder name: described the folder path in which files will be saved.

Genome used: This setting is used to calculate the genome length. Options are 'NP', 'GG', 'Rd', 'other'. If NP, GG or Rd option is chosen, then script assumed you used the strain sequences used in our manuscript. In which NP genome length is 1914386 bp, GG genome length is 1887046 bp and Rd length is 1831585 bp. If another genome is used then choose the 'other' option and set genome length at the 'Custom genome length'.

Custom genome length: as integer. Allows the user to choose the length of a customized genome with a different size from the ones described above. This setting will be used by the script only if the user has chosen 'other' in the 'Genome used' settings.     

Segment start: as integer. Position 1 is the start of the genome.

Segment end:  as integer. NP full length is 1914386 bp, GG full length is 1887046 bp. 

Uptake_table:  Specifies the uptake probability associated with each USS score.  The name should indicate the uptake function used, and the USS cutoff.  (e.g. ‘Linear_9.5’ or ‘Sigmoidal_10’)

Baseline_binding:  This is a single number less than 1.  ‘0’ means no specific binding when there is no USS. (e.g. 0.15)

Baseline_uptake:  This is a single number less than 1. ‘0’ means no DNA uptake when there is no USS.  (e.g. 0.15)

Length_of_end_overlap:  This is the number of USS positions to be duplicated at the end of the USS list or list segment being used.  A length of 10 is a reasonable number for even long DNA fragments.

Fragment_distribution_file:  This file gives the size distribution of the DNA whose uptake is being modeled.  The name should be as informative as possible (e.g. ‘Fake_50-500bp’ or ‘GG_long’).

Run_ID:  A file listing the settings will be generated at the start of the run, with this identifier.  The identifier can also be used for any data files saved from the run.  You might want to include the date in the identifier.
