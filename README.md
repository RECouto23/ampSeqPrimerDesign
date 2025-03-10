# ampSeqPrimerDesign

A primer design script that relies on the Primer3 python package to design amplicon sequence primers adapted with the SP5 and SP7 handle sequences. Inputs should be supplied as the intended target site of the amplicon in BED format, with column headers of "chr", "start", "stop" and "name". Script formats Primer3 outputs into an excel sheet which is stored in the specified directory. Use the -h command line argument to see the other arguments for the script.
