# Favipiravir
The code used for the paper entitled - The Mechanism of Resistance to Favipiravir in Influenza - can be found at github.com/pclangat/barcoded-flu-seq and the code used was in the branch ‘short-dg-amplicon-farm’.  This creates files ending in .barcode_filtered.fas.  The files which you can download in this branch are the barcode_filtered.fas files which have then been mapped to the reference and they end in ga.fasta.  These files were made using the map to reference function in Geneious v6.0.  If you don't have access to this I would recommend you use a different aligner for example msa which will download at the beginning of the program.

To make this work you need to have R.  Run favipiravir_mutationbias_figure2 from a directory with all the ga.fasta files in and the program should hopefully make the figure.

To run the program: favipiravir_mutationbias_figure2()

Please ask if you have any questions.

If you want to run code like this for your own project, I would recommend dms_tools2 as a complete package.  https://jbloomlab.github.io/dms_tools2/index.html
