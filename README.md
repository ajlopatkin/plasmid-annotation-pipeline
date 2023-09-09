# plasmid-annotation-pipeline
Annotates completed plasmid sequence fasta files

**LEGEND:**
  1) Set-Up Computer
     - Required Packages
     - Computer Settings
  2) Script Information
     - Order
     - Scripts
  3) Running PlasAnn
  4) Error Handling

**SET-UP YOUR COMPUTER:**
  1) REQUIRED PACKAGES:
     - [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
     - Main environment for running PlasAnn (containing biopython, blast, goatools, cogclassifier, reportlab, requests)
       - The dependencies between packages can be difficult to solve; we recommend using Conda + Mamba to install this environment by running:
         ```
         conda create -n plasann_env -c conda-forge mamba
         conda activate plasann_env
         mamba install -c bioconda biopython blast goatools cogclassifier reportlab requests altair 
         ```
     - Prokka
       - Use either Conda or Conda + Mamba, i.e.
         ```
         conda create -n prokka_env -c conda-forge mamba
         conda activate prokka_env
         mamba install -c bioconda prokka
         ```
     - Abricate
       - Use either Conda or Conda + Mamba, i.e.
         ```
         conda create -n abricate_env -c conda-forge mamba
         conda activate abricate_env
         mamba install -c bioconda abricate
         ```
     - CD-HIT and PSI-CD-HIT
       - Folder containing programs for CD-HIT is already stored in the scripts folder

  2) COMPUTER SETTINGS:
  The full program does take some time to run due to how much data needs to be processed for each plasmid. To make sure that there are no issues running the code in its entirety, set up your computer so it doesn’t turn off due to inactivity. To do this please go into your computer's settings and set up your computer so it  doesn’t turn off your screen or turn on the screen saver when inactive.
  Otherwise, your computer will enter its inactive state and shut down the program, leading to connection errors.

**SCRIPT INFORMATION:**
  1) ORDER: 
    PlasAnn.py 
    ---bgbm.sh
    ------annotate.py
    ------abricate_dw.py
    ------abricate_test.py
    ------prokka_default.py
    ------combining_hit_info.py
    ------oric_blast.py
    ------orit_blast.py
    ------TnFinder_loc.py
    ------COG_identify.py
    ------prokka2go.py
    ------update_tsv.py
    ------abr_genes.py
    ------update_gff.py
    ------final_edits.py
    ------visual_map.py
    ---lr_file.py
    ---linear-regression.py 

**SCRIPTS:**
  1) PlasAnn.py
     - INPUT: A 0, 1, or 2 depending on what parts of the program need to be run (Will eventually need location information for the folder containing fasta files for each plasmid. This information is currently built into the pipeline but it will be changed before the final version is released). PlasAnn will also ask the user to select whether it wants PlasAnn to break ties between genes or use the user defined naming convention for genes for different INC groups. 
     - OUTPUT: TSV final annotation, plasmid map, and acquisition cost linear regression analysis for all plasmids 
     - DESCRIPTION: PlasAnn is the script a user must call to start annotating plasmids. It reads through my fasta folder, which contains a fasta file for each plasmid I want analysed, one fasta file at a time, taking the plasmid name from the fasta file and passing it on to bgbm.sh, a shell script calling each script in the pipeline in their respective order in order to annotate the given plasmid.
You will be prompted to enter either a 0, 1, or 2 depending on whether you want to run the entire pipeline, just the linear regression script, or run BLAST.

  2) bgbm.sh
     - INPUT: Plasmid name
     - OUTPUT: TSV final plasmid annotation and GFF final plasmid annotation and plasmid map
     - DESCRIPTION: bgbm.sh is a shell script that calls each subsequent script in the annotation and visualization pipelines in their respective order in order to annotate and visualize the given plasmid. It utilizes conda to help isolate each script that needs its own environment in order to run properly.
Order in which these scripts are called:
    1) annotate.py
    2) abricate_dw.py
    3) abricate_test.py
    4) prokka_default.py
    4) combining_hit_info.py
    5) oric_blast.py
    6) orit_blast.py
    7) TnFinder_loc.py
    8) COG_identify.py
    9) prokkka2GO.py
    10) update_tsv.py
    11) abr_genes.py
    12) update_gff.py
    13) final_edits.py
    14) visual_map.py

  3) annotate.py
     - INPUT: Plasmid name
     - OUTPUT: Plasmids folder, which will contain a folder for each plasmid. gb_match folder in the plasmid folder for the plasmid being processed (Ie. ./plasmids/plasmid_name/gb_match) that will contain the results of prokka and abricate. matches.csv=csv file that contains the accession ID’s of the initial BLAST gbk files that match our standards for including them in our final annotation of the given plasmid (located in plasmids folder). Gbk files for each accession ID that was an initial match to our standards for including them in our annotation. Located in 
     - DESCRIPTION: Through a BLAST Search and input fasta file, program finds matching accession ID's to given sequence and tests their reliability. Passing accession ID's are compiled into a CSV file. (1) Run BLAST (2) Read in BLAST file and parse through each accession ID, downloading their gbk files and determining if the given accession ID is fit for use (3) Compile good accession ID's into CSV file.

  4) abricate_dw.py 
     - INPUT: Plasmid name
     - OUTPUT: ab_results folder for each accession ID listed in matches.csv for that particular plasmid, and 2 text files per database (total of 8). NCBI, CARD, and     RESFINDER list matching antibiotic resistance genes. PLASMIDFINDER lists the matching incompatibility group for the given plasmid.
     - DESCRIPTION: Run abricate w/ncbi, card, resfinder, and plasmidfinder databases for each accession ID match 

  5) abricate_test.py 
      - INPUT: Plasmid name
      - OUTPUT: Edited final.tsv file that contains the final annotation of the given plasmid and final.gff file that contains the final annotation of the given plasmid + incompatibility group + fasta sequence
      - DESCRIPTION: After running abricate with 3 different antibiotic resistance databases for each accession ID, check that 2/3 databases found the same resistance genes, and write the names of these confirmed consensus resistance genes to the final annotation for the plasmid. Also add the incompatibility group to the final gff for the plasmid found by plasmidfinder

  6) prokka_default.py
     - INPUT: Plasmid name
     - OUTPUT: Prokka generated files in pk_results folder in the default folder in gb_match.
     - DESCRIPTION: Utilizing conda, this python script calls on the prokka_env to run prokka without providing a gbk file (ie. default mode). Once prokka is done running, it deactivates the prokka_env environment. 

  7) combining_hit_info.py
     - INPUT: Plasmid Name and run selection
     - OUTPUT: Edited the_final.tsv file
     - DESCRIPTION: Combine gene information from BLAST, Prokka, and Abricate to determine concensus gene at each locus. If the User selected to use naming convention files, the program will ask the user to break ties between concensus genes and save the selection to a INC-group specific fasta file, which the program parses to see if tie breakers (or naming conventions) have already been spefied. Otherwise, the program breaks ties by selecting the best gene (ie. which one appears most often, has a higher BLAST score, etc). It then writes the updated names to the_final.tsv file.
   
  8) oric_blast.py
     - INPUT: Plasmid Name
     - OUTPUT: oric.csv
     - DESCRIPTION: Blast fasta against oric.fna (obtained from BAKTA database) to obtain csv of matching oriC sites. 
    
  9) orit_blast.py
     - INPUT: Plasmid Name
     - OUTPUT: orit.csv
     - DESCRIPTION: Blast fasta against orit.fna (obtained from BAKTA database) to obtain csv of matching oriT sites. 

  10) TnFinder_loc.py
     - INPUT: Plasmid Name
     - OUTPUT: blastn file
     - DESCRIPTION: Run TnComp_finder.py using the plasmid fasta to get output blastn file of matching transposon sites.

  11) COG_identify.py
     - INPUT: Plasmid Name
     - OUTPUT: COGCassifier output .tsv files
     - DESCRIPTION: Obtain matching COG ID's for each gene using COGClassifier to obtain tsv files that specify what COGID was found for each locus.
 
  12) COG_identify.py
     - INPUT: Plasmid Name
     - OUTPUT: COGCassifier output .tsv files
     - DESCRIPTION: Obtain matching COG ID's for each gene using COGClassifier to obtain tsv files that specify what COGID was found for each locus.

  13) prokka2go.py
     - INPUT: Input location of gbk file created by prokka, output file location, and plasmid name
     - OUTPUT: prokka2go.txt
     - DESCRIPTION: Code is based upon prokka2kegg.py, which was built by Heyu Lin. Currently, prokka2go parses through the given gbk file, looking at one gene at a time. If a gene is identified, the script looks to see if it has a cited UNIPROT ID. If not, it will do a query search for that gene in UNIPROT to get a list of matching UNIPROT ID’s, pulling out the first reviewed UNIPROT match or first unreviewed match if there are no reviewed matches, saving the results for each gene onto a list. This list is then passed onto get_go_id, which accesses the given UNIPROT entry for that gene and parses through it to obtain the listed GO ID’s, prioritizing the GO ID’s labeled as being biological. If no GO ID’s are listed, the keywords listed under that UNIPROT entry are looked at for labeling. The matching GO ID for each gene is concatenated onto a list, which is passed onto the output function that writes the results to prokka2go.txt, which is located in the folder for the accession file that was just processed by prokka.

  14) abr_genes.py
     - INPUT: Plasmid name
     - OUTPUT: updated the_final.tsv
     - DESCRIPTION: Update final tsv by adding COG and GO functional information from prokka2go and COGClassifier

  15) abr_genes.py
     - INPUT: Plasmid name
     - OUTPUT: updated the_final.tsv
     - DESCRIPTION: Write ABR genes to their corresponding locus in the tsv file.

  16) update_gff.py
     - INPUT: Plasmid name
     - OUTPUT: updated the_final.tsv
     - DESCRIPTION: Update gff with new abr gene and combining_hit_info information.

  17) final_edits.py
     - INPUT: Plasmid name
     - OUTPUT: updated the_final.tsv
     - DESCRIPTION: Update oriC, oriT, and TnFinder information. 

  18) visual_map.py
     - INPUT: Plasmid name
     - OUTPUT: PDF of plasmid map for given plasmid
     - DESCRIPTION: Using GenomeDiagram, visual_map.py parses through the final annotation for the given plasmid to create a plasmid map of the specified plasmid. Each gene on the plasmid is labeled to fit one of the following categories of interest: Conjugation, Antibiotic Resistance, Stress-Response, Toxin, Metabolism, Other, Unknown (More will be added, such as Replication, Background/Backbone, and Antitoxin)

  19) lr_file.py
      - INPUT: File containing acquisition cost information for all or some of the plasmids analysed previously in the script
      - OUTPUT: Linear Regression file, which will be used in the R-script to run a linear regression analysis on the plasmids and the acquisition cost data, with the following information: 1) Binary plasmid gene presence information for each plasmid 2) Respective number of gene functions for each function 3) Acquisition cost data 4) Plasmid data (ie. source, etc.)
      - DESCRIPTION: Using PSI CD-Hit, I create a list of unique homologous genes found across all of the plasmids analysed for the purposes of adding binary information about gene presence for each plasmid in the final linear regression file. I append acquisition cost data, plasmid data, gene presence, and gene function data for each plasmid to this file so the file contains a table containing data for each plasmid in each corresponding row.  

  20) linear-regression.py
      - INPUT: PlasmidRegression.csv created by lr_file.py
      - OUTPUT: Print linear regression results + results.txt file containing linear regression output
      - DESCRIPTION: Print linear regression results based on gene presence data, plasmid incompatibility group, gene groups, and acquisition cost data.

**RUNNING PLASANN:**
1) Download the folder uploaded on GitHub onto your computer. It contains all scripts involved in the program + folders needed at this current stage to run it
2) Install conda + abricate + prokka (see SET-UP YOUR COMPUTER)
3) Make sure all fasta files you want to analyse are in the fasta folder.
4) Make sure that for every plasmid/fasta file, you have a corresponding folder in the plasmids folder with the NCBI BLAST .xml file containing 1000 accession ID hits for that plasmid (will be done automatically by the program in the future).
5) On your Mac Terminal, use the cd command to access the “PlasAnn” folder (you can right click on the folder and click "Get Info" to see the file path for the folder to access it)
6) Run python PlasAnn.py: 
   - python PlasAnn.py
   - Once you run it, the code will prompt you to enter a binary number (0 or 1).
      1) Enter 1 : For running ONLY the linear regression pipeline. Do this ONLY if you have all the annotation data OR you ran the FULL plasmid annotation pipeline. 
      2) Enter 0 : Run the FULL pipeline (annotation + linear regression).
      3) Enter 2 : Run BLAST.
         - WARNING: Depending on how many plasmids you are analyzing, running the full pipeline will take a while. If you decide to do this, please go into your computer's settings and set up your computer so it doesn’t turn off your screen or turn on the screen saver when inactive. Do not turn off your computer while the program runs; this could lead to errors or could cause the program to stall!
5) Wait for the program to finish! Do check on it from time to time to make sure that it’s running and let me know if any issues/questions arise.

**ERROR HANDLING**
1) Prokka/Abricate
   - If any errors occur with Prokka or Abricate (ie. tbl2asn not running OR prokka is not outputting any files) terminate the program. You will now need to uninstall prokka and/or abricate and then reinstall prokka and/or abricate. To reinstall, you will need to use the .yml files in the scripts folder. Thus, make sure that in your terminal you are in the PlasAnn directory. Below is the code you need:
     1) Make sure you’re in the PlasAnn folder in your terminal
        - cd file_path/PlasAnn 
     2) For uninstalling
        - conda env remove -n prokka_env
        - conda env remove -n abricate_env
     3) For reinstalling
        - conda env create --file ./scripts/prokka_env.yml
        - conda env create --file ./scripts/abricate_env.yml
   - If the following does not work, email me ASAP (jc4919@barnard.edu). It is possible that another program may be needed to correctly install prokka or abricate. If the issue persists, check that your virus blocking software or other applications are not installing or blocking prokka/abricate.

2) Connection Issues
   - If you are having connection issues or the error message “ConnectionResetError: [Errno 54] Connection reset by peer” appears on your terminal, you will need to terminate the program and rerun it. I can also help you set up the program so you can start from where it left off so you do not have to rerun it.

3) Other
   - If at any point the program crashes or some of the scripts begin to fail, quit the program (either terminate or quit using Command-C until no scripts are running) and send me a message (jc4919@barnard.edu) about what script was failing with the error message. I will look over the issue, fix it, and send you an updated file of the faulty script. I can also help you set up the program so you can start from where it left off so you do not have to rerun it.
