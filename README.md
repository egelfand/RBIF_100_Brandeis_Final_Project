## RBIF_100_Brandeis_Final_Project
Final project for RBIF 100 Intro to Bioinformatics course from Brandeis

Emily Gelfand


### DESCRIPTION: 

Using RESTful APIs, this script collects Insulin protein data from the UniProt database and analyzes whether sequence length or 
long simple amino acid repeat length impact the number of interactions a given Insulin protein has with others. The script runs 
the following steps, sequentially:

- Step 1 – Pull down Insulin proteins and their interactions in batches from UniProt, then subset that database to just 
                the proteins with at least one interaction. For those proteins, it then calls UniProt's API again to gather gene 
                name and sequence information before saving the results as a final database file.
- Step 2 – Calculates sequence lengths and number of interactions for each protein, then creates and saves a scatterplot of 
                the interactions between the two. 
- Step 3 – For each protein, calculates the longest chain of repeated amino acid polypeptides and saves the maximum length 
                value. Next, calculates the mean and standard deviation of the maximum amino acid repeat (AAR) values to label 
                all proteins with significantly long maximum AARs, and creates and saves a boxplot of the number of protein 
                interactions vs maximum AAR lengths, with the significantly long maximum AAR values highlighted. 
- Step 4 - Finally, labels the subset of proteins that have both the sequence lengths within which the peak number of 
                interactions fall in addition to the significantly high maximum AAR lengths. Creates and saves a 3d interactive 
                plot of the sequence lengths vs maximum AAR lengths vs number of interactions with this intersection subset of 
                proteins highlighted.

### EXTRAS:
- Script and results can also be downloaded from https://github.com/egelfand/RBIF_100_Brandeis_Final_Project.git
- Summary report PDF containing description of analyses and embedded images can be found in the results folder, titled 
        "final_project_summary_report_emilyg.pdf"

### HOW TO RUN:
- Installation requirements:
    - pandas
    - requests
    - re
    - seaborn
    - matplotlib
    - requests
    - plotly
    - kaleido
- To run, type the following into your CLI (save_path parameter optional).

      /bin/python3 final_project.py save_path='./results/'
  
- Optional Input Parameters: 
    - save_path: location to which files should be saved. default save path (no parameters passed in) is './results/'

OUTPUTS:
- /results/data.csv: Insulin proteins and their interactions, downloaded from the UniProt database, filtered to only include proteins 
                with at least one interaction. Sequences and gene names have been added, also from the UniProt database.
- /results/figure1.png: 2D scatterplot of the number of protein interactions vs the sequence lengths. Band of sequences within which the 
                    peak number of interactions fall highlighted for emphasis. 
- /results/figure2.png: Boxplot of the number of protein interactions for each of the maximum AAR lengths. Significantly long maximum AAR 
                    sequences highlighted separately. 
- /results/figure3.html: 3D interactive scatterplot of the sequence lengths, maximum AAR lengths, and number of interactions for each protein. 
                    Proteins with both the sequence lengths within which the peak number of interactions fall as well as significantly 
                    high maximum AAR lengths are highlighted for emphasis. 
