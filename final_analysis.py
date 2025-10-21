import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from requests.adapters import HTTPAdapter, Retry
from itertools import groupby
import plotly.express as px
from pathlib import Path
import requests, re, plotly

BASE_SAVEPATH = ''

def main(base_savepath = './results/'):
    """
    DESCRIPTION: Primary analysis function. Calls, in order, functions to collect and save dataset from UniProt, calculate sequence 
                    lengths and protein interactions for proteins in the dataset, calculate Amino Acid Repeat polypeptide chains and 
                    calculate their significance, and multiple plots to demonstrate the results. 
    INPUTS: 
        base_savepath: String file path to which all files should be saved. 
    OUTPUTS:
        None
    """
    # update base savepath with parameter
    global BASE_SAVEPATH
    BASE_SAVEPATH = base_savepath
    if not Path(BASE_SAVEPATH).exists():
        Path(BASE_SAVEPATH).mkdir(parents=True, exist_ok=True)
    # Skip downloading if data already downloaded, otherwise download dataset
    if not Path(f'{BASE_SAVEPATH}data.csv').exists():
        #collect dataset of all insulin proteins from uniprot
        insulin_batch_path = download_insulin_batch()
        # narrow those down to proteins with interactions, then collect their sequences and gene names
        with_interact = populate_proteins_with_interactions(insulin_batch_path)
    # Analyze sequence length vs number of interactions
    protein_df = count_lengths(with_interact)
    protein_df = count_interactions(protein_df)
    #plot results
    length_vs_interact(protein_df)
    # Calculate amino acid repeat lengths
    protein_df = get_max_AAR_lengths(protein_df)
    protein_df = sig_AAR_lengths(protein_df)
    # plot results
    plot_int_vs_aar(protein_df)
    # finally plot length vs interaction vs AAR length together
    protein_df = label_peak_seq_len_and_aar(protein_df)
    plot_triple_interaction(protein_df)
    

def get_next_link(re_next_link, headers):
    """
    DESCRIPTION: Get the next download link information from the passed in response headers
    INPUTS:
        re_next_link: regex pattern for finding the next link in the headers
        headers: list of response headers from previous batch of requests
    OUTPUTS:
        group: group of matching headers
    """
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(session, re_next_link, batch_url):
    """
    DESCRIPTION: calls each batch of URLs as specified by calling function and returns, in order, the contents of the result headers.
    INPUTS:
        re_next_link: regex pattern for finding the next link in the headers
        batch_url: url for the next batch of requests
    OUTPUTS:
        response: http response from querying the batch
        total: header contents
    """
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(re_next_link, response.headers)

def download_insulin_batch():
    """
    DESCRIPTION: From Uniprot, downloads full batch of Insulin and related proteins which interact with it, and saves a dataframe of 
                    the results.
    INPUTS:
        None
    OUTPUTS:
        result_path: path to which downloaded result dataframe was saved
    """
    #pattern for determining appropriate next link to select
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    # begins the batch download process
    session.mount("https://", HTTPAdapter(max_retries=retries))
    url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=Insulin%20AND%20%28reviewed%3Atrue%29&size=500'
    progress = 0
    result_path = f'{BASE_SAVEPATH}data.csv'
    with open(result_path, 'w') as f:
        for batch, total in get_batch(session, re_next_link, url):
            # calculates results from each batch using above functions 
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])
            print(f'{progress} / {total}')
    return result_path

def populate_proteins_with_interactions(insulin_batch_path):
    """
    DESCRIPTION: Filters the Insulin and related proteins dataframe saved from the insulin batch query function to only proteins with 
                    some interactions, and collects the sequences and gene names for each of those proteins. Overwrites the original 
                    saved dataframe with the updated version. 
    INPUTS:
        insulin_batch_path: String path to saved insulin and related proteins dataframe from batch query function
    OUTPUTS:
        with_interact: pandas dataframe of proteins, interactions, sequences, and gene names for Insulin dataset.
    """
    # filter the input dataframe to only proteins with at least one interaction
    insulin_df = pd.read_csv(insulin_batch_path, sep='\t')
    with_interact = insulin_df[insulin_df['Interacts with'].isna()==False]
    all_proteins = with_interact['Entry'].unique().tolist()

    # collects sequence and gene name information for each of the remaining proteins
    for i, prot in enumerate(all_proteins):
        # takes about 20 minutes to run
        print(f"downloading and adding {i} out of {len(all_proteins)}")
        url = f'https://rest.uniprot.org/uniprotkb/{prot}'
        protein_entry = requests.get(url).json()
        # skip those with no gene name
        try:
            gene_name = protein_entry['genes'][0]['geneName']['value']
        except KeyError:
            continue
        seq_val = protein_entry['sequence']['value']
        curr_idx = with_interact[with_interact['Entry']==prot].index.values.tolist()
        # adds data back to original dataframe
        with_interact.loc[curr_idx, ['sequence', 'gene_name']] = [seq_val, gene_name]
    # save results
    with_interact.to_csv(insulin_batch_path)
    return with_interact

def count_lengths(protein_df):
    """
    DESCRIPTION: Calculates the sequence lengths for each sequence in the insulin and interactions dataframe. Additionally, labels all
                sequences within sequence length range of peak interactions, specifically lengths of 500 to 1500 amino acids. 
    INPUTS:
        protein_df: pandas dataframe of proteins, interactions, sequences, and gene names for Insulin dataset.
    OUTPUTS:
        protein_df: same dataframe as input, but now with sequence length (seq_len) and peak interaction sequence length (peak_int_seq_len) columns added.
    """
    protein_df['seq_len'] = protein_df['sequence'].str.len()
    protein_df['peak_int_seq_len'] = (protein_df['seq_len'] >= 500)&(protein_df['seq_len'] <= 1500)
    return protein_df

def count_interactions(protein_df):
    """
    DESCRIPTION: Calculates the number of proteins with which a given protein interacts
    INPUTS:
        protein_df: pandas dataframe of proteins, interactions, sequences, gene names, sequence lengths, and peak intereaction sequence 
                    lengths for Insulin dataset.
    OUTPUTS:
        protein_df: Same dataframe as input, but now with column of interaction counts (num_interacts)
    """
    protein_df['num_interacts'] = protein_df["Interacts with"].str.count(";")
    return protein_df

def length_vs_interact(protein_df):   
    """
    DESCRIPTION: Plots and saves scatterplot of number of interactions vs sequence lengths for all proteins within protein interaction 
                    dataset. Additionally, highlights band of lengths with proteins that have greatest number of interactions 
                    (500-1500 amino acids long).
    INPUTS:
        protein_df:pandas dataframe of proteins, interactions, sequences, gene names, sequence lengths, peak intereaction sequence 
        lengths, and counts of number of interactions
    OUTPUTS:
        None
    """
    sns.scatterplot(data=protein_df, y='num_interacts', x='seq_len', s=5, alpha=0.5)
    #highlight area of lengths with greatest number of interactions
    plt.axvspan(xmin=500, 
            xmax=1500, 
            color='grey', 
            label = 'Peak interactions lengths', 
            alpha=0.25)
    # add labels
    plt.xlabel('Length of sequence')
    plt.ylabel('Number of interactions')
    plt.title('Sequence lengths vs number of gene interactions')
    plt.legend()
    plt.savefig(f'{BASE_SAVEPATH}/figure1.png')

def get_max_AAR_lengths(protein_df):
    """
    DESCRIPTION: Calculates longest polypeptide chain of repeated single Amino Acids across each protein sequence and logs those values 
                to the protein dataframe
    INPUTS:
        protein_df: pandas dataframe of proteins, interactions, sequences, gene names, sequence lengths, peak intereaction sequence 
                    lengths, and counts of number of interactions
    OUTPUTS:
        protein_df: same dataframe as input, but now with the length of the longest polypeptide chain of repeated single Amino Acids 
                    stored (max_AAR_length)
    """
    all_max_AARs = []
    for seq in protein_df['sequence'].values.tolist():
        try:
            groups = groupby(seq) # group each amino acid in the sequence and get a count
        except TypeError:
            all_max_AARs.append(0) #some proteins erroneously don't have sequences - skip these
            continue
        
        # iterate over group tuples of (amino acid name, count of sequential repeats) and get the highest value
        result = [(label, sum(1 for _ in group)) for label, group in groups]
        max_AAR = -1
        for j,k in result:
            if k>max_AAR:
                max_AAR = k
        all_max_AARs.append(max_AAR)

    # add the highest values for all proteins back to the dataframe
    protein_df['max_AAR_length'] = all_max_AARs
    return protein_df

def sig_AAR_lengths(protein_df):
    """
    DESCRIPTION: Calculates the mean and standard deviation of the max AAR lengths, and labels all proteins with max AAR lengths 2 
                standard deviations above the mean. 
    INPUTS: 
        protein_df: pandas dataframe of proteins, interactions, sequences, gene names, sequence lengths, peak intereaction sequence 
                    lengths, counts of number of interactions, and longest polypeptide chain of repeated single Amino Acids
    OUTPUTS:
        protein_df: same dataframe as above, but now with all significantly long AAR lengths tagged in a new column (sig_AAR_length)
    """
    # initialize column
    protein_df['sig_AAR_length'] = [False]*len(protein_df)
    
    # calculate mean and standard deviation
    aar_mean = protein_df['max_AAR_length'].mean()
    aar_stdev = protein_df['max_AAR_length'].std()

    # find all proteins with maximum AAR lengths greater than 2 standard deviations above the mean and label those in the original dataframe
    sig_AAR_idxs = protein_df[protein_df['max_AAR_length'] > aar_mean+2*(aar_stdev)].index.values.tolist()
    protein_df.loc[sig_AAR_idxs, 'sig_AAR_length'] = True
    return protein_df

def plot_int_vs_aar(protein_df):
    """
    DESCRIPTION: Creates and saves boxplot of the number of protein interactions for each max AAR length, with significantly long max 
                AAR lengths highlighted separately. 
    INPUTS:
        protein_df: pandas dataframe of proteins, interactions, sequences, gene names, sequence lengths, peak intereaction sequence 
                    lengths, counts of number of interactions, longest polypeptide chain of AARs, and labels of AARs of significant 
                    length
    OUTPUTS:
        None
    """
    sns.boxplot(data=protein_df, y='num_interacts', x='max_AAR_length', hue='sig_AAR_length')

    #label the plot clearly and save
    plt.xlabel('Maximum number of chain Amino Acid Repeats')
    plt.ylabel('Number of interactions')
    plt.title('Number of interactions vs length of Amino Acid Repeats')
    plt.legend(title='AAR length significantly high')
    plt.savefig(f'{BASE_SAVEPATH}/figure2.png')

def label_peak_seq_len_and_aar(protein_df):
    """
    DESCRIPTION: Label proteins in the Insulin protein dataframe that both have sequence lengths that fall within the peak number of 
                    interactions sequence length band of 500-1500 and also have significantly high AAR polypeptide lengths. 
    INPUTS:
        protein_df: pandas dataframe of proteins, interactions, sequences, gene names, sequence lengths, peak intereaction sequence 
                    lengths, counts of number of interactions, longest polypeptide chain of AARs, and labels of AARs of significant 
                    length
    OUTPUTS:
        protein_df: same as input dataframe, but now with sequences that both have sequence lengths that fall within the peak number 
                    of interactions sequence length band of 500-1500 and also have significantly high AAR polypeptide lengths 
                    (peak_len_and_peak_aar)
    """
    both_peaks_idxs = protein_df[(protein_df['peak_int_seq_len']==True)&(protein_df['sig_AAR_length']==True)].index.values.tolist()
    protein_df['peak_len_and_peak_aar'] = [False]*len(protein_df)
    protein_df.loc[both_peaks_idxs, 'peak_len_and_peak_aar'] = True
    return protein_df

def plot_triple_interaction(protein_df):
    """
    DESCRIPTION: Create and save (as interactive html) a 3d interactive plotly scatterplot of the max aar lengths, sequence lengths, 
                and number of interactions for all proteins in the protein dataframe. Proteins which both have sequence lengths that 
                fall within the peak number of interactions sequence length band of 500-1500 and also have significantly high AAR 
                polypeptide lengths are highlighted separately from the remaining sequences to highlight whether the two components 
                of the proteins impact the number of interactions.
    INPUTS:
        protein_df: andas dataframe of proteins, interactions, sequences, gene names, sequence lengths, peak intereaction 
                    sequence lengths, counts of number of interactions, longest polypeptide chain of AARs, labels of AARs of 
                    significant length, and peak sequence and AAR lengths highlighted
    OUTPUTS:
        None
    """
    #create 3d plot of all 3 comparison points
    fig = px.scatter_3d(protein_df, x='max_AAR_length', y='seq_len', z='num_interacts',
                        color='peak_len_and_peak_aar', 
                        labels={"peak_len_and_peak_aar": "Peak sequence length and max AAR length"}, 
                        size_max=5, opacity=0.7)
    
    # organize layout of plot with labels and sizing
    fig.update_layout(
        title='Impact of max AAR length and sequence length on number of protein interactions', 
        autosize=False,
        width=1000, 
        height=1000,
        scene=dict(
            xaxis_title='Maximum AAR length',
            yaxis_title='Sequence length',
            zaxis_title='Number of Interactions',
        ),
    )

    # save interactive plot
    plotly.offline.plot(fig, filename=f'{BASE_SAVEPATH}figure3.html')

if __name__ == "__main__":
    main()
