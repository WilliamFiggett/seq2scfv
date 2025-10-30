#!/usr/bin/env python3

"""
This script performs iterative IgBLAST analysis on DNA sequences.
It splits sequences based on V, D, and J gene alignments to identify potential V-domains.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd 
import os
import argparse
import logging
import subprocess
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Iterative IgBLAST analysis on DNA sequences to identify V-domains."
    )
    parser.add_argument("--query", required=True, type=str,
                        help="Input file name for DNA sequences in fasta format.")
    parser.add_argument("--germline_V", dest="dbv", required=True, type=str,
                        help="Germline database name for V gene.")
    parser.add_argument("--germline_D", dest="dbd", required=True, type=str,
                        help="Germline database name for D gene.")
    parser.add_argument("--germline_J", dest="dbj", required=True, type=str,
                        help="Germline database name for J gene.")
    parser.add_argument("--organism", required=True, type=str,
                        help="Organism of the query sequence: 'human', 'mouse', 'rat', 'rabbit' or 'custom' (see IgBLAST web site).")
    parser.add_argument("--ig_seqtype", dest="igtype", default="Ig", type=str,
                        help="Receptor sequence type: 'Ig' or 'TCR'.")
    parser.add_argument("--domain_system", dest="domainsys", default="imgt", type=str,
                        help="Domain system to be used for segment annotation: 'imgt' or 'kabat'.")
    parser.add_argument("--auxiliary_data", dest="auxdata", required=True, type=str,
                        help="File containing the coding frame start positions for sequences in germline J database.")
    parser.add_argument("--num_threads", dest="nthreads", default=1, type=int,
                        help="Number of threads (CPUs) to use in the BLAST search. Default: 1.")
    parser.add_argument("--minlen", default=150, type=int,
                        help="Minimum length (bp) required to continue searching for a V-domain. Default: 150 bp.")
    parser.add_argument("--escore", default=1.0e-02, type=float,
                        help="Minimum V and J alignment support (E, expectation value) to consider a hit fully characterizing a V-(D)-J.")
    parser.add_argument("--log", dest="loglevel", default="info", type=str,
                        choices=["debug", "info", "warning", "error", "critical"],
                        help="Logging level (debug, info, warning, error, critical). Default: info.")
    return parser.parse_args()

relative_pos = [
    "v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end", 
    "j_sequence_start", "j_sequence_end", "cdr1_start", "cdr1_end", "cdr2_start", 
    "cdr2_end", "cdr3_start", "cdr3_end", "fwr1_start", "fwr1_end", "fwr2_start", 
    "fwr2_end", "fwr3_start", "fwr3_end", "fwr4_start", "fwr4_end"
]

def check_file_exists(filepath, description):
    if not os.path.isfile(filepath):
        logging.error(f"{description} file '{filepath}' does not exist.")
        sys.exit(1)

def run_igblast(germlineV, germlineD, germlineJ, organism, igseqtype, threads, domain, auxiliary, queryfile, outfile):
    """
    Run IgBLAST to align sequences to germline genes.
    """
    threads = str(threads)
    command = [
        'igblastn', '-germline_db_V', germlineV, '-germline_db_D', germlineD, 
        '-germline_db_J', germlineJ, '-organism', organism, '-ig_seqtype', igseqtype, 
        '-num_threads', threads, '-domain_system', domain, '-num_alignments_V', '1', 
        '-num_alignments_D', '1', '-num_alignments_J', '1', '-outfmt', '19', 
        '-extend_align5end', '-extend_align3end', '-auxiliary_data', auxiliary, 
        '-query', queryfile, '-out', outfile
    ]
    try:
        logging.info(f"Running IgBLAST: {' '.join(command)}")
        subprocess.check_output(command, universal_newlines=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        logging.error(f"IgBLAST failed with exit code {e.returncode}: {e.output}")
        sys.exit(1)
    except FileNotFoundError:
        logging.error("igblastn executable not found in PATH. Please ensure IgBLAST is installed and accessible.")
        sys.exit(1)

def igtable(input):
    """
    Read IgBLAST output and add coordinates to the sequence ID.

    Args:
        input (str): Path to IgBLAST output file.

    Returns:
        pd.DataFrame: Formatted IgBLAST output.
    """
    try:
        blast = pd.read_csv(input, sep='\t', header=0)
    except FileNotFoundError:
        logging.error(f"IgBLAST output file '{input}' not found.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        logging.error(f"IgBLAST output file '{input}' is empty or corrupt.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading IgBLAST output file '{input}': {e}")
        sys.exit(1)

    blast['sequence'] = blast['sequence'].astype(str)
    if not blast["sequence_id"].str.contains(r'\|').any():
        blast["sequence_id"] = blast.apply(lambda row: row["sequence_id"] + "|0:" + str(len(row["sequence"])), axis=1)
        blast["origin"] = "uncut"
    return blast

def blastfilt(blast, escore, output2):
    """
    Filter IgBLAST output based on E score to select well-delimited V-domains.

    Args:
        blast (pd.DataFrame): IgBLAST output.
        escore (float): Minimum E score.
        output2 (str): Path to save filtered output.

    Returns:
        pd.DataFrame: Filtered IgBLAST output.
    """
    try:
        filtered_blast = blast[(blast["v_support"] < float(escore)) & (blast["j_support"] < float(escore))].copy()
        filtered_blast["short"] = filtered_blast["sequence_id"].str.split('|').str[0]
        filtered_blast["coord"] = filtered_blast["sequence_id"].str.split('|').str[1].str.split(':').str[0]
        filtered_blast.reset_index(drop=True, inplace=True)
        filtered_blast.to_csv(output2, sep="\t")
    except KeyError as e:
        logging.error(f"Missing expected column in IgBLAST output: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error filtering IgBLAST output: {e}")
        sys.exit(1)
    return(filtered_blast) 

def scan5p(filtered_blast, minlen): 
    """
    Scan the 5' end for sufficient sequence length. Returns sequences meeting length requirement on 5' end.
    """
    try:
        seq5p = filtered_blast[filtered_blast["v_sequence_start"] > int(minlen)].copy()
        if not seq5p.empty: 
            seq5p["sequence_id"] = seq5p.apply(lambda row: row["short"] + "|" + str(row["coord"]) + ":" + str(round(row["v_sequence_start"])), axis=1)  
            seq5p["sequence"] = seq5p.apply(lambda row: row["sequence"][:int(row["v_sequence_start"])], axis=1)
        else: 
            seq5p = pd.DataFrame()
    except KeyError as e:
        logging.error(f"Missing expected column in filtered IgBLAST output: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error scanning 5' end: {e}")
        sys.exit(1)
    return seq5p

def scan3p(filtered_blast, minlen): 
    """
    Scan the 3' end for sufficient sequence length. Returns sequences meeting length requirement on 3' end.
    """
    try:
        seq3p = filtered_blast[(filtered_blast["sequence"].str.len() - filtered_blast["j_sequence_end"]) > int(minlen)].copy()
        if not seq3p.empty: 
            seq3p["sequence_id"] = seq3p.apply(lambda row: row["short"] + "|" + str(round(int(row["j_sequence_end"]) + int(row["coord"]))) + ":" + str(len(row["sequence"])+ int(row["coord"])), axis=1)
            seq3p["sequence"] = seq3p.apply(lambda row: row["sequence"][(int(row["j_sequence_end"])-1):], axis=1)
        else: 
            seq3p = pd.DataFrame()
    except KeyError as e:
        logging.error(f"Missing expected column in filtered IgBLAST output: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error scanning 3' end: {e}")
        sys.exit(1)
    return seq3p 

def vscan(filtered_blast, minlen):  
    """
    Combine 5' and 3' end scans.

    Args:
        filtered_blast (pd.DataFrame): Filtered IgBLAST output.
        minlen (int): Minimum length required.

    Returns:
        pd.DataFrame: Combined sequences from 5' and 3' scans.
    """
    seq5p = scan5p(filtered_blast, minlen)
    seq3p = scan3p(filtered_blast, minlen)
    if seq5p.empty and seq3p.empty: 
        newset = pd.DataFrame()
    elif seq5p.empty: 
        newset = seq3p[["sequence_id", "sequence"]]
    elif seq3p.empty: 
        newset = seq5p[["sequence_id", "sequence"]]
    else:
        newset = pd.concat([seq5p[["sequence_id", "sequence"]], seq3p[["sequence_id", "sequence"]]], axis=0)
    return newset

def df_to_fasta(newset):
    """
    Convert DataFrame to FASTA format.

    Args:
        newset (pd.DataFrame): DataFrame of sequences.

    Returns:
        list: List of SeqRecord objects.
    """
    try:
        if not newset.empty: 
            newset.reset_index(drop=True, inplace=True)
            newsetlist = list(newset.itertuples(index=False, name=None))
            newfasta = []
            for i in newsetlist:
                newfasta.append(SeqRecord(Seq(i[1]), id=i[0], name="", description=""))
        else: 
            newfasta = []
    except Exception as e:
        logging.error(f"Error converting DataFrame to FASTA format: {e}")
        sys.exit(1)
    return newfasta

def splitfasta(input, minlen, escore, output1, output2):
    """
    Process IgBLAST output and prepare new fragments for search.

    Args:
        input (str): Path to IgBLAST output file.
        minlen (int): Minimum length required.
        escore (float): Minimum E score.
        output1 (str): Path to save new fragments in FASTA format.
        output2 (str): Path to save filtered output.

    Returns:
        None
    """   
    df = igtable(input)
    filt_df = blastfilt(df, escore, output2)
    vscan_df = vscan(filt_df, minlen)
    try:
        SeqIO.write(df_to_fasta(vscan_df), output1, "fasta")
    except Exception as e:
        logging.error(f"Error writing FASTA file '{output1}': {e}")
        sys.exit(1)

def rescanner(germlineV, germlineD, germlineJ, organism, igseqtype, threads, domain, auxiliary, queryfile, minlen, escore): 
    """
    Perform iterative IgBLAST analysis and split sequences. Takes as input IgBLAST args and returns list of paths to filtered output files. 
    """    
    counter = 0 
    res = []
    while True:
        try:
            if not os.path.isfile(queryfile):
                logging.info(f"Query file '{queryfile}' does not exist. Stopping iteration.")
                break
            if os.stat(queryfile).st_size == 0:
                logging.info(f"Query file '{queryfile}' is empty. Stopping iteration.")
                break
        except Exception as e:
            logging.error(f"Error checking query file '{queryfile}': {e}")
            sys.exit(1)
        counter +=1
        igblasttsv = f"ddDNA_igblast_{counter}.tsv"
        run_igblast(germlineV, germlineD, germlineJ, organism, igseqtype, threads, domain, auxiliary, queryfile, igblasttsv)
        splitfastatsv = f"ddDNA_splitfasta_{counter}.tsv"
        queryfile = f"ddDNA_splitfasta_{counter}.fasta"
        splitfasta(igblasttsv, minlen, escore, queryfile, splitfastatsv) 
        res.append(splitfastatsv)
    return res

def gatherdf(result, changepos): 
    """
    Combine all filtered results and adjust coordinates.

    Args:
        result (list): List of paths to filtered output files.
        changepos (list): List of columns to adjust coordinates.

    Returns:
        pd.DataFrame: Combined and adjusted results.
    """    
    dfs = []
    for i in result: 
        try:
            df = pd.read_csv(i, sep='\t', header=0, index_col=0)
        except FileNotFoundError:
            logging.error(f"Filtered result file '{i}' not found.")
            sys.exit(1)
        except Exception as e:
            logging.error(f"Error loading filtered result file '{i}': {e}")
            sys.exit(1)
        dfs.append(df)
    try:
        result_df = pd.concat(dfs, ignore_index=True)
        for j in changepos: 
            result_df.loc[result_df[j].notna(), j] += result_df.loc[result_df[j].notna(), "coord"]
        # sequence field represents the full length sequence
        fullseq = result_df.groupby("short").filter(lambda x: any(x['origin'] == 'uncut'))
        uncut_mapping = (
            fullseq.loc[fullseq['origin'] == 'uncut']
            .set_index('short')
            ['sequence']
            .to_dict()
        )
        result_df.loc[(result_df['origin'] != 'uncut') & (result_df['short'].isin(uncut_mapping)), 'sequence'] = (
            result_df['short'].map(uncut_mapping)
        )
        # remove coordinate data (original position now taken into account in table)
        result_df = result_df.drop(columns=["sequence_id", "coord", "origin"])
        newcol = result_df.pop("short")
        result_df.insert(0, "sequence_id", newcol)
        result_df.to_csv("igBLAST.tsv", sep="\t")
    except Exception as e:
        logging.error(f"Error during DataFrame merge and coordinate adjustment: {e}")
        sys.exit(1)
    return result_df

def main(args, relative_pos): 
    # Check required files before starting
    check_file_exists(args.query, "Query")
    check_file_exists(args.dbv, "V germline database")
    check_file_exists(args.dbd, "D germline database")
    check_file_exists(args.dbj, "J germline database")
    check_file_exists(args.auxdata, "Auxiliary data")
    res = rescanner(args.dbv, args.dbd, args.dbj, args.organism, args.igtype, args.nthreads, args.domainsys, args.auxdata, args.query, args.minlen, args.escore)
    finaldf = gatherdf(res, relative_pos)    

if __name__ == '__main__':
    args = parse_args()
    numeric_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_level, int):    
        raise ValueError('Invalid log level: %s' % args.loglevel)
    logging.basicConfig(level=numeric_level)
    logging.debug(args)
    main(args, relative_pos)
