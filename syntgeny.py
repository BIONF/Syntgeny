# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2023 Hannah Muelbaier
#
#  This script is used to run Syntgeny which performs targeted ortholog
#  searches with fDOG or fDOG-Assembly to investigate if the gene order
#  is conserved over a set of taxa
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: hannah.muelbaier@gmail.com
#
#######################################################################

############################ imports ###########################################
import os
#import os.path
#import sys
from Bio import SeqIO
#from Bio.Phylo.TreeConstruction import DistanceCalculator
#from Bio import AlignIO
import argparse
from collections import deque
#import yaml
import subprocess
import time
#import shutil
#import multiprocessing as mp
import re
from fdog import setupfDog
import pandas as pd


########################### objects ############################################

def starting_subprocess(cmd, mode ='silent', time_out = None):
    try:
        if mode == 'debug':
            result = subprocess.run(cmd, shell=True, timeout = time_out)
        elif mode == 'silent':
            result = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True, timeout = time_out)
        elif mode == 'normal':
            result = subprocess.run(cmd, stdout = subprocess.PIPE, shell=True, timeout = time_out)
    except subprocess.TimeoutExpired:
        return 1
class gene:
    def __init__(self, name, contig, start, end, strand,position, longest_isoform):
        self.name = name
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.position = position
        self.longest_isoform = longest_isoform
        self.protein_ids = set()
        self.mrna_ids = set()
    def set_position(self, position):
        self.position = position
    def set_iso(self, longest_isoform):
        self.longest_isoform = longest_isoform
    def set_protein_ids(self, protein_ids):
        self.protein_ids = protein_ids
    def set_mrna_ids(self, mrns_ids):
        self.mrna_ids = mrns_ids

########################### functions ##########################################

def get_nb_from_gff(gff_path, neighbours, seed):
    gene_stack = deque()
    N = int(neighbours)
    stop = False
    print(seed)
    protein_names = set()
    mrna_names = set()
    #gene = False
    with open(gff_path) as gff:
        for line in gff:
            try:
                contig, source, type, start, end, score, strand, phase, att = line.split('\t')
            except ValueError:
                continue
            if type == 'gene':
                if gene_stack:
                    gene_stack[0].set_protein_ids(protein_names)
                    gene_stack[0].set_mrna_ids(mrna_names)
                    protein_names = set()
                    mrna_names = set()
                name = (att.split(';')[0]).split('ID=')[1]
                if name == seed:
                    N = neighbours*2 +1
                    stop = True
                if len(gene_stack) == N:
                    if stop == False:
                        gene_stack.pop()
                    else:
                        return gene_stack
                gene_stack.appendleft(gene(name, contig, start, end, strand, None, None))
            elif type == 'CDS':
                try:
                    protein_name = re.search(r'Name=(.*?);', att).group(1)
                except AttributeError:
                    protein_name = None
                protein_names.add(protein_name)
            elif type == 'mRNA':
                try:
                    mRNA_name = re.search(r'Name=(.*?);', att).group(1)
                except AttributeError:
                    mRNA_name = None
                mrna_names.add(mRNA_name)

def nb_validation(genes, ng):
    contig_seed = genes[ng].contig
    to_pop= []
    counter = 0
    for obj in genes:
        obj.set_position(ng - counter)
        if obj.contig != contig_seed:
            to_pop.append(obj)
        counter +=1
    if to_pop == []:
        return genes
    else:
        for i in to_pop:
            genes.remove(i)
        return genes
def extract_seq(fasta_path, gene, out_folder):
    isoforms = max(len(gene.mrna_ids), len(gene.protein_ids))
    print(gene.mrna_ids)
    with open(out_folder + gene.name + ".fa", "w") as f:
        max_length = 0
        id = None
        seq = None
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            if isoforms > 0:
                if seq_record.id in gene.mrna_ids:
                    if len(seq_record.seq) > max_length:
                        id = seq_record.id
                        seq = seq_record.seq
                    isoforms -=1
                elif seq_record.id in gene.protein_ids:
                    if len(seq_record.seq) > max_length:
                        id = seq_record.id
                        seq = seq_record.seq
                    isoforms -= 1
            else:
                f.write(">" + str(id) + "\n")
                f.write(str(seq) + "\n")
                return id

def parse_mapping_file():
    pass

def default_data_path():
    sp = setupfDog.get_source_path()
    dp = setupfDog.get_data_path(sp)
    print(dp)
    return dp

def start_fdog(seed_folder, refSpec, searchTaxaDir, annodir, coreTaxaDir, name, out, taxa):
    if taxa == '':
        cmd = f"fdogs.run --seqFolder {seed_folder} --jobName {name} --refspec {refSpec} --outpath {out} --corepath {coreTaxaDir} --searchpath {searchTaxaDir} --annopath {annodir}"
    else:
        cmd = f"fdogs.run --seqFolder {seed_folder} --jobName {name} --refspec {refSpec} --outpath {out} --corepath {coreTaxaDir} --searchpath {searchTaxaDir} --annopath {annodir} --searchTaxa {taxa}"
    print(cmd)
    starting_subprocess(cmd, 'normal')

def parse_profile(fdog_out):
    profile_file = open(fdog_out, 'r')
    profile= profile_file.readlines()

    seed_presence = {}
    seed_orthologs = {}
    genes_to_extract = {}
    for line in profile[1:]:
        seed, ortho_species, ortho_gene, co_ortho = line.split('\t')[2].split('|')
        try:
            seed_presence[seed].add(ortho_species)
        except KeyError:
            seed_presence[seed] = {ortho_species}
            seed_orthologs[seed] = {}
        try:
            seed_orthologs[seed][ortho_species].add(ortho_gene)
        except KeyError:
            seed_orthologs[seed][ortho_species] = {ortho_gene}
        try:
            genes_to_extract[ortho_species].add(ortho_gene)
        except KeyError:
            genes_to_extract[ortho_species] = {ortho_gene}

    return seed_presence, seed_orthologs, genes_to_extract

def get_positions_from_gff(gff_path, genes):
    pass
def main ():
    #################### handle user input #####################################

    start = time.time()
    version = '0.0.1'
    ################### initialize parser ######################################
    parser = argparse.ArgumentParser(description='You are running Syntgeny version ' + str(version) + '.')
    parser.add_argument('--version', action='version', version=str(version))
    ################## required arguments ######################################
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--seed', help='Seed gene name',action='store', default='', required=True)
    required.add_argument('--refSpec', help='Reference taxon', action='store', default='', required=True)
    required.add_argument('--neighbours', help='Number of neighbouring genes up and downstream that should be checked', action ='store', default=3, required=True, type=int)
    required.add_argument('--gff', help='Path to gff file of reference species', action='store', default='', required=True)
    required.add_argument('--jobName', help='Job name', required=True, default='')
    ################## optional arguments ######################################
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--fasta', help='Path to protein fasta file', action='store', default ='')
    optional.add_argument('--searchTaxa', help='File containing search species line by line', action='store', default ='')
    optional.add_argument('--out', help='Output folder', action='store', default='')
    optional.add_argument('--mapping_file', help='Path to file mapping gene ID and protein IDs', action='store', default ='')
    optional.add_argument('--searchpath', help='Path for the search taxa directory', action='store', default='')
    optional.add_argument('--corepath', help='Path for the core taxa directory', action='store', default='')
    optional.add_argument('--annopath', help='Path for the pre-calculated feature annotion directory', action='store', default='')
    args = parser.parse_args()
    # required
    refSpec = args.refSpec
    seed = args.seed
    ng = args.neighbours
    gff_path = args.gff
    jobName = args.jobName
    # not required
    fasta_path = args.fasta
    taxa_path = args.searchTaxa
    out_folder = args.out
    mapping_path = args.mapping_file
    corepath = args.corepath
    annopath= args.annopath
    searchpath= args.searchpath


    #default_path = default_data_path() # in case I want to use the QfO data for core compilation, have ti find a away to add the reference species data

    genes = get_nb_from_gff(gff_path, ng, seed)
    genes = nb_validation(genes,ng)

    for obj in genes:
        print(obj.name, obj.contig, obj.start, obj.end, obj.strand, obj.position, sep=' ')

    for gene in genes:
        seed_folder = out_folder + "/tmp/seeds/"
        cmd = "mkdir -p " + seed_folder
        starting_subprocess(cmd)
        longest_isoform = extract_seq(fasta_path, gene, seed_folder)
        gene.set_iso(longest_isoform)


    #start_fdog(seed_folder, refSpec, searchpath, annopath, corepath, jobName, out_folder, taxa_path)


    fdog_out = f"{out_folder}/{jobName}.phyloprofile"
    seed_presence, seed_orthologs, genes_to_extract = parse_profile(fdog_out)
    print(f'{seed_presence}\n')
    print(f'{seed_orthologs}\n')
    print(f'{genes_to_extract}\n')

    if taxa == '':
        searchTaxa = os.listdir(searchpath)
    else:
        searchTaxa = open(taxa_path,'r').readlines()

    for taxa in searchTaxa:
        gff_file = f"{gff_path}/{taxa}/{taxa}.gff"
        get_positions_from_gff(gff_file, genes_to_extract[taxa])

    end = start - time.time()
    print(end)






main()