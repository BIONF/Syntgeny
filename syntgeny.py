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
import os.path
import sys
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
    def __init__(self, name, contig, start, end, strand,position):
        self.name = name
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.position = position
    def set_position(self, position):
        self.position = position

########################### functions ##########################################

def parse_gff(gff_path, neighbours, seed):
    gene_stack = deque()
    N = int(neighbours)
    stop = False
    print(seed)
    with open(gff_path) as gff:
        for line in gff:
            try:
                contig, source, type, start, end, score, strand, phase, att = line.split('\t')
            except ValueError:
                #print(line)
                continue
            if type == 'gene':
                #print(line)
                name = (att.split(';')[0]).split('ID=')[1]
                if name == seed:
                    N = neighbours*2 +1
                    stop = True

                if len(gene_stack) == N:
                    if stop == False:
                        gene_stack.pop()
                    else:
                        #gene_stack.pop()
                        return gene_stack
                gene_stack.appendleft(gene(name, contig, start, end, strand, None))
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
    with open(out_folder + gene + ".fa", "w") as f:
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            if seq_record.id == gene:
                f.write(str(seq_record.id) + "\n")
                f.write(str(seq_record.seq) + "\n")
                break


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
    required.add_argument('--neighbours', help='Number of neighbouring genes up and downstream that should be checked', action='store', default=3, required=True, type=int)
    required.add_argument('--gff', help='Path to gff file of reference species', action='store', default='', required=True)
    ################## optional arguments ######################################
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--fasta', help='Path to protein fasta file', action='store', default='')
    optional.add_argument('--searchTaxa', help='File containing search species line by line', action='store', default='')
    optional.add_argument('--coreGroups', help='Path to folder that contains precompiled core groups', action='store', default='')
    optional.add_argument('--out', help='Output folder', action='store', default='')
    args = parser.parse_args()
    # required
    refSpec = args.refSpec
    seed = args.seed
    ng = args.neighbours
    gff_path = args.gff
    fasta_path = args.fasta
    taxa_path = args.searchTaxa
    out_folder = args.out

    genes = parse_gff(gff_path, ng, seed)
    genes = nb_validation(genes,ng)

    #for obj in genes:
        #print(obj.name, obj.contig, obj.start, obj.end, obj.strand, obj.position, sep=' ')

    for gene in genes:
        seed_folder = out_folder + "/tmp/seeds/"
        cmd = "mkdir -p " + seed_folder
        starting_subprocess(cmd)
        extract_seq(fasta_path, gene.name, seed_folder)

    end = start - time.time()
    print(end)






main()