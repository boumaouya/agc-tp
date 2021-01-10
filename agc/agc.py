#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
#from operator import itemgetter

# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Maouya Bou"
__copyright__ = "EISTI / CY Tech"
__credits__ = ["Maouya Bou"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Maouya Bou"
__email__ = "boumaouya@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """ Prend deux arguments correspondant au fichier fasta ou fasta.gz et à la longueur minimale des
    séquences et retourne un ​ générateur de séquences de longueur l >= ​ minseqlen ​ : ​ yield
    sequence.
    """
    if (isfile(amplicon_file)):
        with gzip.open(amplicon_file, "rt") as fpt:
            contig =''
            for apt in fpt:
                if len(apt.strip())==0 or apt[0]=='>' :
                    if len(contig)>=1:
                        if len(contig) >= minseqlen:
                            yield contig
                        contig =''
                else:
                    contig = contig + apt.strip()

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """ Prend trois arguments correspondant au fichier fasta, la longueur minimale des séquences et
    leur comptage minimum. Retourne un générateur des séquences uniques ayant une occurrence O>=​ mincount ainsi que leur
    occurrence
    """
    sequences = list(read_fasta(amplicon_file, minseqlen))
    liste = Counter(sequences).most_common()
    for sequence in liste:
        if sequence[1]>=mincount:
            yield sequence

def get_chunks(sequence, chunk_size):
    """ Prend une séquence et une longueur de segment l: ​ chunk_size et retourne une liste de
    sous-séquences de taille l non chevauchantes
    """
    liste_segment = []
    cpt = 0
    for i in range(4):
        liste_segment.append(sequence[cpt:cpt+chunk_size])
        cpt = cpt+chunk_size

    return liste_segment

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """Prend une séquence et une longueur de séquence k et retourne un générateur de tous les mots
    de longueur k présents dans cette séquence
    """
    for i, chaine in enumerate(sequence):
        if(len(sequence)-i >= kmer_size):
            seq = sequence[i:i+kmer_size]
            yield seq


def get_identity(alignment_list):
    """Prend un alignement (sous forme de liste) et calcule le pourcentage d’identité entre les deux
    nb nucléotides identiques séquences selon la formule
    """
    count=0
    for kpt in range(len(alignment_list[0])):
        if alignment_list[0][kpt] == alignment_list[1][kpt]: count=count+1
    id_= count *100 / len(alignment_list[0])
    return id_

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Prend un dictionnaire ayant pour clé un dictionnaire de kmer (vide au départ)
    et pour valeur une liste d’identifiant des séquences dont ils proviennent.
    Retourne un dictionnaire de kmer contenant les kmers uniques présents dans chaque séquence
    pour une longueur de donnée de kmer
    """

    kmer_list = list(cut_kmer(sequence, kmer_size))

    for kmer in kmer_list:
        if kmer in kmer_dict:
            kmer_dict["kmer"].append(id_seq)
        else:
            kmer_dict["kmer"] = [id_seq]

    return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
    """Prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste d’identifiant
    des séquences dont ils proviennent, une séquence et une longueur de kmer.
    Retourne une liste de 8 séquences les plus similaires à la séquence d'entrée"""

    liste_ids = []
    liste_seq=[]
    kmers=list(cut_kmer(sequence, kmer_size))

    for kmer in kmers:
        if kmer in kmer_dict:
            for ipt in kmer_dict[kmer]:
                liste_ids.append(ipt)

    for apt in Counter(liste_ids).most_common(8):
        liste_seq.append(apt[0])
    return liste_seq

def detect_chimera(perc_identity_matrix):
    """Prend une matrice donnant par segment le taux d’identité entre la séquence candidate et deux
    séquences parentes et retourne un booléen indiquant si la séquence candidate est une chimère.
    Si l’écart type moyen des pourcentages est supérieur à ​ 5 et que 2 segments minimum
    de notre séquence montrent une similarité différente à un des deux parents, nous
    identifierons cette séquence comme chimérique.
    """
    ecarts_types = []
    sequence1 = False
    sequence2 = False

    for apt in perc_identity_matrix:
        ecarts_types.append(statistics.stdev(apt))
        if apt[0] > apt[1]:
            sequence1 = True
        if apt[0] < apt[1]:
            sequence2 = True
        if statistics.mean(ecarts_types) > 5 and sequence1 and sequence2:
            return True
        return False

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Fait appel au générateur fourni par ​ dereplication_fulllength et retourne un générateur des
    séquences non chimérique au format: yield [sequence, count]
    """
    #### Cette fonction n'est pas compléte ##########
    kmer_dict = {}
    sequences = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))

    for sequence in sequences:
        chunk_list = get_chunks(sequence[0], chunk_size)
        mates_list = []
        for chunk in chunk_list:
            mates_list.append(search_mates(kmer_dict, chunk, kmer_size))

        yield sequence

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Fait appel à chimera removal et retourne une liste d’OTU, cette liste indiquera pour chaque
    séquence son occurrence (count).
    """
    liste=[]
    for apt in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
        liste.append(apt)
    return liste

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """Prend une liste d’OTU et le chemin vers un fichier de sortie et affiche les OTU au format:
    >OTU_{numéro partant de 1} occurrence:{nombre d’occurrence à la déréplication}
    {séquence au format fasta}
    """
    # Ecrire dans le fichier
    with open(output_file, "wt") as fi:
        for k,contigs in enumerate(OTU_list):
            fi.write(f">OTU_{k+1} occurrence:{contigs[1]}\n" + fill(contigs[0], width=80) + "\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = args.chunk_size
    kmer_size = args.kmer_size
    output_file = args.output_file

    OTU_list = abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    write_OTU(OTU_list, output_file)

if __name__ == '__main__':
    main()