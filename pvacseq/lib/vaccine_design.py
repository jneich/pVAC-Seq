#  Input would be a protein fasta
# ie:
# python spacer_junction.py peptides.fa output.csv ann  H-2-Kb -l 8


import sys
import os
from pathlib import Path

root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
import operator
import pandas
import pyprind
import time
import networkx as nx
import matplotlib.pyplot as plt

import itertools
from Bio import SeqIO
import lib
import csv
from prediction_class import *
from parse_output import *

import time

from simanneal import Annealer
import random


def define_parser():
    parser = argparse.ArgumentParser('pvacseq vaccine_design')
    parser.add_argument('input_file', type=argparse.FileType('r'),
                        help="Input FASTA file")
    parser.add_argument(
        "run_name",
        help="The name of the run being processed. This will be used as a prefix for output files"
    )
    parser.add_argument("-k", "--keep-tmp", default=False,
                        help="Option to store tmp files in tmp directory. Default: False")
    parser.add_argument('method',
                        choices=PredictionClass.iedb_prediction_methods(),
                        help="The iedb analysis method to use")
    parser.add_argument('allele',
                        help="Allele for which to make prediction")
    parser.add_argument(
        "-e", "--epitope-length", type=lambda s: [int(epl) for epl in s.split(',')],
        help="Length of subpeptides (neoepitopes) to predict. "
             + "Multiple epitope lengths can be specified using a comma-separated list. "
             + "Typical epitope lengths vary between 8-11. "
             + "Required for Class I prediction algorithms",
    )
    parser.add_argument("-c", "--cutoff", type=int,
                        default=500,
                        help="Optional ic50 cutoff value. Neoepitopes with IC50 values below this " +
                             " value will be labeled as problematic. Default: 500")
    parser.add_argument(
        "-r", "--iedb-retries", type=int,
        default=5,
        help="Number of retries when making requests to the IEDB RESTful web interface. Must be less than or equal to 100."
             + "Default: 5"
    )
    parser.add_argument("-p", "--shortest-paths", type=int, default=20,
                        help="Number of peptide paths to return. Default: 20")
    # parser.add_argument("-s", "--sort-by-lowest", type=bool, default=False,
    #                     help="If True, sorts results by lowest junction IC50 score in peptide rather than " +
    #                          "median score across all junctions. Default: False")
    return parser


class OptimalPeptide(Annealer):

    def __init__(self, state, distance_matrix):
        self.distance_matrix = distance_matrix
        super(OptimalPeptide, self).__init__(state)  # important!

    def move(self):
        """Swaps two peptides in the path."""
        a = random.randint(0, len(self.state) - 1)
        b = random.randint(0, len(self.state) - 1)
        self.state[a], self.state[b] = self.state[b], self.state[a]

    def energy(self):
        """Calculates the length of the route."""
        e = 0
        for i in range(len(self.state)):
            e += self.distance_matrix[self.state[i - 1]][self.state[i]]
        return e



def main(args_input=sys.argv[1:]):

    if "." in args.run_name:
        sys.exit("Run name cannot contain '.'")

    if os.stat(args.input_file).st_size == 0:
        sys.exit("The input fasta is empty.")

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")



    #FIXME add option to keep tmp files. We need to better manage our intermediate files so they are out of ./lib

    #FIXME add option for k shortest paths and if we want to return all shortest paths (warn user that this is O(n!)

    parser = define_parser()
    args = parser.parse_args(args_input)

    input_file = args.input_file
    output_file = args.output_file
    iedb_method = args.method
    ic50_cutoff = args.cutoff
    alleles = args.allele.split(',')
    epl = args.epitope_length
    sort_by_lowest = args.sort_by_lowest
    print("IC50 cutoff: " + str(ic50_cutoff))

    # base_output_dir = os.path.abspath(args.output_dir)
    # tmp_dir = os.path.join(args.output_dir, 'tmp')
    # os.makedirs(tmp_dir, exist_ok=True)
    # self.tmp_dir = tmp_dir

    #   Get fasta info
    peptides = SeqIO.parse(input_file, "fasta")
    seq_dict = dict()
    for record in peptides:
        seq_dict[record.id] = str(record.seq)
    # Get a list of sequence keys
    seq_keys = sorted(seq_dict)

    seq_tuples = list(itertools.combinations_with_replacement(seq_keys, 2))
    combinations = list()
    for key in seq_tuples:
        if (key[0] != key[1]):
            combinations.append((key[0], key[1]))
            combinations.append((key[1], key[0]))
    seq_tuples = combinations

    epitopes = dict()
    rev_lookup = dict()
    for comb in seq_tuples:
        seq1 = comb[0]
        seq2 = comb[1]
        for length in range(8, 11):
            seq_ID = seq1 + "|" + seq2
            trunc_seq1 = seq_dict[seq1][(len(seq_dict[seq1]) - length):len(seq_dict[seq1])]
            trunc_seq2 = seq_dict[seq2][0:(length - 1)]
            epitopes[seq_ID] = trunc_seq1 + trunc_seq2
            rev_lookup[(trunc_seq1 + trunc_seq2)] = seq_ID

            spacers = ["HH", "HHC", "HHH", "HHHD", "HHHC", "AAY", "HHHH", "HHAA", "HHL", "AAL"]
            for this_spacer in spacers:
                seq_ID = seq1 + "|" + this_spacer + "|" + seq2
                epitopes[seq_ID] = (trunc_seq1 + this_spacer + trunc_seq2)
                rev_lookup[(trunc_seq1 + this_spacer + trunc_seq2)] = seq_ID

    with open("epitopes.fa", "w") as tmp:
        for each in epitopes:
            tmp.write(">" + each + "\n" + epitopes[each] + "\n")

    # now we call iedb on the tmp fasta
    output_file = open('iedb_out.csv', 'r+')
    output_file.close()
    outfile = os.path.basename('iedb_out_2.csv')
    file_path = os.path.basename('epitopes.fa')
    split_out = []

    for a in alleles:
        for l in epl:
            print ("Calling iedb for " + a + " of length " + str(l))
            #FIXME need to provide option for local iedb install if available
            lib.call_iedb.main([
                file_path,
                outfile,
                iedb_method,
                a,
                '-l', str(l),
                '-r', str(args.iedb_retries)
            ])
            with open('iedb_out_2.csv', 'rU') as sheet:
                split_out.append(pandas.read_csv(sheet, delimiter='\t'))

    print("IEDB calls complete. Merging data.")
    with open('iedb_out_2.csv', 'rU') as sheet:
        split_out.append(pandas.read_csv(sheet, delimiter='\t'))
    epitope_binding = pandas.concat(split_out)


    # Split into good and bad junctions by filtering the junction table
    problematic_neoepitopes = epitope_binding[epitope_binding.ic50 < ic50_cutoff]
    merged = pandas.DataFrame(pandas.merge(epitope_binding, problematic_neoepitopes, how='outer',
                                           indicator=True).query('_merge == "left_only"').drop(['_merge'], axis=1))
    merged = merged.sort_values('ic50', ascending=False)
    # set the index of the merged values to be the peptide
    peptides = merged.set_index('peptide').T.to_dict('dict')


    keyErrorCount = 0
    successCount = 0
    iedb_results = dict()
    for seqID in epitopes:
        #FIXME eventually get substrings for all lengths of epitopes we tested
        for l in epl:
            #epitopes[seqID][i:i+l]
            for i in range (0, len(epitopes[seqID]) - (l-1)):
                key = epitopes[seqID][i:i+l]
                try:
                    peptides[key]
                except KeyError:
                    keyErrorCount += 1
                    continue

                if seqID not in iedb_results:
                    iedb_results[seqID] = {}
                allele = peptides[key]['allele']
                if allele not in iedb_results[seqID]:
                    iedb_results[seqID][allele] = {}
                    #FIXME not sure we need this, because we want the worst case but with different epl well need it
                    if 'total_score' not in iedb_results[seqID][allele]:
                        iedb_results[seqID][allele]['total_score'] = list()
                        iedb_results[seqID][allele]['total_score'].append(peptides[key]['ic50'])
                    else:
                        iedb_results[seqID][allele]['total_score'].append(peptides[key]['ic50'])

                if 'min_score' in iedb_results[seqID][allele]:
                    iedb_results[seqID][allele]['min_score'] = min(iedb_results[seqID][allele]['min_score'], peptides[key]['ic50'])
                else:
                    iedb_results[seqID][allele]['min_score'] = peptides[key]['ic50']
                    successCount += 1

    print("Successful ic50 mappings: " + str(successCount) + " errors: " + str(keyErrorCount))

    with open('iedb_results_file.txt', 'w') as f:
        f.write(str(iedb_results))


    ## this one takes worst score and maximizes
    Paths = nx.DiGraph()
    spacers = [None, "HH", "HHC", "HHH", "HHHD", "HHHC", "AAY", "HHHH", "HHAA", "HHL", "AAL"]
    for ep in combinations:
        ID_1 = ep[0]
        ID_2 = ep[1]
        Paths.add_node(ID_1)
        Paths.add_node(ID_2)
        best_spacer = ''
        best_allele = None
        for space in spacers:
            if space is None:
                key = str(ID_1 + "|" + ID_2)
            else:
                key = str(ID_1 + "|" + space + "|" + ID_2)
            worst_case = sys.maxsize;
            for allele in iedb_results[key]:
                if iedb_results[key][allele]['min_score'] < worst_case:
                    worst_case = iedb_results[key][allele]['min_score']
            if Paths.has_edge(ID_1, ID_2) and Paths[ID_1][ID_2]['weight'] < worst_case:
                Paths[ID_1][ID_2]['weight'] = worst_case
                if space is not None:
                    Paths[ID_1][ID_2]['spacer'] = space
                else:
                    Paths[ID_1][ID_2]['spacer'] = ''
            elif not Paths.has_edge(ID_1, ID_2):
                if space is not None:
                    Paths.add_edge(ID_1, ID_2, weight=worst_case, spacer=space)
                else:
                    Paths.add_edge(ID_1, ID_2, weight=worst_case, spacer='')


    with open('iedb_results_file.txt', 'r+') as f:
        f.write("\n\n\n\n\n\n")
        f.write(str(Paths))

    print ("Graph contains " + str(len(Paths)) + " nodes and " + str(Paths.size()) + " edges.")
    print("Finding vector paths.")
    print(str(Paths))


    def k_shortest_paths(G, source, target, k, weight=None):
        return list(itertools.islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))

    simple_paths = list()
    for node_s in seq_keys:
        print ("Using node " + node_s + " as starting node. ")
        for node_e in pyprind.prog_bar(seq_keys):
            if node_s is node_e:
                continue
            #print ("Finding path to  " + node_e)
            #paths = nx.shortest_simple_paths(Paths, source=node_s, target=node_e, weight='weight')
            try:
                for path in nx.shortest_simple_paths(Paths, source=node_s, target=node_e, weight='weight'):
                    if (len(path) >= (nx.number_of_nodes(Paths))):
                        # print(path)
                        simple_paths.append(path)
            except nx.NetworkXNoPath:
                print ("No path between" + node_s + " and " + node_e)
                continue

    pos = nx.spring_layout(Paths)  # positions for all nodes

    # nodes
    nx.draw_networkx_nodes(Paths, pos,
                           node_color='b',
                           node_size=500,
                           alpha=0.8, with_labels=True)

    # edges
    nx.draw_networkx_edges(Paths, pos, width=1.0, alpha=0.5)
    nx.draw_networkx_edge_labels(Paths, pos, label_pos=1,font_size=4)
    plt.savefig("simple_path.png", figsize=(7,7))
    #plt.show()



    init_state = seq_keys
    random.shuffle(init_state)

    distance_matrix = {}
    for ID_1 in Paths:
        try:
            distance_matrix[ID_1]
        except KeyError:
            distance_matrix[ID_1] = {}
        for ID_2 in Paths[ID_1]:
            distance_matrix[ID_1][ID_2] = Paths[ID_1][ID_2]['weight']


    peptide = OptimalPeptide(init_state, distance_matrix)

    peptide.copy_strategy = "slice"
    state, e = peptide.anneal()


    while state[0] != seq_keys[0]:     # just use the first key as the end
        state = state[1:] + state[:1]  # rotate key to start
    print("%i distance :" % e)

    for id in state:
        print("\t", id)

    with open('results_temp.txt', 'w') as f:
        output = list()
        name = list()
        min_score = Paths[state[0]][state[1]]['weight']
        cummulative_weight = 0
        all_scores = list()

        for i in range(0, len(state)):
            name.append(state[i])
            try:
                min_score = min(min_score, Paths[state[i]][state[i + 1]]['weight'])
                cummulative_weight += Paths[state[i]][state[i + 1]]['weight']
                all_scores.append(str(Paths[state[i]][state[i + 1]]['weight']))
                spacer = Paths[state[i]][state[i + 1]]['spacer']
                if spacer is not '':
                    name.append(spacer)
            except IndexError:
                #print("Finished traversing path.")
                continue
        median_score = str(cummulative_weight/len(all_scores))
        peptide_id_list = ','.join(name)
        score_list = ','.join(all_scores)
        output = list()
        output.append(">")
        output.append(peptide_id_list)
        output.append("|Median_Junction_Score:")
        output.append(median_score)
        output.append("|Lowest_Junction_Score:")
        output.append(str(min_score))
        output.append("|All_Junction_Scores:")
        output.append(score_list)
        output.append("\n")
        for id in name:
            try:
                output.append(seq_dict[id])
            except KeyError:
                output.append(id)
            output.append("\n")
        # if (sort_by_lowest):
        #     results_dict[min_score] = output
        # else:
        #     results_dict[median_score] = output
        #f.write(''.join(output))
        f.write(''.join(output))


    # Step 7
    #   reformat output into something readable
    # since its almost impossible to get two different peptides with exactly the same MT score, well store the
    # paths as a dictionary of MT scores keys with value 'output'  and as a dictionary of lowest junction score
    # with value 'output'. well sort by that key on user specification

    results_dict = {}
    path_count = 0
    with open('results.txt', 'w') as f:
        #writer = csv.writer(f)
        for path in simple_paths:
            if path_count >= args.shortest_paths:
                break
            name = list()
            cummulative_weight = 0
            #print (path)
            min_score = Paths[path[0]][path[1]]['weight'];
            all_scores = list()
            for i in range(0, len(path)):
                name.append(path[i])
                try:
                    min_score = min(min_score, Paths[path[i]][path[i + 1]]['weight'])
                    cummulative_weight += Paths[path[i]][path[i + 1]]['weight']
                    all_scores.append(str(Paths[path[i]][path[i + 1]]['weight']))
                    spacer = Paths[path[i]][path[i + 1]]['spacer']
                    if spacer is not '':
                        name.append(spacer)
                except IndexError:
                    #print("Finished traversing path.")
                    continue
            median_score = str(cummulative_weight/len(all_scores))
            peptide_id_list = ','.join(name)
            score_list = ','.join(all_scores)
            output = list()
            output.append(">")
            output.append(peptide_id_list)
            output.append("|Median_Junction_Score:")
            output.append(median_score)
            output.append("|Lowest_Junction_Score:")
            output.append(str(min_score))
            output.append("|All_Junction_Scores:")
            output.append(score_list)
            output.append("\n")
            for id in name:
                try:
                    output.append(seq_dict[id])
                except KeyError:
                    output.append(id)
                output.append("\n")
            if (sort_by_lowest):
                results_dict[min_score] = output
            else:
                results_dict[median_score] = output
            f.write(''.join(output))

        # for key, value in sorted(results_dict.items(), reverse=True):
        #     for item in value:
        #         f.write(item)


def removekey(d, key):
    r = dict(d)
    del r[key]
    return r


if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))