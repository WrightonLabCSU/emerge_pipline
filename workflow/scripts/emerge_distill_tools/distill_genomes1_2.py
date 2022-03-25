"""Make a distillate for DRAM"""
from os import path, mkdir
from itertools import tee
import re
import pandas as pd
from collections import Counter, defaultdict
import networkx as nx
import numpy as np
from datetime import datetime
import click

from utils1_2 import get_ids_from_annotation, get_ids_from_row, get_ordered_uniques

# TODO: add RBH information to output
# TODO: add flag to output table and not xlsx
# TODO: add flag to output heatmap table

FRAME_COLUMNS = ['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader']
RRNA_TYPES = ['5S rRNA', '16S rRNA', '23S rRNA']
HEATMAP_MODULES = ['M00001', 'M00004', 'M00008', 'M00009', 'M00012', 'M00165', 'M00173', 'M00374', 'M00375',
                   'M00376', 'M00377', 'M00422', 'M00567']
HEATMAP_CELL_HEIGHT = 10
HEATMAP_CELL_WIDTH = 10
KO_REGEX = r'^K\d\d\d\d\d$'
ETC_COVERAGE_COLUMNS = ['module_id', 'module_name', 'complex', 'genome', 'path_length', 'path_length_coverage',
                        'percent_coverage', 'genes', 'missing_genes', 'complex_module_name']
TAXONOMY_LEVELS = ['d', 'p', 'c', 'o', 'f', 'g', 's']

DATA_FOLDER = path.join(path.dirname(path.realpath(__file__)), 'data')

def fill_genome_summary_frame(annotations, genome_summary_frame, groupby_column):
    genome_summary_id_sets = [set([k.strip() for k in j.split(',')]) for j in genome_summary_frame['gene_id']]
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        id_dict = get_ids_from_annotation(frame)
        counts = list()
        for i in genome_summary_id_sets:
            identifier_count = 0
            for j in i:
                if j in id_dict:
                    identifier_count += id_dict[j]
            counts.append(identifier_count)
        genome_summary_frame[genome] = counts
    return genome_summary_frame


def fill_genome_summary_frame_gene_names(annotations, genome_summary_frame, groupby_column):
    genome_summary_id_sets = [set([k.strip() for k in j.split(',')]) for j in genome_summary_frame['gene_id']]
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        # make dict of identifiers to gene names
        id_gene_dict = defaultdict(list)
        for gene, row in frame.iterrows():
            ids = get_ids_from_row(row)
            for id_ in ids:
                id_gene_dict[id_].append(gene)
        # fill in genome summary_frame
        values = list()
        for id_set in genome_summary_id_sets:
            this_value = list()
            for id_ in id_set:
                this_value += id_gene_dict[id_]
            values.append(','.join(this_value))
        genome_summary_frame[genome] = values
    return genome_summary_frame


def summarize_rrnas(rrnas_df, groupby_column='fasta'):
    genome_rrna_dict = dict()
    for genome, frame in rrnas_df.groupby(groupby_column):
        genome_rrna_dict[genome] = Counter(frame['type'])
    row_list = list()
    for rna_type in RRNA_TYPES:
        row = [rna_type, '%s ribosomal RNA gene' % rna_type.split()[0], 'rRNA', 'rRNA', '', '']
        for genome, rrna_dict in genome_rrna_dict.items():
            row.append(genome_rrna_dict[genome].get(rna_type, 0))
        row_list.append(row)
    rrna_frame = pd.DataFrame(row_list, columns=FRAME_COLUMNS + list(genome_rrna_dict.keys()))
    return rrna_frame


def summarize_trnas(trnas_df, groupby_column='fasta'):
    # first build the frame
    combos = {(line.Type, line.Codon, line.Note) for _, line in trnas_df.iterrows()}
    frame_rows = list()
    for combo in combos:
        if combo[2] == 'pseudo':
            gene_id = '%s, pseudo (%s)'
            gene_description = '%s pseudo tRNA with %s Codon'
        else:
            gene_id = '%s (%s)'
            gene_description = '%s tRNA with %s Codon'
        gene_id = gene_id % (combo[0], combo[1])
        gene_description = gene_description % (combo[0], combo[1])
        module_description = '%s tRNA' % combo[0]
        frame_rows.append([gene_id, gene_description, module_description, 'tRNA', 'tRNA', ''])
    trna_frame = pd.DataFrame(frame_rows, columns=FRAME_COLUMNS)
    trna_frame = trna_frame.sort_values('gene_id')
    # then fill it in
    trna_frame = trna_frame.set_index('gene_id')
    for group, frame in trnas_df.groupby(groupby_column):
        gene_ids = list()
        for index, line in frame.iterrows():
            if line.Note == 'pseudo':
                gene_id = '%s, pseudo (%s)'
            else:
                gene_id = '%s (%s)'
            gene_ids.append(gene_id % (line.Type, line.Codon))
        trna_frame[group] = pd.Series(Counter(gene_ids))
    trna_frame = trna_frame.reset_index()
    trna_frame = trna_frame.fillna(0)
    return trna_frame


def make_genome_summary(annotations, genome_summary_frame, trna_frame=None, rrna_frame=None, groupby_column='fasta'):
    summary_frames = list()
    # get ko summaries
    summary_frames.append(fill_genome_summary_frame(annotations, genome_summary_frame.copy(), groupby_column))

    # add rRNAs
    if rrna_frame is not None:
        summary_frames.append(summarize_rrnas(rrna_frame, groupby_column))

    # add tRNAs
    if trna_frame is not None:
        summary_frames.append(summarize_trnas(trna_frame, groupby_column))

    # merge summary frames
    summarized_genomes = pd.concat(summary_frames, sort=False)
    return summarized_genomes


def write_summarized_genomes_to_xlsx(summarized_genomes, output_file):
    # turn all this into an xlsx
    with pd.ExcelWriter(output_file) as writer:
        for sheet, frame in summarized_genomes.groupby('sheet', sort=False):
            frame = frame.sort_values(['header', 'subheader', 'module', 'gene_id'])
            frame = frame.drop(['sheet'], axis=1)
            frame.to_excel(writer, sheet_name=sheet, index=False)


# TODO: add assembly stats like N50, longest contig, total assembled length etc
def make_genome_stats(annotations, rrna_frame=None, trna_frame=None, groupby_column='fasta'):
    rows = list()
    columns = ['genome']
    if 'scaffold' in annotations.columns:
        columns.append('number of scaffolds')
    if 'bin_taxonomy' in annotations.columns:
        columns.append('taxonomy')
    if 'bin_completeness' in annotations.columns:
        columns.append('completeness score')
    if 'bin_contamination' in annotations.columns:
        columns.append('contamination score')
    if rrna_frame is not None:
        columns += RRNA_TYPES
    if trna_frame is not None:
        columns.append('tRNA count')
    if 'bin_completeness' in annotations.columns and 'bin_contamination' in annotations.columns \
       and rrna_frame is not None and trna_frame is not None:
        columns.append('assembly quality')
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        row = [genome]
        if 'scaffold' in frame.columns:
            row.append(len(set(frame['scaffold'])))
        if 'bin_taxonomy' in frame.columns:
            row.append(frame['bin_taxonomy'][0])
        if 'bin_completeness' in frame.columns:
            row.append(frame['bin_completeness'][0])
        if 'bin_contamination' in frame.columns:
            row.append(frame['bin_contamination'][0])
        has_rrna = list()
        if rrna_frame is not None:
            genome_rrnas = rrna_frame.loc[rrna_frame.fasta == genome]
            for rrna in RRNA_TYPES:
                sixteens = genome_rrnas.loc[genome_rrnas.type == rrna]
                if sixteens.shape[0] == 0:
                    row.append('')
                    has_rrna.append(False)
                elif sixteens.shape[0] == 1:
                    row.append('%s (%s, %s)' % (sixteens['scaffold'].iloc[0], sixteens.begin.iloc[0],
                                                sixteens.end.iloc[0]))
                    has_rrna.append(True)
                else:
                    row.append('%s present' % sixteens.shape[0])
                    has_rrna.append(False)
        if trna_frame is not None:
            row.append(trna_frame.loc[trna_frame[groupby_column] == genome].shape[0])  # TODO: remove psuedo from count?
        if 'assembly quality' in columns:
            if frame['bin_completeness'][0] > 90 and frame['bin_contamination'][0] < 5 and np.all(has_rrna) and \
               len(set(trna_frame.loc[trna_frame[groupby_column] == genome].Type)) >= 18:
                row.append('high')
            elif frame['bin_completeness'][0] >= 50 and frame['bin_contamination'][0] < 10:
                row.append('med')
            else:
                row.append('low')
        rows.append(row)
    genome_stats = pd.DataFrame(rows, columns=columns)
    return genome_stats


def build_module_net(module_df):
    """Starts with a data from including a single module"""
    # build net from a set of module paths
    num_steps = max([int(i.split(',')[0]) for i in set(module_df['path'])])
    module_net = nx.DiGraph(num_steps=num_steps, module_id=list(module_df['module'])[0],
                            module_name=list(module_df['module_name'])[0])
    # go through all path/step combinations
    for module_path, frame in module_df.groupby('path'):
        split_path = [int(i) for i in module_path.split(',')]
        step = split_path[0]
        module_net.add_node(module_path, kos=set(frame['ko']))
        # add incoming edge
        if step != 0:
            module_net.add_edge('end_step_%s' % (step - 1), module_path)
        # add outgoing edge
        module_net.add_edge(module_path, 'end_step_%s' % step)
    return module_net


def get_module_step_coverage(kos, module_net):
    # prune network based on what kos were observed
    pruned_module_net = module_net.copy()
    module_kos_present = set()
    for node, data in module_net.nodes.items():
        if 'kos' in data:
            ko_overlap = data['kos'] & kos
            if len(ko_overlap) == 0:
                pruned_module_net.remove_node(node)
            else:
                module_kos_present = module_kos_present | ko_overlap
    # count number of missing steps
    missing_steps = 0
    for node, data in pruned_module_net.nodes.items():
        if ('end_step' in node) and (pruned_module_net.in_degree(node) == 0):
            missing_steps += 1
    # get statistics
    num_steps = pruned_module_net.graph['num_steps'] + 1
    num_steps_present = num_steps - missing_steps
    coverage = num_steps_present / num_steps
    return num_steps, num_steps_present, coverage, sorted(module_kos_present)


def make_module_coverage_df(annotation_df, module_nets):
    kos_to_genes = defaultdict(list)
    ko_col = 'ko_id' if 'ko_id' in annotation_df.columns else 'kegg_id'
    for gene_id, ko_list in annotation_df[ko_col].iteritems():
        if type(ko_list) is str:
            for ko in ko_list.split(','):
                kos_to_genes[ko].append(gene_id)
    coverage_dict = {}
    for i, (module, net) in enumerate(module_nets.items()):
        module_steps, module_steps_present, module_coverage, module_kos = \
            get_module_step_coverage(set(kos_to_genes.keys()), net)
        module_genes = sorted([gene for ko in module_kos for gene in kos_to_genes[ko]])
        coverage_dict[module] = [net.graph['module_name'], module_steps, module_steps_present, module_coverage,
                                 len(module_kos), ','.join(module_kos), ','.join(module_genes)]
    coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index',
                                         columns=['module_name', 'steps', 'steps_present', 'step_coverage', 'ko_count',
                                                  'kos_present', 'genes_present'])
    return coverage_df


def make_module_coverage_frame(annotations, module_nets, groupby_column='fasta'):
    # go through each scaffold to check for modules
    module_coverage_dict = dict()
    for group, frame in annotations.groupby(groupby_column, sort=False):
        module_coverage_dict[group] = make_module_coverage_df(frame, module_nets)
    module_coverage = pd.concat(module_coverage_dict)
    module_coverage.index = module_coverage.index.set_names(['genome', 'module'])
    return module_coverage.reset_index()



def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def first_open_paren_is_all(str_):
    """Go through string and return true"""
    curr_level = 1
    for char in str_[1:-1]:
        if char == ')':
            curr_level -= 1
        elif char == '(':
            curr_level += 1
        if curr_level == 0:
            return False
    return True


def split_into_steps(definition, split_char=' '):
    """Very fancy split on string of chars"""
    curr_level = 0
    step_starts = [-1]
    for i, char in enumerate(definition):
        if char == '(':
            curr_level += 1
        if char == ')':
            curr_level -= 1
        if (curr_level == 0) and (char in split_char):
            step_starts.append(i)
    step_starts.append(len(definition))
    steps = list()
    for a, b in pairwise(step_starts):
        step = definition[a+1:b]
        if step.startswith('(') and step.endswith(')'):
            if first_open_paren_is_all(step):
                step = step[1:-1]
        steps.append(step)
    return steps


def is_ko(ko):
    return re.match(KO_REGEX, ko) is not None


def make_module_network(definition, network: nx.DiGraph = None, parent_nodes=('start',)):
    # TODO: Figure out how to add 'end' step to last step at end
    if network is None:
        network = nx.DiGraph()
    last_steps = []
    for step in split_into_steps(definition, ','):
        prev_steps = parent_nodes
        for substep in split_into_steps(step, '+'):
            if is_ko(substep):
                for prev_step in prev_steps:
                    network.add_edge(prev_step, substep)
                prev_steps = [substep]
            else:
                network, prev_steps = make_module_network(substep, network, prev_steps)
        last_steps += prev_steps
    return network, last_steps


def get_module_coverage(module_net: nx.DiGraph, genes_present: set):
    max_coverage = -1
    max_coverage_genes = list()
    max_coverage_missing_genes = list()
    max_path_len = 0
    for net_path in nx.all_simple_paths(module_net, source='start', target='end'):
        net_path = set(net_path[1:-1])
        overlap = net_path & genes_present
        coverage = len(overlap) / len(net_path)
        if coverage > max_coverage:
            max_coverage = coverage
            max_coverage_genes = overlap
            max_coverage_missing_genes = net_path - genes_present
            max_path_len = len(net_path)
    return max_path_len, len(max_coverage_genes), max_coverage, max_coverage_genes, max_coverage_missing_genes


def make_etc_coverage_df(etc_module_database, annotations, groupby_column='fasta'):
    etc_coverage_df_rows = list()
    for _, module_row in etc_module_database.iterrows():
        definition = module_row['definition']
        # remove optional subunits
        definition = re.sub(r'-K\d\d\d\d\d', '', definition)
        module_net, _ = make_module_network(definition)
        # add end node
        no_out = [node for node in module_net.nodes() if module_net.out_degree(node) == 0]
        for node in no_out:
            module_net.add_edge(node, 'end')
        # go through each genome and check pathway coverage
        for group, frame in annotations.groupby(groupby_column):
            # get annotation genes
            grouped_ids = set(get_ids_from_annotation(frame).keys())
            path_len, path_coverage_count, path_coverage_percent, genes, missing_genes = \
                get_module_coverage(module_net, grouped_ids)
            complex_module_name = 'Complex %s: %s' % (module_row['complex'].replace('Complex ', ''),
                                                      module_row['module_name'])
            etc_coverage_df_rows.append([module_row['module_id'], module_row['module_name'],
                                         module_row['complex'].replace('Complex ', ''), group, path_len,
                                         path_coverage_count, path_coverage_percent, ','.join(sorted(genes)),
                                         ','.join(sorted(missing_genes)), complex_module_name])
    return pd.DataFrame(etc_coverage_df_rows, columns=ETC_COVERAGE_COLUMNS)



def make_functional_df(annotations, function_heatmap_form, groupby_column='fasta'):
    # clean up function heatmap form
    function_heatmap_form = function_heatmap_form.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    function_heatmap_form = function_heatmap_form.fillna('')
    # build dict of ids per genome
    genome_to_id_dict = dict()
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        id_list = get_ids_from_annotation(frame).keys()
        genome_to_id_dict[genome] = set(id_list)
    # build long from data frame
    rows = list()
    for function, frame in function_heatmap_form.groupby('function_name', sort=False):
        for bin_name, id_set in genome_to_id_dict.items():
            presents_in_bin = list()
            functions_present = set()
            for _, row in frame.iterrows():
                function_id_set = set([i.strip() for i in row.function_ids.strip().split(',')])
                present_in_bin = id_set & function_id_set
                functions_present = functions_present | present_in_bin
                presents_in_bin.append(len(present_in_bin) > 0)
            function_in_bin = np.all(presents_in_bin)
            row = frame.iloc[0]
            rows.append([row.category, row.subcategory, row.function_name, ', '.join(functions_present),
                         '; '.join(get_ordered_uniques(frame.long_function_name)),
                         '; '.join(get_ordered_uniques(frame.gene_symbol)), bin_name, function_in_bin,
                         '%s: %s' % (row.category, row.function_name)])
    return pd.DataFrame(rows, columns=list(function_heatmap_form.columns) + ['genome', 'present',
                                                                             'category_function_name'])


# TODO: refactor this to handle splitting large numbers of genomes into multiple heatmaps here
def fill_liquor_dfs(annotations, module_nets, etc_module_database, function_heatmap_form, groupby_column='fasta'):
    module_coverage_frame = make_module_coverage_frame(annotations, module_nets, groupby_column)

    # make ETC frame
    etc_coverage_df = make_etc_coverage_df(etc_module_database, annotations, groupby_column)

    # make functional frame
    function_df = make_functional_df(annotations, function_heatmap_form, groupby_column)

    return module_coverage_frame, etc_coverage_df, function_df


def rename_genomes_to_taxa(function_df, labels, mag_order):
    function_df = function_df.copy()
    new_genome_column = [labels[i] for i in function_df['genome']]
    function_df['genome'] = pd.Series(new_genome_column, index=function_df.index)
    mag_order = [labels[i] for i in mag_order]
    return function_df, mag_order




def make_liquor_df(module_coverage_frame, etc_coverage_df, function_df):
    liquor_df = pd.concat([module_coverage_frame.pivot(index='genome', columns='module_name', values='step_coverage'),
                           etc_coverage_df.pivot(index='genome', columns='complex_module_name', values='percent_coverage'),
                           function_df.pivot(index='genome', columns='category_function_name', values='present')],
                          axis=1, sort=False)
    return liquor_df


def get_phylum_and_most_specific(taxa_str):
    taxa_ranks = [i[3:] for i in taxa_str.split(';')]
    phylum = taxa_ranks[1]
    most_specific_rank = TAXONOMY_LEVELS[sum([len(i) > 0 for i in taxa_ranks])-1]
    most_specific_taxa = taxa_ranks[sum([len(i) > 0 for i in taxa_ranks]) - 1]
    if most_specific_rank == 'd':
        return 'd__%s;p__' % most_specific_taxa
    if most_specific_rank == 'p':
        return 'p__%s;c__' % most_specific_taxa
    else:
        return 'p__%s;%s__%s' % (phylum, most_specific_rank, most_specific_taxa)


def make_strings_no_repeats(genome_taxa_dict):
    labels = dict()
    seen = Counter()
    for genome, taxa_string in genome_taxa_dict.items():
        final_taxa_string = '%s_%s' % (taxa_string, str(seen[taxa_string]))
        seen[taxa_string] += 1
        labels[genome] = final_taxa_string
    return labels


@click.command()
@click.option("-i", "--input_file", help="Annotations path")
@click.option("-o", "--output_dir", help="Directory to write summarized genomes")
@click.option("--rrna_path", help="rRNA output from annotation")
@click.option("--trna_path", help="tRNA output from annotation")
@click.option("--groupby_column", help="Column from annotations to group as organism units",
              default='fasta')
@click.option("--custom_distillate", help="Custom distillate form to add your own modules")
@click.option("--custom_pathways", help="Custom pathway form to add your own modules")
@click.option("--custom_function", help="Custom function form to add your own modules")
@click.option("--custom_steps", help="Custom steps form to add your own modules")
@click.option("--replace_forms/--combine_forms", default=False,
              help="In some instances you may want to completely replace the"
              " set of forms used for the distillate with your own. Note"
              " that default forms will be used for any you do not specify"
              " with additional arguments")
@click.option("--genomes_per_product",
              help="Number of genomes per product.html output. Decrease "
                   "value if getting JavaScript Error: Maximum call stack "
                   "size exceeded when viewing product.html in browser.",
              default=1000, type=int)
@click.option("--distillate_gene_names/--distillate_gene_counts", default=False,
              help="Give names of genes instead of counts in genome metabolism summary")
def distill_genomes(input_file:str, trna_path:str=None, rrna_path:str=None,
                    output_dir:str='.', groupby_column:str='fasta',
                    custom_distillate:str=None,
                    custom_pathways:str=None,
                    custom_function:str=None,
                    custom_steps:str=None,
                    replace_forms:bool=False,
                    distillate_gene_names:bool=False,
                    genomes_per_product:bool=1000):
    start_time = datetime.now()

    # read in data
    annotations = pd.read_csv(input_file, sep='\t', index_col=0)
    if 'bin_taxnomy' in annotations:
        annotations = annotations.sort_values('bin_taxonomy')

    if trna_path is None:
        trna_frame = None
    else:
        trna_frame = pd.read_csv(trna_path, sep='\t')
    if rrna_path is None:
        rrna_frame = None
    else:
        rrna_frame = pd.read_csv(rrna_path, sep='\t')

    # read in dbs
    # TODO this is repetitive
    if custom_distillate is None:
        genome_summary_form = pd.read_csv(path.join(DATA_FOLDER, 'genome_summary_form.tsv'), sep='\t')
    else:
        if replace_forms:
            genome_summary_form = pd.read_csv(custom_distillate, sep='\t')
        else:
            genome_summary_form = pd.concat([
                pd.read_csv(path.join(DATA_FOLDER, 'genome_summary_form.tsv'), sep='\t'),
                pd.read_csv(custom_distillate, sep='\t')])

    if custom_steps is None:
        module_step_form = pd.read_csv(path.join(DATA_FOLDER, 'module_step_form.tsv'), sep='\t')
    else:
        if replace_forms:
            module_step_form = pd.read_csv(custom_steps, sep='\t')
        else:
            module_step_form = pd.concat([
                pd.read_csv(path.join(DATA_FOLDER, 'module_step_form.tsv'), sep='\t'),
                pd.read_csv(custom_steps, sep='\t')])

    if custom_function is None:
        function_heatmap_form = pd.read_csv(path.join(DATA_FOLDER, 'function_heatmap_form.tsv'), sep='\t')
    else:
        if replace_forms:
            function_heatmap_form = pd.read_csv(custom_function, sep='\t')
        else:
            function_heatmap_form = pd.concat([
                pd.read_csv(path.join(DATA_FOLDER, 'function_heatmap_form.tsv'), sep='\t'),
                pd.read_csv(custom_function, sep='\t')])

    # TODO Rename etc_module_database to some thing more generic
    if custom_pathways is None:
        etc_module_database = pd.read_csv(path.join(DATA_FOLDER, 'etc_module_database.tsv'), sep='\t')
    else:
        if replace_forms:
            etc_module_database = pd.read_csv(custom_pathways, sep='\t')
        else:
            etc_module_database = pd.concat([
                pd.read_csv(path.join(DATA_FOLDER, 'etc_module_database.tsv'), sep='\t'),
                pd.read_csv(custom_pathways, sep='\t')])

    genome_summary_form = genome_summary_form.drop('potential_amg', axis=1)

    print('%s: Retrieved database locations and descriptions' % (str(datetime.now() - start_time)))

    # make output folder
    # mkdir(output_dir)
    

    # make genome stats
    genome_stats = make_genome_stats(annotations, rrna_frame, trna_frame, groupby_column=groupby_column)
    genome_stats.to_csv(path.join(output_dir, 'genome_stats.tsv'), sep='\t', index=None)
    print('%s: Calculated genome statistics' % (str(datetime.now() - start_time)))

    # make genome metabolism summary
    genome_summary = path.join(output_dir, 'metabolism_summary.xlsx')
    if distillate_gene_names:
        summarized_genomes = fill_genome_summary_frame_gene_names(annotations, genome_summary_form, groupby_column)
    else:
        summarized_genomes = make_genome_summary(annotations, genome_summary_form, trna_frame, rrna_frame,
                                                 groupby_column)
    write_summarized_genomes_to_xlsx(summarized_genomes, genome_summary)
    print('%s: Generated genome metabolism summary' % (str(datetime.now() - start_time)))

    # make liquor
    if 'bin_taxonomy' in annotations:
        genome_order = get_ordered_uniques(annotations.sort_values('bin_taxonomy')[groupby_column])
        # if gtdb format then get phylum and most specific
        if all([i[:3] == 'd__' and len(i.split(';')) == 7 for i in annotations['bin_taxonomy'].fillna('')]):
            taxa_str_parser = get_phylum_and_most_specific
        # else just throw in what is there
        else:
            taxa_str_parser = lambda x: x
        labels = make_strings_no_repeats({row[groupby_column]: taxa_str_parser(row['bin_taxonomy'])
                                         for _, row in annotations.iterrows()})
    else:
        genome_order = get_ordered_uniques(annotations.sort_values(groupby_column)[groupby_column])
        labels = None

    # make module coverage frame
    module_nets = {module: build_module_net(module_df)
                   for module, module_df in module_step_form.groupby('module') if module in HEATMAP_MODULES}

    if len(genome_order) > genomes_per_product:
        module_coverage_dfs = list()
        etc_coverage_dfs = list()
        function_dfs = list()
        # generates slice start and slice end to grab from genomes and labels from 0 to end of genome order
        pairwise_iter = pairwise(list(range(0, len(genome_order), genomes_per_product)) + [len(genome_order)])
        for i, (start, end) in enumerate(pairwise_iter):
            genomes = genome_order[start:end]
            annotations_subset = annotations.loc[[genome in genomes for genome in annotations[groupby_column]]]
            dfs = fill_liquor_dfs(annotations_subset, module_nets, etc_module_database, function_heatmap_form,
                                  groupby_column='fasta')
            module_coverage_df_subset, etc_coverage_df_subset, function_df_subset = dfs
            module_coverage_dfs.append(module_coverage_df_subset)
            etc_coverage_dfs.append(etc_coverage_df_subset)
            function_dfs.append(function_df_subset)
        liquor_df = make_liquor_df(pd.concat(module_coverage_dfs), pd.concat(etc_coverage_dfs), pd.concat(function_dfs))
        liquor_df.to_csv(path.join(output_dir, 'product.tsv'), sep='\t')
    else:
        module_coverage_df, etc_coverage_df, function_df = fill_liquor_dfs(annotations, module_nets,
                                                                           etc_module_database,
                                                                           function_heatmap_form,
                                                                           groupby_column=groupby_column)
        liquor_df = make_liquor_df(module_coverage_df, etc_coverage_df, function_df)
        liquor_df.to_csv(path.join(output_dir, 'product.tsv'), sep='\t')
    print('%s: Generated product heatmap and table' % (str(datetime.now() - start_time)))
    print("%s: Completed distillation" % str(datetime.now() - start_time))


if __name__ == '__main__':
    distill_genomes()
