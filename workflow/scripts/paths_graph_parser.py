"""
Tool to parse the rules tsv, and make a graph.

use:

python3 ./workflow/scripts/paths_graph_parser.py \
    -p ./resources/pathways_refined_test.csv \
    -r ./resources/custom_input_modules/EMERGE_pathways_module.tsv -d \
    ./results/emerge_distilled_filtered/EMERGE_170222_distillate.tsv \
        evaluate  \
            -o results/emerge_distilled_filtered/EMERGE_170222_product_refined.tsv  \
            -t 7

python3 ./workflow/scripts/paths_graph_parser.py \
    -p ./resources/pathways_refined_test.csv \
    -r ./resources/custom_input_modules/EMERGE_pathways_module.tsv \
    -d ./results/emerge_distilled_filtered/EMERGE_170222_distillate.tsv \
        test
"""
import re
from functools import reduce, partial
from multiprocessing import Pool
import pandas as pd
import numpy as np
import networkx as nx
import graphviz
import click
import pytest


"""
import os
os.system("python3 ./workflow/scripts/paths_graph_parser.py -p ./resources/pathways_refined_test.csv -r ./resources/custom_input_modules/EMERGE_pathways_module.tsv -d ./results/emerge_distilled_filtered/EMERGE_170222_distillate.tsv evaluate  -o results/emerge_distilled_filtered/EMERGE_170222_product_refined.tsv  -t 7")
os.system("python3 ./workflow/scripts/paths_graph_parser.py -p ./resources/pathways_refined_test.csv -r ./resources/custom_input_modules/EMERGE_pathways_module.tsv -d ./results/emerge_distilled_filtered/EMERGE_170222_distillate.tsv test")
#"""



def parse_ands(logic:str):
    depth = 0
    for i, c in enumerate(logic):
        if c == "(":
            depth += 1
        if c == ")":
            depth -= 1
        if c == "+" and depth == 0:
            return parse_ands(logic[:i]) + parse_ands(logic[i+1:])
    return [logic]


def parse_ors(logic:str):
    depth = 0
    for i, c in enumerate(logic):
        if c == "(":
            depth += 1
        if c == ")":
            depth -= 1
        if c == "," and depth == 0:
            return parse_ors(logic[:i]) + parse_ors(logic[i+1:])
    return [logic]


def node_color(val):
    if pd.isnull(val):
        return 'grey'
    if val:
        return 'green'
    return 'red'


def ko_set(G, node):
    if G.nodes[node]['type'] == 'ko':
        return node
        # return G.nodes[node]['display']
    if G.nodes[node]['type'] in ['path_or', 'react_or']:
        kids = [ko_set(G, i) for i in G.successors(node)
                if G.nodes[i]['present']]
        lists = [j for i in kids if isinstance(i, list) for j in i]
        kids = lists + [i for i in kids if not isinstance(i, list)]
        kos = [{i} if isinstance(i, str) else i for i in kids ]
        return kos
    if G.nodes[node]['type'] in ['path_and', 'react_and']:
        kids = [ko_set(G, i) for i in G.successors(node)]
        kos = [{i for i in kids if  isinstance(i, str)}]
        kids = [i for i in kids if  isinstance(i, list)]
        kids.append(kos)
        sets = reduce(lambda x, y: [i | j for i in x for j in y], kids)
        return sets
    if len(nexts := [i for i in G.successors(node)]) == 1:
        if G.nodes[nexts[0]]['present']:
            kids = ko_set(G, nexts[0])
            if isinstance(kids, list):
                return kids
            elif isinstance(kids,str):
                return [{kids}]
        else:
            return [set()]
    else:
        raise ValueError('Only ands and ors can have multiple successors')


class RuleParser():


    def __init__(
        self,
        dist:pd.Series,
        present_pct:float=0.7,
        ko_pct:float=0.6
        ):
        self.dist = dist
        self.genome = self.dist.name
        self.sG = None
        self.pG = None
        self.present_pct = present_pct
        self.ko_pct = ko_pct


    def check_signature(self):
       return self.sG.nodes[self.genome]['value']


    def make_signiture_graph(self, signature_definition:str):
        self.sG = nx.DiGraph()
        self.sG.add_node(self.genome, display=self.genome, type='root', value=None)
        self._parse_signature(signature_definition, parent=self.genome)
        only_node = [i for i in self.sG.successors(self.genome)][0]
        self.sG.nodes[self.genome]['value'] = self.sG.nodes[only_node]['value']


    def check_presance(self):
        return self.pG.nodes[self.genome]['present'] >= self.present_pct


    def get_presance(self):
        return self.pG.nodes[self.genome]['present']


    def check_ko_percent(self):
        if not self.pG.nodes[self.genome]['present']:
            return False
        ko_opt = ko_set(self.pG, self.genome)
        for i in ko_opt:
            ave = np.average([self.pG.nodes[k]['present'] for k in i])
            if ave >= self.ko_pct:
                return True
        return False


    def get_ko_percent(self):
        if self.pG.nodes[self.genome]['present'] == 0 :
            return 0
        ko_opt = ko_set(self.pG, self.genome)
        return np.max([
            np.average([self.pG.nodes[k]['present'] for k in i])
            for i in ko_opt])


    def get_ko_set(self):
        if self.pG.nodes[self.genome]['present'] == 0 :
            return []
        ko_opt = ko_set(self.pG, self.genome)
        return [{self.pG.nodes[k]['display'] for k in i} for i in ko_opt]


    def make_path_graph(self, pathway_definition:str, reactn:pd.DataFrame):
        self.reactn = reactn
        if not reactn.index.is_unique:
            raise ValueError("Index is not unique, make it so!")
        self.pG = nx.DiGraph()
        self.pG.add_node(self.genome, display=self.genome, type='root', present=None)
        self.parse_paths(pathway_definition, 0, parent=self.genome)
        if len(nxt:=[i for i in self.pG.successors(self.genome)]) != 1:
            raise ValueError("there is something wrong with the graph")
        only_node = [i for i in self.pG.successors(self.genome)][0]
        self.pG.nodes[self.genome]['present'] = self.pG.nodes[only_node]['present']


    def parse_paths(self, logic:str, num, parent:str=None):
        if logic.isnumeric():
            self._parse_path_id(logic, num, parent)
            return
        if len(ors := parse_ors(logic)) > 1:
            self._parse_path_or(ors, num, parent)
            return
        if len(ands := parse_ands(logic))> 1:
            self._parse_path_and(ands, num, parent)
            return
        if logic.startswith('(') and logic.endswith(')'):
            self.parse_paths(logic[1:-1], num, parent=parent)
            return
        raise ValueError("Parsing Fail")



    def _ko_func(self, ko:str):
        return self.dist.loc[ko] > 0


    def _parse_signature(self, logic:str, parent:str):
        if re.match(r'^K[0-9,\.]+$', logic):
            id = logic
            self.sG.add_node(id, display=id, type='ko',
                           value=self._ko_func(logic))
            self.sG.add_edge(parent, id)
            return
        if len(ors := parse_ors(logic)) > 1:
            id = f"{parent}-or"
            self.sG.add_node(id, display="or", type='or', value=None)
            self.sG.add_edge(parent, id)
            for i  in ors:
                self._parse_signature(i, parent=id)
            self.sG.nodes[id]['value'] = np.any(
                [self.sG.nodes[i]['value'] for i in self.sG.successors(id)])
            return
        if len(ands := parse_ands(logic))> 1:
            id = f"{parent}-and"
            self.sG.add_node(id, display="and", type='and', value=None)
            self.sG.add_edge(parent, id)
            for i  in ands:
                self._parse_signature(i, parent=id)
            self.sG.nodes[id]['value'] = np.all(
                [self.sG.nodes[i]['value'] for i in self.sG.successors(id)])
            return
        if logic.startswith('(') and logic.endswith(')'):
            self._parse_signature(logic[1:-1], parent=parent)
            return
        self.sG.add_node(f"error = {logic}", display="error", type='error',
                       value=None)
        self.sG.add_edge(parent, f"error = {logic}")


    def _parse_reactions(self, logic:str, num:int, parent:str=None):
        if re.match(r'^K[0-9,\.]+$', logic):
            id = f"{parent}-{logic}[{num}]"
            self.pG.add_node(id, display=logic, type='ko',
                           present=self._ko_func(logic))
            self.pG.add_edge(parent, id)
            return
        if len(ors := parse_ors(logic)) > 1:
            id = f"{parent}-or[{num}]"
            self.pG.add_node(id, display="or",
                                 type='react_or', present=None)
            self.pG.add_edge(parent, id)
            for n, i  in enumerate(ors):
                self._parse_reactions(i, n, parent=id)
            self.pG.nodes[id]['present'] = np.any(
                [self.pG.nodes[i]['present']
                 for i in self.pG.successors(id)])
            return
        if len(ands := parse_ands(logic))> 1:
            id = f"{parent}-and[{num}]"
            self.pG.add_node(id, display="and",
                                 type='react_and', present=None)
            self.pG.add_edge(parent, id)
            for n, i  in enumerate(ands):
                self._parse_reactions(i, n, parent=id)
            # Note that this is a rather large assumption!!!!
            # TODO should this be np.all???
            # DONE yes it should be np.any
            self.pG.nodes[id]['present'] = np.any(
                [self.pG.nodes[i]['present']
                 for i in self.pG.successors(id)])
            return
        if logic.startswith('(') and logic.endswith(')'):
            self._parse_reactions(logic[1:-1], num, parent=parent)
            return
        raise ValueError("Parsing Fail")


    def _parse_path_or(self, ors, num, parent):
        id = f"{parent}-or[{num}]"
        self.pG.add_node(id, display="or",
                             type='path_or', present=None)
        self.pG.add_edge(parent, id)
        for n, i  in enumerate(ors):
            self.parse_paths(i, n, parent=id)
        nexts = list(self.pG.successors(id))
        self.pG.nodes[id]['present'] = np.max([
            self.pG.nodes[i]['present'] for i in nexts])


    def _parse_path_and(self, ands, num, parent):
        id = f"{parent}-and[{num}]"
        self.pG.add_node(id,
                             display="and",
                             type='path_and',
                             present=None)
        self.pG.add_edge(parent, id)
        for n, i  in enumerate(ands):
            self.parse_paths(i, n, parent=id)
        nexts = list(self.pG.successors(id))
        if np.all([self.pG.nodes[i]['type'] == 'react_id' for i in nexts]):
            self.pG.nodes[id]['present'] = \
                np.average([self.pG.nodes[i]['present'] for i in nexts])
        else:
            raise ValueError('We have no way to deal with nested ands')


    def _parse_path_id(self, react, num, parent):
            id = f"{parent}-{react}[{num}]"
            self.pG.add_node(id, display=react,
                                 type='react_id', present=None)
            self.pG.add_edge(parent, id)
            self._parse_reactions(self.reactn.loc[react], 0, parent=id)
            self.pG.nodes[id]['present'] = np.any(
                [self.pG.nodes[i]['present']
                 for i in self.pG.successors(id)])


    def _add_to_dot(self, node, graph):
        # self.dot.node(node, self.G[node].get('display'))
        for i in graph.successors(node):
            self.dot.node(
                i,
                f"{graph.nodes[i]['display']}",
                color=node_color(graph.nodes[i].get('value')))
            self.dot.edge(i, node)
            self._add_to_dot(i, graph)


    # def show(self):
    #     if not os.path.exists('figures'):
    #         os.mkdir('figures')
    #     self._show_signature_graph()
    #     self._show_path_graph()


    def _show_signature_graph(self):
        self.dot = graphviz.Digraph(strict=True)
        self.dot.node(
            self.genome,
            f"{self.genome}",
            color=node_color(self.sG.nodes[self.genome]['value']))
        self._add_to_dot(self.genome, self.sG)
        self.dot.render(view=True, filename=f"{self.genome}_signiture", directory='figures')


    def _add_to_path_dot(self, node):
        # self.dot.node(node, self.G[node].get('display'))
        for i in self.pG.successors(node):
            self.dot.node(
                i,
                f"{self.pG.nodes[i]['display']}",
                color=node_color(self.pG.nodes[i].get('present')))
            self.dot.edge(i, node)
            self._add_to_path_dot(i)


    def _show_path_graph(self):
        self.dot = graphviz.Digraph(strict=True)
        self.dot.node(
            self.genome,
            f"{self.genome}",
            color=node_color(self.pG.nodes[self.genome]['present']))
        self._add_to_path_dot(self.genome)
        self.dot.render(view=True, filename=f"{self.genome}_path", directory='figures')

@click.group()
@click.option('-p', '--pathways_file', type=click.Path(exists=True))
@click.option('-r', '--reaction_file', type=click.Path(exists=True))
@click.option('-d', '--distillate_file', type=click.Path(exists=True))
@click.pass_context
def met_paths(ctx, pathways_file:str, reaction_file:str, distillate_file:str):
    ctx.ensure_object(dict)
    click.echo(f"Reading data")
    pathwy = (pd.read_csv(pathways_file, sep='\t',
                          dtype=str)
              .drop([
                  'note',
                  'signature_gene_name'
              ], axis=1)
              )
    pathwy.set_index(
        pathwy['pathway'] + '-' + pathwy['subpathway'],
        inplace=True)
    reactn = (
        pd.read_csv(
            reaction_file,
            sep='\t',
            dtype=str)
        .assign(module_num=(
            lambda x: x['module_name']
            .str
            .split(' -', expand=True)[1]))
        .drop(['module_id', 'module_name'], axis=1)
        .set_index(['complex', 'module_num'])
    )
    dist = pd.read_csv(
        distillate_file,
        sep='\t', index_col='gene_id').drop([
            'gene_description',
            'module',
            'header',
            'subheader'
        ], axis=1)
    ctx.obj['DATA'] = pathwy, reactn, dist


def evaluate_pathway(pathway:pd.Series, tree:RuleParser,
                     reactions:pd.DataFrame):
    if pd.notna(pathway['signature_definition']):
        tree.make_signiture_graph(pathway['signature_definition'])
        if not tree.check_signature():
            return False
    tree.make_path_graph(
        pathway['reaction'],
        reactions.loc[pathway['pathway'], 'definition'])
    if not tree.check_presance():
        return False
    return tree.check_ko_percent()


def evaluate_distilate(dist_col, reactions, pathways):
    tree = RuleParser(dist_col)
    apply_funk = partial(evaluate_pathway, tree=tree, reactions=reactions)
    results = pd.Series(
        pathways.apply(apply_funk, axis=1), name=dist_col.name)
    return results


@met_paths.command()
@click.pass_context
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-t', '--threads', type=int)
def evaluate(ctx, output:str, threads:int):
    pool = Pool(threads)
    pathways, reactions, distillate = ctx.obj['DATA']
    pool_funk = partial(
        evaluate_distilate,
        pathways=pathways,
        reactions=reactions
    )
    with Pool(threads) as pool:
        results = pd.concat(pool.map(
            pool_funk,
            [col for _, col in distillate.iteritems()]
        ), axis=1).T
    results.index.name = 'genome'
    results.to_csv(output, sep="\t")


@met_paths.command()
@click.pass_context
def test(ctx):
    pathways, reactions, distillate = ctx.obj['DATA']
    pool_funk = partial(
        evaluate_distilate,
        pathways=pathways,
        reactions=reactions
    )
    [pool_funk(col) for _, col in distillate.iteritems()]


if __name__ == '__main__':
    met_paths()


def test_signiture():
    test_dist =pd.DataFrame({
        'K01198': [4, 1, 0],
        'K01625': [0, 1, 0],
        'K00008': [3, 0, 0],
        'K01783': [1, 1, 1]},
       index=['PNBF01', 'PMDY01', 'PKVC01']
       ).T
    test_n_out = [
            ['K01198,K01625',               'PNBF01', True ],
            ['K01198+K01625',               'PNBF01', False],
            ['K01198,K01625',               'PMDY01', True ],
            ['K01198+K01625',               'PMDY01', True ],
            ['K01198,K01625',               'PKVC01', False],
            ['K01783+K00008,K01198+K01625', 'PKVC01', False],
            ['K01783,K01198+K01625+K00008', 'PKVC01', True ]]
    for logic, genom, output in test_n_out:
        tree = RuleParser(test_dist[genom])
        tree.make_signiture_graph(logic)
        assert tree.check_signature() == output, "problem in key parser"

def test_pathways():
    test_dist = pd.DataFrame({
        'K00001': [4, 1, 0],
        'K00002': [0, 1, 0],
        'K00003': [3, 2, 0],
        'K00004': [1, 1, 1],
        'K00005': [0, 4, 0]},
       index=['G01', 'G02', 'G03']
       ).T
    reactions = (pd.DataFrame(
        [
            ["EMERGE0", "nitrate->nitrite -1", "nitrogen_redox",
            "(K00001+K00003),(K00005+K00005)"],
           ["EMERGE1", "nitrite->ammonium -2", "nitrogen_redox",
            "(K00005+K00005),(K00005+K00005)"],
           ["EMERGE6", "dinitrogen->ammonium -3", "nitrogen_redox",
            "(K00005+K00005+K00005+K00005+K00005+K00005+K00005+K00005),"
            "(K00001+K00003+K00004+K00001+K00002+K00005+K00005),"
            "(K00002+K00001+K00003+K00005+K00005+K00005)"],
            ["EMERGE0", "nitrate->nitrite -4", "nitrogen_redox",
            "(K00005+K00005),(K00005+K00005)"],
           ["EMERGE1", "nitrite->ammonium -5", "nitrogen_redox",
            "(K00004),(K00003)"],
        ],
        columns=["module_id", "module_name", "complex", "definition"]
        ).assign(
            module_num=(lambda x: x['module_name'].str.split(' -', expand=True)[1])
        )
        .drop(['module_id', 'module_name'], axis=1)
        .set_index(['complex', 'module_num'])
    )
    reaction_n = reactions.loc['nitrogen_redox', 'definition']
    tree = RuleParser(test_dist['G01'])
    tree.make_path_graph('1', reaction_n)
    assert tree.check_presance()
    assert tree.get_presance() == 1
    assert tree.check_ko_percent()
    assert len(tree.get_ko_set()) == 1
    assert tree.get_ko_percent() == 1
    tree.make_path_graph('1+1+1+2', reaction_n)
    assert tree.check_presance()
    assert tree.get_presance() == 0.75
    tree.make_path_graph('1+2', reaction_n)
    assert not tree.check_presance()
    assert tree.get_presance() == 0.5
    assert tree.check_ko_percent()
    assert len(tree.get_ko_set()) == 1
    assert tree.get_ko_set() == [{'K00001', 'K00003'}]
    assert tree.get_ko_percent() == 1
    tree.make_path_graph('2', reaction_n)
    assert not tree.check_presance()
    assert tree.get_presance() == 0
    tree.make_path_graph('1+5', reaction_n)
    assert tree.get_ko_set() == [{'K00001', 'K00003', 'K00004'}, {'K00001', 'K00003'}]
    assert tree.get_ko_percent() == 1
    tree.make_path_graph('4', reaction_n)
    assert tree.get_ko_set() == []
    assert tree.get_ko_percent() == 0
    tree = RuleParser(test_dist['G02'])
    tree.make_path_graph('3+1+5', reaction_n)
    assert len(tree.get_ko_set()) == 12
    assert tree.get_ko_percent() == 1

