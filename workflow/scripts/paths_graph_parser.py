"""Tool to parse the rules tsv, and make a graph."""
import re
import os
from functools import reduce, partial
from multiprocessing import Pool
import pandas as pd
import numpy as np
import networkx as nx
import graphviz
import click



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


def get_ko_set(G, node):
    if G.nodes[node]['type'] == 'ko':
        return node
        # return G.nodes[node]['display']
    if G.nodes[node]['type'] in ['path_or', 'react_or']:
        kids = [get_ko_set(G, i) for i in G.successors(node)]
        lists = [j for i in kids if isinstance(i, list) for j in i]
        kids = lists + [i for i in kids if not isinstance(i, list)]
        kos = [{i} if isinstance(i, str) else i for i in kids ]
        return kos
    if G.nodes[node]['type'] in ['path_and', 'react_and']:
        kids = [get_ko_set(G, i) for i in G.successors(node)]
        kos = [{i for i in kids if  isinstance(i, str)}]
        kids = [i for i in kids if  isinstance(i, list)]
        kids.append(kos)
        sets = reduce(lambda x, y: [i | j for i in x for j in y], kids)
        return sets
    if len(nexts := [i for i in G.successors(node)]) == 1:
         if G.nodes[nexts[0]]['present']:
             kids = get_ko_set(G, nexts[0])
             if isinstance(kids, list):
                 return kids
             elif isinstance(kids,str):
                 return [{kids}]
         else:
             return []
    else:
        raise ValueError('Only ands and ors can have multiple successors')


class RuleParser():


    def __init__(
        self,
        dist:pd.Series,
        genome:str,
        pathwy:pd.Series,
        reactn:pd.DataFrame,
        present_pct:float=0.7,
        ko_pct:float=0.6
        ):
        self.dist = dist
        self.genome = genome
        self.pathwy = pathwy
        self.reactn = reactn
        self.signature_g = nx.DiGraph()
        self.path_g = nx.DiGraph()
        self.present_pct = present_pct
        self.ko_pct = ko_pct
        # parse the rules
        # parse the paths
        if pd.notnull(sig := pathwy['signature_definition']):
            self.make_signature_graph(self.signature_g,
                                     sig,
                                     self.genome)
        self.make_path_graph(pathwy['reaction'],
                                 self.genome)


    def make_signature_graph(self, graph, logic, id):
        graph.add_node(id, display="and", type='and', value=None)
        self.parse_signature(logic,
                             graph,
                             parent=id)
        only_node = [i for i in graph.successors(id)][0]
        graph.nodes[id]['value'] = graph.nodes[only_node]['value']


    def make_path_graph(self, logic, id):
        self.path_g.add_node(id, display="and", type='root', present=None)
        self.parse_paths(logic, 0, parent=id)
        only_node = [i for i in self.path_g.successors(id)][0]
        self.path_g.nodes[id]['present'] = self.path_g.nodes[only_node]['present']


    def ko_func(self, ko:str):
        return self.dist.loc[ko] > 0


    def parse_signature(self, logic:str, graph:nx.DiGraph,
                        parent:str = None):
        if re.match(r'^K[0-9,\.]+$', logic):
            id = logic
            graph.add_node(id, display=id, type='ko',
                           value=self.ko_func(logic))
            graph.add_edge(parent, id)
            return
        if len(ors := parse_ors(logic)) > 1:
            id = f"{parent}-or"
            graph.add_node(id, display="or", type='or', value=None)
            graph.add_edge(parent, id)
            for i  in ors:
                self.parse_signature(i, graph, parent=id)
            graph.nodes[id]['value'] = np.any(
                [graph.nodes[i]['value'] for i in graph.successors(id)])
            return
        if len(ands := parse_ands(logic))> 1:
            id = f"{parent}-and"
            graph.add_node(id, display="and", type='and', value=None)
            graph.add_edge(parent, id)
            for i  in ands:
                self.parse_signature(i, graph, parent=id)
            graph.nodes[id]['value'] = np.all(
                [graph.nodes[i]['value'] for i in graph.successors(id)])
            return
        if logic.startswith('(') and logic.endswith(')'):
            self.parse_signature(logic[1:-1], graph, parent=parent)
            return
        graph.add_node(f"error = {logic}", display="error", type='error',
                       value=None)
        graph.add_edge(parent, f"error = {logic}")


    def parse_reactions(self, logic:str, num:int, parent:str=None):
        if re.match(r'^K[0-9,\.]+$', logic):
            id = f"{parent}-{logic}[{num}]"
            self.path_g.add_node(id, display=logic, type='ko',
                           present=self.ko_func(logic))
            self.path_g.add_edge(parent, id)
            return
        if len(ors := parse_ors(logic)) > 1:
            id = f"{parent}-or[{num}]"
            self.path_g.add_node(id, display="or",
                                 type='react_or', present=None)
            self.path_g.add_edge(parent, id)
            for n, i  in enumerate(ors):
                self.parse_reactions(i, n, parent=id)
            self.path_g.nodes[id]['present'] = np.any(
                [self.path_g.nodes[i]['present']
                 for i in self.path_g.successors(id)])
            return
        if len(ands := parse_ands(logic))> 1:
            id = f"{parent}-and[{num}]"
            self.path_g.add_node(id, display="and",
                                 type='react_and', present=None)
            self.path_g.add_edge(parent, id)
            for n, i  in enumerate(ands):
                self.parse_reactions(i, n, parent=id)
            # Note that this is a rather large assumption!!!!
            # TODO should this be np.all???
            # DONE yes it should be np.any
            self.path_g.nodes[id]['present'] = np.any(
                [self.path_g.nodes[i]['present']
                 for i in self.path_g.successors(id)])
            return
        if logic.startswith('(') and logic.endswith(')'):
            self.parse_reactions(logic[1:-1], num, parent=parent)
            return
        raise ValueError("Parsing Fail")

    def parse_path_or(self, ors, num, parent):
        id = f"{parent}-or[{num}]"
        self.path_g.add_node(id, display="or",
                             type='path_or', present=None)
        self.path_g.add_edge(parent, id)
        for n, i  in enumerate(ors):
            self.parse_paths(i, n, parent=id)
        nexts = self.path_g.successors(id)
        self.path_g.nodes[id]['present'] = np.max([
            self.path_g.nodes[i]['present'] for i in nexts])


    def parse_path_and(self, ands, num, parent):
        id = f"{parent}-and[{num}]"
        self.path_g.add_node(id,
                             display="and",
                             type='path_and',
                             present=None)
        self.path_g.add_edge(parent, id)
        for n, i  in enumerate(ands):
            self.parse_paths(i, n, parent=id)
        nexts = self.path_g.successors(id)
        if np.all([self.path_g.nodes[i]['type'] == 'react_id'
                   for i in nexts]):
            self.path_g.nodes[id]['present'] = \
                np.average([i['present'] == 'react_id' for i in nexts]) >= self.present_pct
        else:
            raise ValueError('We have no way to deal with nested ands')

    def parse_path_id(self, react, num, parent):
            id = f"{parent}-{react}[{num}]"
            self.path_g.add_node(id, display=react,
                                 type='react_id', present=None)
            self.path_g.add_edge(parent, id)
            self.parse_reactions(self.reactn.loc[react], 0, parent=id)
            self.path_g.nodes[id]['present'] = np.any(
                [self.path_g.nodes[i]['present']
                 for i in self.path_g.successors(id)])

    def parse_paths(self, logic:str, num, parent:str=None):
        if logic.isnumeric():
            self.parse_path_id(logic, num, parent)
            return
        if len(ors := parse_ors(logic)) > 1:
            self.parse_path_or(ors, num, parent)
            return
        if len(ands := parse_ands(logic))> 1:
            self.parse_path_and(ands, num, parent)
            return
        if logic.startswith('(') and logic.endswith(')'):
            self.parse_paths(logic[1:-1], num, parent=parent)
            return
        raise ValueError("Parsing Fail")


    def add_to_dot(self, node, graph):
        # self.dot.node(node, self.G[node].get('display'))
        for i in graph.successors(node):
            self.dot.node(
                i,
                f"{graph.nodes[i]['display']}",
                color=node_color(graph.nodes[i].get('value')))
            self.dot.edge(i, node)
            self.add_to_dot(i, graph)

    def show(self, name):
        if name == 'signature':
            self._show_signature_graph()
        if name == 'path':
            self._show_path_graph()


    def _show_signature_graph(self):
        self.dot = graphviz.Digraph(strict=True)
        self.dot.node(
            self.genome,
            f"{self.genome}",
            color=node_color(self.signature_g.nodes[self.genome]['value']))
        self.add_to_dot(self.genome, self.signature_g)
        self.dot.render(view=True)

    def check_presance(self):
        return self.path_g.nodes[self.genome]['present']

    def check_signature(self):
        if pd.isnull(self.pathwy['signature_definition']):
            return True
        return self.signature_g.nodes[self.genome]['value']

    def check_ko_percent(self):
        if not self.path_g.nodes[self.genome]['present']:
            return False
        ko_opt = get_ko_set(self.path_g, self.genome)
        for i in ko_opt:
            ave = np.average([self.path_g.nodes[k]['present'] for k in i])
            if ave >= self.ko_pct :
                return True


    def add_to_path_dot(self, node):
        # self.dot.node(node, self.G[node].get('display'))
        for i in self.path_g.successors(node):
            self.dot.node(
                i,
                f"{self.path_g.nodes[i]['display']}",
                color=node_color(self.path_g.nodes[i].get('present')))
            self.dot.edge(i, node)
            self.add_to_path_dot(i)


    def _show_path_graph(self):
        self.dot = graphviz.Digraph(strict=True)
        self.dot.node(
            self.genome,
            f"{self.genome}",
            color=node_color(self.path_g.nodes[self.genome]['present']))
        self.add_to_path_dot(self.genome)
        self.dot.render(view=True)


# @met_paths.command()
# @met_paths.option('-a', '--adjectives', type=click.Path(exists=True))
# def show():
#     pass


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
              ], axis=1))
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

def check_rules(args, reactn, dist):
    row, col = args
    # print(row, col)
    test = RuleParser(
        dist[col],
        col,
        row,
        reactn.loc[row['pathway']]['definition'],
    )
    if test.check_ko_percent() \
            and test.check_signature() \
            and test.check_presance():
        return pd.DataFrame({
            'genome': [col],
            'path': [f"{row['pathway']}-{row['subpathway']}"],
            'value': True})
    else:
        return pd.DataFrame({
            'genome': [col],
            'path': [f"{row['pathway']}-{row['subpathway']}"],
            'value': False})

@met_paths.command()
@click.pass_context
def evaluate(ctx):
    pool = Pool(8)
    pathwy, reactn, dist = ctx.obj['DATA']
    results = (
        pd.concat(
            pool.map(
                partial(check_rules, reactn=reactn, dist=dist),
                [(row, col)
                 for _, row in pathwy.iterrows()
                 for col in dist.columns]
            )
        )
    .reset_index()
    .groupby(['path', 'genome'])['value'].max().unstack()
    )
    results.to_csv("")
    # results = [partial(check_rules, reactn=reactn, dist=dist)((row, col))
    #            for _, row in pathwy.iterrows() for col in dist.columns]


if __name__ == '__main__':
    met_paths()
#test.show('signature')

"""
import os
os.system(
    "python3 ./graph_parser.py -p ./stage2_paths/pathways2.tsv"
    " -r ./custom_input_modules/EMERGE_pathways_module.tsv"
    " -d ./emerge_filtered_output/EMERGE_170222_distillate.tsv evaluate "
)


[len(i) for i in kids]
len(kids)
ko_opt
len(ko_opt)
continue
exit
#"""
