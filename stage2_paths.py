import re
import pandas as pd
import numpy as np
import pytest
from datetime import date
from multiprocessing import Pool


def get_value(pathway:str, step:str, values:pd.Series):
    column = [i for i in values.index
            if re.match(f"^Complex {pathway}:.*-{step}$", i)]
    if len(column) != 1:
        raise ValueError("The column name is ambiguous for"
                         f" pathway={pathway}, step={step}")
    return [values[column[0]] >= 0.6]



def evaluat(logic:str, pathway:str,  values:pd.Series):
    l = 0 # level
    for i, c in enumerate(logic):
        if l == 0:
           if c == ',':
               return [np.max([
                   np.average(evaluat(logic[:i], pathway, values)),
                   np.average(evaluat(logic[i+1:], pathway, values))
                             ])]
        if c == '(':
            l += 1
        if c == ')':
            l -= 1
    if l != 0:
        raise ValueError('Bad use of (), you may have an unclosed set')
    for i, c in enumerate(logic):
        if l == 0:
           if c == '+':
               return evaluat(logic[:i], pathway, values) + \
                   evaluat(logic[i+1:], pathway, values)
        if c == '(':
            l += 1
        if c == ')':
            l -= 1
    if l != 0:
        raise ValueError('Bad use of (), you may have an unclosed set')
    if logic.startswith('(') and logic.endswith(')'):
        return evaluat(logic[1:-1], pathway, values)
    return get_value(pathway=pathway, step=logic, values=values)

def evaluat_signature(logic, dist, genome):
    if ',' in logic:
        return np.any([evaluat_signature(i, dist, genome) for i in logic.split(',')])
    if '+' in logic:
        return np.all([evaluat_signature(i, dist, genome) for i in logic.split('+')])
    return dist.loc[logic,genome] > 0

def evaluat_at_percent(logic:str, pathway:str, values:pd.Series,
                       percent:float, sig_ko, dist:pd.DataFrame=None):
    evaluation = np.average(evaluat(logic, pathway, values)) >= percent
    if evaluation and dist is not None and not pd.isnull(sig_ko):
        return evaluat_signature(str(sig_ko), dist, values['genome'])
    return evaluation

def test_evaluat_signiture():
    test_dist =pd.DataFrame({
        'K01198': [4, 1, 0],
        'K01625': [0, 1, 0],
        'K00008': [3, 0, 0],
        'K01783': [1, 1, 1]},
       index=['PNBF01', 'PMDY01', 'PKVC01']
       ).T
    assert evaluat_signature('K01198,K01625', test_dist, 'PNBF01') == True
    assert evaluat_signature('K01198+K01625', test_dist, 'PNBF01') == False
    assert evaluat_signature('K01198,K01625', test_dist, 'PMDY01') == True
    assert evaluat_signature('K01198+K01625', test_dist, 'PMDY01') == True
    assert evaluat_signature('K01198,K01625', test_dist, 'PKVC01') == False
    assert evaluat_signature('K01783+K00008,K01198+K01625', test_dist, 'PKVC01') == False
    assert evaluat_signature('K01783,K01198+K01625+K00008', test_dist, 'PKVC01') == True

def test_get_value():
    t_values = pd.Series(
         {"Complex nitrogen_redox: nitrate->nitrite -0": 0.0,
          "Complex nitrogen_redox: nitrite->ammonium -1": 0.6},
         name="20100900_E1D_1")
    assert get_value(pathway='nitrogen_redox', step='1',
                         values=t_values) == [True]
    assert get_value(pathway='nitrogen_redox', step='0',
                         values=t_values) == [False]

def test_evaluat():
    values = pd.Series({'Complex nitrogen_redox: any -1': 0.6,
                    'Complex nitrogen_redox: any -2': 0.25,
                    'Complex nitrogen_redox: any -3': 0.8,
                    })
    assert evaluat('2', 'nitrogen_redox', values) == [False]
    assert evaluat('1', 'nitrogen_redox', values) == [True]
    assert evaluat('1+2', 'nitrogen_redox', values) == [True, False]
    assert evaluat('1,2', 'nitrogen_redox', values) == [True]
    assert isinstance(evaluat('1+2', 'nitrogen_redox', values),
                      list)
    assert evaluat('(1+3)+(1+2)', 'nitrogen_redox',
                   values) == [True, True, True, False]
    assert evaluat('(1+2),(1+2)', 'nitrogen_redox',
                   values) == [0.5]
    assert isinstance(evaluat('(1+2),(1+2)', 'nitrogen_redox',
                   values), list)

def test_evaluat_signiture():
test_dist =pd.DataFrame({
    'K01198': [4, 1, 0],
    'K01625': [0, 1, 0],
    'K00008': [3, 0, 0],
    'K01783': [1, 1, 1]},
   index=['PNBF01', 'PMDY01', 'PKVC01']
   ).T
assert evaluat_signature('K01198,K01625', test_dist, 'PNBF01') == True
assert evaluat_signature('K01198+K01625', test_dist, 'PNBF01') == False
assert evaluat_signature('K01198,K01625', test_dist, 'PMDY01') == True
assert evaluat_signature('K01198+K01625', test_dist, 'PMDY01') == True
assert evaluat_signature('K01198,K01625', test_dist, 'PKVC01') == False
assert evaluat_signature('K01783+K00008,K01198+K01625', test_dist, 'PKVC01') == False
assert evaluat_signature('K01783,K01198+K01625+K00008', test_dist, 'PKVC01') == True


def test_evaluat_at_percent():
    values = pd.Series({
       'Complex nitrogen_redox: any -1': 0.6,
       'Complex nitrogen_redox: any -2': 0.25,
       'Complex nitrogen_redox: any -3': 0.8,
       'genome': 'PNBF01'
       })
    test_dist =pd.DataFrame({
        'K01198': [4, 1, 0],
        'K01625': [0, 1, 0],
        'K00008': [3, 0, 0],
        'K01783': [1, 1, 1]},
       index=['PNBF01', 'PMDY01', 'PKVC01']
       ).T
    assert evaluat_at_percent('2', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == False
    assert evaluat_at_percent('1', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == True
    assert evaluat_at_percent('1+2', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == False
    assert evaluat_at_percent('1,2', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == True
    assert evaluat_at_percent('1,2', 'nitrogen_redox', values, 0.7, 'K01198+K01625', test_dist) == False
    assert evaluat_at_percent('1+3+2', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == False
    assert evaluat_at_percent('(1+2),(1+2)', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == False
    assert evaluat_at_percent('(1+2),(1+2)', 'nitrogen_redox', values, 0.7, 'K01198+K01625', test_dist) == False
    assert evaluat_at_percent('1+3', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == True
    assert evaluat_at_percent('(1+3)+(1+2)', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == True
    assert evaluat_at_percent('(1,2)+(1,2)', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == True
    assert evaluat_at_percent('((1+3)+(3+2)),((1,2)+(1,2))', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == True
    assert evaluat_at_percent('((1+3)+(1+2))+((1+3)+(1+2)),((1,2)+(1,2))', 'nitrogen_redox', values, 0.7, 'K01198,K01625', test_dist) == True
    assert evaluat_at_percent('(1+3)+(1+2)', 'nitrogen_redox', values, 0.7, 'K01198+K01625', test_dist) == False
    assert evaluat_at_percent('(1,2)+(1,2)', 'nitrogen_redox', values, 0.7, 'K01198+K01625', test_dist) == False
    assert evaluat_at_percent('((1+3)+(3+2)),((1,2)+(1,2))', 'nitrogen_redox', values, 0.7, 'K01198+K01625', test_dist) == False
    assert evaluat_at_percent('((1+3)+(1+2))+((1+3)+(1+2)),((1,2)+(1,2))', 'nitrogen_redox', values, 0.7, 'K01198+K01625', test_dist) == False

test_evaluat_signiture()
test_evaluat_at_percent
if __name__ == '__main__':
    l2_path = pd.read_csv('./stage2_paths/pathways2.tsv', sep='\t',
                          dtype=str)
    l1_product = pd.read_csv(
        './emerge_filtered_output/EMERGE_170222_product.tsv',
        sep='\t')
    l1_dist = pd.read_csv(
        './emerge_filtered_output/EMERGE_170222_distillate.tsv',
        sep='\t', index_col='gene_id')
    l2_product = l2_path.apply(
                    lambda x:
                    l1_product.apply(lambda y:
                                     evaluat_at_percent(x['reaction'], x['pathway'], y,
                                                        0.7, x['signature_definition'], l1_dist),
                                    axis=1),
                        axis=1)
    l2_product.columns = l1_product['genome']
    l2_product.index = (l2_path['pathway'] + " - " + l2_path['subpathway'])
    l2_product.T.to_csv("./stage2_paths/EMERGE_030222_product2.tsv", sep='\t')

