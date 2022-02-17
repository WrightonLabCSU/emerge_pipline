import pandas as pd
import pytest
import re
from datetime import date


def get_value(pathway:str, step:str, values:pd.Series):
    column = [i for i in values.index
            if re.match(f"^Complex {pathway}:.*-{step}$", i)]
    if len(column) != 1:
        raise ValueError("The column name is ambiguous for"
                         f" pathway={pathway}, step={step}")
    return values[column[0]] > 0



def evaluat(logic:str, pathway:str,  values:pd.Series):
    l = 0 # level
    for i, c in enumerate(logic):
        if l == 0:
           if c == ',':
               return evaluat(logic[:i], pathway, values) or \
                   evaluat(logic[i+1:], pathway, values)
           if c == '+':
               return evaluat(logic[:i], pathway, values) and \
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


if __name__ == '__main__':
    l2_path = pd.read_csv('./stage2_paths/pathways2.tsv', sep='\t')
    l2_path.columns
    l1_output = pd.read_csv(
        './emerge_filtered_output/EMERGE_170222_product.tsv',
        sep='\t')
    l2_output = l2_path.apply(
                    lambda x:
                    l1_output.apply(lambda y:
                                    evaluat(x['reaction'], x['pathway'], y)
                                    , axis=1)
                              , axis=1)
    l2_path.iloc[0].T
    l2_path.columns
    l2_output.columns = l1_output['genome']
    l2_output.index = (l2_path['pathway'] + " - " + l2_path['subpathway'])
    l2_output.T.to_csv("./stage2_paths/EMERGE_170222_product2.tsv", sep='\t')


def test_get_value():
    t_values = pd.Series(
         {"Complex nitrogen_redox: nitrate->nitrite -0": 0.0,
          "Complex nitrogen_redox: nitrite->ammonium -1": 0.5},
         name="20100900_E1D_1")
    assert get_value(pathway='nitrogen_redox', step='1',
                         values=t_values) == True
    assert get_value(pathway='nitrogen_redox', step='0',
                         values=t_values) == False


def test_evaluat():
    values = pd.Series({'Complex nitrogen_redox: any -1': True,
                        'Complex nitrogen_redox: any -2': False,
                        'Complex nitrogen_redox: any -3': True,
                        })
    assert evaluat('2', 'nitrogen_redox', values) == False
    assert evaluat('1', 'nitrogen_redox', values) == True
    assert evaluat('1+2', 'nitrogen_redox', values) == False
    assert evaluat('1,2', 'nitrogen_redox', values) == True
    assert evaluat('1+3', 'nitrogen_redox', values) == True
    assert evaluat('1+3+2', 'nitrogen_redox', values) == False
    assert evaluat('(1+3)+(1+2)', 'nitrogen_redox',
                   values) == False
    assert evaluat('(1+2),(1+2)', 'nitrogen_redox',
                   values) == False
    assert evaluat('(1,2)+(1,2)', 'nitrogen_redox',
                   values) == True
    assert evaluat('((1+3)+(1+2)),((1,2)+(1,2))', 'nitrogen_redox',
                   values) == True
    assert evaluat('((1+3)+(1+2))+((1+3)+(1+2)),((1,2)+(1,2))',
                   'nitrogen_redox', values) == False

