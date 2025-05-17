import re, os
from collections import defaultdict
import pandas as pd
import urllib.request

def fetch_nmrstar_file(bmrb_id):
    path=f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_3.str"
    return urllib.request.urlretrieve(path, f"bmr{bmrb_id}_3.str")

def parse_nmr_star(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    data = {}
    saveframe_name = None
    in_saveframe, in_loop = False, False
    loop_tags, loop_data = [], []
    current_saveframe, current_tags = {},{}
    current_loops = defaultdict(list)

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('save_'):
            if in_saveframe:
                current_saveframe.update(current_tags)
                current_saveframe.update(dict(current_loops))
                if saveframe_name:
                    data[saveframe_name] = current_saveframe

            saveframe_name = line[5:].strip() or None
            current_saveframe = {}
            current_tags = {}
            current_loops = defaultdict(list)
            loop_tags = []
            loop_data = []
            in_loop = False
            in_saveframe = bool(saveframe_name)
            i += 1
            continue

        if in_saveframe:
            if line == 'loop_':
                in_loop = True
                loop_tags, loop_data = [], []
                i += 1
                while i < len(lines):
                    tag_line = lines[i].strip()
                    if tag_line.startswith('_'):
                        loop_tags.append(tag_line)
                        i += 1
                    else:
                        break
                while i < len(lines):
                    data_line = lines[i].strip()
                    if data_line == 'stop_':
                        break
                    if data_line:
                        values = re.findall(r'(?:"[^"]*"|\'[^\']*\'|[^\s]+)', data_line)
                        if len(values) == len(loop_tags):
                            base_prefix = loop_tags[0].split('.')[0] + '.'
                            row = {
                                tag[len(base_prefix):] if tag.startswith(base_prefix) else tag: val.strip('"\'')
                                for tag, val in zip(loop_tags, values)
                            }
                            loop_category = base_prefix.rstrip('.')
                            current_loops[loop_category].append(row)
                    i += 1
                in_loop = False
                i += 1
                continue

            elif line == 'stop_':
                current_saveframe.update(current_tags)
                current_saveframe.update(dict(current_loops))
                if saveframe_name:
                    data[saveframe_name] = current_saveframe

                current_saveframe, current_tags = {}, {}
                current_loops = defaultdict(list)
                loop_tags, loop_data = [], []
                in_loop, in_saveframe = False, False
                saveframe_name = None
                i += 1
                continue

            elif line.startswith('_'):
                tag = line.split()[0]
                key = tag.split('.')[-1] if '.' in tag else tag
                # Check for multiline value
                if len(line.split()) == 1 and i + 1 < len(lines) and lines[i + 1].strip() == ';':
                    i += 2
                    value_lines = []
                    while i < len(lines):
                        val_line = lines[i].rstrip('\n')
                        if val_line.strip() == ';':
                            break
                        value_lines.append(val_line)
                        i += 1
                    value = ''.join(value_lines)
                    current_tags[key] = value
                    i += 1  # move past closing ';'
                else:
                    parts = line.split(None, 1)
                    value = parts[1] if len(parts) > 1 else ''
                    current_tags[key] = value.strip('"\'')
                    i += 1
                continue
            else:
                i += 1
                continue
        else:
            i += 1
            continue

    if in_saveframe and saveframe_name:
        current_saveframe.update(current_tags)
        current_saveframe.update(dict(current_loops))
        data[saveframe_name] = current_saveframe

    return data

def convert_loop_to_dataframe(loop):
    dct = {k:[] for k in loop[0].keys()}
    for entr in loop:
        for k, v in entr.items():
            dct[k].append(v)
    return pd.DataFrame.from_records(dct)

def clean_cs_dataframe(df):
    df['Val'] = df['Val'].astype(float)
    df['Val_err'] = df['Val_err'].astype(float)
    return df[['Entity_ID', 'Seq_ID', 'Auth_seq_ID','Comp_ID',
                'Atom_ID','Atom_type','Val','Val_err', ]]

def get_sequences(parsed):
    out={}
    tags=['ID', 'Polymer_type', 'Polymer_seq_one_letter_code']

    for k in parsed.keys():
        if parsed[k]['Sf_category'] == 'entity':
            out[k] = {i: parsed[k][i] for i in tags}

    return pd.DataFrame.from_records(out)

def get_sample_info(parsed):
    out = {}
    tags=['ID', 'Mol_common_name','Entity_ID','Isotopic_labeling','Concentration_val','Concentration_val_units']
    for k in parsed.keys():
        if parsed[k]['Sf_category'] == 'sample':
            for x in parsed[k]['_Sample_component']:
                out[k] = {i: x[i] for i in tags}

    return pd.DataFrame.from_records(out)

def get_chem_shifts(parsed):
    out_dfs=[]
    for k in parsed.keys():
        if parsed[k]['Sf_category'] == 'assigned_chemical_shifts':
            cs = convert_loop_to_dataframe(parsed[k]['_Atom_chem_shift'])
            cs = clean_cs_dataframe(cs)
            try:
                cs['name'] = parsed[k]['Name']
            except:
                cs['name'] = '.'

            cs['cs_saveframe_id'] = k
            out_dfs.append(cs)
    out = pd.concat(out_dfs)
    out = out.reset_index()
    return out
