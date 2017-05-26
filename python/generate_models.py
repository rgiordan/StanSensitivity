#!/usr/bin/python

import argparse
from collections import OrderedDict
import os
import re

parser = argparse.ArgumentParser(
    description='Generate stan models for sensitivity analysis.')
parser.add_argument('--base_model',
                    help='Directory and base filename for stan scripts.')

def parse_model_blocks(script):
    # Parse the script into blocks.
    contents = OrderedDict()
    last_i = 0
    in_block = False
    for i in range(len(script)):
        if in_block:
            if script[i] == '{':
                num_parens += 1
            elif script[i] == '}':
                num_parens -= 1
            if num_parens == 0:
                # ...then we have closed the opening block bracket.
                block = script[last_i:i]
                last_i = i + 1
                in_block = False
                contents[tag] = block
        else:
            if script[i] == '{':
                # ...then we have found the opening bracked of a stan block.
                in_block = True
                tag = script[last_i:i].strip()
                last_i = i + 1
                num_parens = 1

    return contents


def generate_base_script(contents):
    script = []
    for k, v in contents.iteritems():
        if k == 'hyperparameters':
            continue
        if k == 'data':
            script.append('data {')
            script.append(contents['data'])
            script.append('\n  // Hyperparameters:')
            script.append(contents['hyperparameters'])
            script.append('}\n')
        else:
            script.append('{} {{'.format(k))
            script.append(v)
            script.append('}\n')

    script.append('\n')
    return ''.join(script)


def generate_sensitivity_script(contents):
    script = []
    for k, v in contents.iteritems():
        if k == 'hyperparameters':
            continue
        if k == 'parameters':
            script.append('parameters {')
            script.append(contents['parameters'])

            # TODO: script the constraints.
            script.append('\n  // Hyperparameters:')
            script.append(contents['hyperparameters'])
            script.append('}\n')

        else:
            script.append('{} {{'.format(k))
            script.append(v)
            script.append('}\n')

    script.append('\n')
    return ''.join(script)


def check_contents(contents):
    for block in [ 'data', 'parameters', 'hyperparameters' ]:
        if not block in contents:
            raise ValueError('The model is missing a {} block.'.format(block))


if __name__ == "__main__":
    args = parser.parse_args()
    if args.base_model[-5:] != '.stan':
        raise ValueError('The base model file must have the suffix .stan')

    f = open(args.base_model, 'r')
    script = f.read()
    f.close()

    model_name = re.sub('\.stan$', '', args.base_model)

    contents = parse_model_blocks(script)
    check_contents(contents)

    f = open(os.path.join(model_name + '_generated.stan'), 'w')
    f.write(generate_base_script(contents))
    f.close()

    f = open(os.path.join(model_name + '_sensitivity.stan'), 'w')
    f.write(generate_sensitivity_script(contents))
    f.close()
