#!/usr/bin/python

import argparse
import os
from shutil import copyfile


parser = argparse.ArgumentParser(
    description='Generate stan models for sensitivity analysis.')
parser.add_argument('--model_name', help='A model name to split up.')

def parse_model_blocks(script):
    # Parse the script into blocks.
    contents = {}
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

if __name__ == "__main__":
    args = parser.parse_args()

    # Make a backup copy.
    #copyfile(args.model_name + '.stan', args.model_name + '_original.stan')
    f = open(args.model_name + '.stan', 'r')
    script = f.read()
    f.close()

    contents = parse_model_blocks(script)

    # Write the blocks individually.
    assert 'data' in contents
    f = open(args.model_name + '_data_block.stanblock', 'w')
    f.write(contents['data'])
    f.close()

    assert 'model' in contents
    f = open(args.model_name + '_model_block.stanblock', 'w')
    f.write(contents['model'])
    f.close()

    assert 'parameters' in contents
    f = open(args.model_name + '_parameters_block.stanblock', 'w')
    f.write(contents['parameters'])
    f.close()

    f = open(args.model_name + '_hyperparameters_block.stanblock', 'w')
    f.write('// Manually move hyperparameters from the data block here.\n' +
            '// Note: all hyperparameters must be unconstrained real values.\n')
    f.close()

    extra_blocks = set.difference(set(contents.keys()), set(['data', 'model', 'parameters']))
    extra_blocks_filename = args.model_name + '_extra_blocks.stanblock'
    if len(extra_blocks) > 0:
        f = open(extra_blocks_filename, 'w')
        for block in extra_blocks:
            f.write(block + '{\n')
            f.write(contents[block])
            f.write('}\n')
        f.close()
    else:
        # If present (e.g. from a past run), remove it, because now there are
        # no extra blocks.
        if os.path.isfile(extra_blocks_filename):
            os.remove(extra_blocks_filename)

    print 'Done.  Found the following model blocks: '
    print contents.keys()
