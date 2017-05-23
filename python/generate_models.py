#!/usr/bin/python

import argparse
import os

parser = argparse.ArgumentParser(
    description='Generate stan models for sensitivity analysis.')
parser.add_argument('--model_name', help='Directory and base filename for stan scripts.')


class StanBlocks(object):
    def __init__(self, model_name):
        # TODO: add tabs if they're not there.
        f = open(model_name + '_data_block.stanblock', 'r')
        self.data_block = f.read()
        f.close()

        f = open(model_name + '_model_block.stanblock', 'r')
        self.model_block = f.read()
        f.close()

        f = open(model_name + '_parameters_block.stanblock', 'r')
        self.parameters_block = f.read()
        f.close()

        f = open(model_name + '_hyperparameters_block.stanblock', 'r')
        self.hyperparameters_block = f.read()
        f.close()

        extra_blocks_fname = model_name + '_extra_blocks.stanblock'
        if os.path.isfile(extra_blocks_fname):
            f = open(extra_blocks_fname, 'r')
            self.extra_blocks = f.read()
            f.close()
        else:
            self.extra_blocks = ''


def write_base_script(model_name, blocks):
    f = open(os.path.join(model_name + '_generated.stan'), 'w')

    f.write('data {\n')
    f.write(blocks.data_block)
    f.write(blocks.hyperparameters_block)
    f.write('\n}\n')

    f.write('parameters {\n')
    f.write(blocks.parameters_block)
    f.write('\n}\n')

    f.write('model {\n')
    f.write(blocks.model_block)
    f.write('\n}\n')

    f.write(blocks.extra_blocks)
    f.write('\n')

    f.close()



def write_sensitivity_script(model_name, blocks):
    f = open(os.path.join(model_name + '_sensitivity.stan'), 'w')

    f.write('data {\n')
    f.write(blocks.data_block)
    f.write('\n}\n')

    f.write('parameters {\n')
    f.write(blocks.parameters_block)
    # TODO: make sure to strip constraints out of the hyperparameter block
    # in the sensitivity script.
    f.write(blocks.hyperparameters_block)
    f.write('\n}\n')

    f.write('model {\n')
    f.write(blocks.model_block)
    f.write('\n}\n')

    f.write(blocks.extra_blocks)
    f.write('\n')

    f.close()


if __name__ == "__main__":
    args = parser.parse_args()
    blocks = StanBlocks(args.model_name)
    write_base_script(args.model_name, blocks)
    write_sensitivity_script(args.model_name, blocks)
