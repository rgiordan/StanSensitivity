#!/usr/bin/python

from generate_models import *
import string
import unittest

data_block = """
data {
  int<lower=0> y;
}
"""

hp_block = """
hyperparameters {
  real cauchy_par;
}
"""

par_block = """
parameters {
  real<lower=0> alpha;
}
"""

model_block = """
model {
  alpha ~ cauchy(cauchy_par, 10);
  target += neg_binomial_lpmf(y | alpha, 1.0);
}
"""

extra_block = """
extra {
  real beta;
}
"""

all_blocks = [ data_block, hp_block, par_block, model_block, extra_block ]


def join_blocks(contents):
    script = []
    for k, v in contents.iteritems():
        script.append('{} {{'.format(k))
        script.append(v)
        script.append('}')

    script.append('\n')
    return ''.join(script)


def rm_ws(str):
    return str.translate(None, string.whitespace)


class TestGenerators(unittest.TestCase):
    def test_parse_model_blocks(self):
        model = '\n'.join(all_blocks)
        contents = parse_model_blocks(model)
        check_contents(contents) # Should not raise an error
        self.assertEqual(rm_ws(join_blocks(contents)), rm_ws(model))

        # Make sure errors are raised if blocks are missing.
        contents = parse_model_blocks(
            '\n'.join([ data_block, par_block, model_block, extra_block ]))
        with self.assertRaises(ValueError):
            check_contents(contents)

        contents = parse_model_blocks(
            '\n'.join([ data_block, hp_block, model_block, extra_block ]))
        with self.assertRaises(ValueError):
            check_contents(contents)

        contents = parse_model_blocks(
            '\n'.join([ hp_block, par_block, model_block, extra_block ]))
        with self.assertRaises(ValueError):
            check_contents(contents)

    def test_generate_scripts(self):
        model = '\n'.join(all_blocks)
        contents = parse_model_blocks(model)
        gen_model = generate_base_script(contents)
        sens_model = generate_sensitivity_script(contents)

        gen_contents = parse_model_blocks(gen_model)
        sens_contents = parse_model_blocks(sens_model)

        expected_blocks = set([ 'data', 'parameters', 'model', 'extra' ])
        self.assertTrue(set(gen_contents.keys()) == expected_blocks)
        self.assertTrue(set(sens_contents.keys()) == expected_blocks)

        self.assertTrue('real cauchy_par;' in gen_contents['data'])
        self.assertTrue('real cauchy_par;' in sens_contents['parameters'])

        self.assertFalse('real cauchy_par;' in sens_contents['data'])
        self.assertFalse('real cauchy_par;' in gen_contents['parameters'])


if __name__ == '__main__':
    unittest.main()
