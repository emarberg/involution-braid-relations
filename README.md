# Braid relations for involution words in Coxeter groups

The code in this repository implements the algorithms described in the paper
"Braid relations for involution words in affine Coxeter groups," available online at

https://arxiv.org/abs/1703.10437.

The main program takes a twisted Coxeter system as input, and tries to compute
a minimal set of "braid relations" which span and preserve all sets of involution words
for twisted involutions in the provided Coxeter group.

## Setup
1. Create a virtual environment: `virtualenv -p python3.5 py3`
1. Active the virtual environment: `source py3/bin/activate`
1. Install requirements: `pip3 install -r requirements/base.txt`

## Tests
1. After setting up, run the tests at the command line: `tox`
1. (Optional) Review the generated coverage by opening `tests/htmlcov/index.html`

## Run
To see the full list of command line options, run:
```
python3.5 solve.py --help
```
The main things to specify when running the program are the type and rank of the 
(twisted) Coxeter system to consider:
```
python3.5 solve.py --type <Coxeter graph type> \
                   --rank <Coxeter graph rank> \
                  [--s <optional simple generator>] \
                  [--t <optional simple generator>] \
                  [--verify] \
                  [--limit <maximum braid length>] \
                  [--log=<logging level>]
```
The last five arguments are optional.
Including the `--verify` option will make the program explicitly check that the relations
it finds actually span and preserve all sets of involution words for the group under consideration.
This option is required if a numeric value for `--limit` is provided.
When `--limit` is set, the program will only search for relations of length up to the given limit.
This is a useful option if the Coxeter group is infinite. The `--log` option controls
how much output the program logs during execution.

Some examples:
```
python3.5 solve.py --type D --rank 4 --verify
```
```
python3.5 solve.py --type 2F --rank 4 --log=debug
```
```
python3.5 solve.py --type D~ --rank 5 --verify --limit 2
```
```
python3.5 solve.py --type G~ --rank 3
```

## Computations
To generate Table 1 in the accompanying paper, run these commands:

* Type 2A9: 
```
python3.5 solve.py --type 2A --rank 9 --s 4 --t 5
```
* Type B5:
```
python3.5 solve.py --type B --rank 5 --s 4 --t 5
```
* Types D7 and D'7:
```
python3.5 solve.py --type D --rank 7 --s 4 --t 5
```
```
python3.5 solve.py --type D --rank 7 --s 5 --t 6
```
* Types 2D7 and 2D'7:
```
python3.5 solve.py --type 2D --rank 7 --s 4 --t 5
```
```
python3.5 solve.py --type 2D --rank 7 --s 5 --t 6
```
