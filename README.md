# Involution braid solver


## Setup

```
pip3 install -r requirements/base.txt
```

## Run

```
python3.5 solve.py --help
```

```
python3.5 solve.py --type <Coxeter graph type> \
				   --rank <Coxeter graph rank> \
				   [--verify] \
				   [--verbosity <verbosity level>] \
				   [--limit <maximum braid length>]
```

For example:
```
python3.5 solve.py --type D --rank 4 --verify
```
```
python3.5 solve.py --type 2F --rank 4
```
```
python3.5 solve.py --type A~ --rank 4 --verify --limit 6
```

## Test

```
pip3 install -r requirements/test.txt
tox
```
Review coverage by opening `htmlcov/index.html`.
