# tinymasker -- ultrasensitive alignment for masking repeats on assembly contigs

## Features


## Requirements

* CPU: x86\_64 (Intel or AMD)
* OS: Linux or macOS
* Compiler: one of gcc, clang, or icc

## Usage

Installation:

```
$ make CC=gcc -j4
$ make install PREFIX=$HOME/local
```

Repeat masking:

```
$ tinymasker -t4 -d repDB.tmi repDB.fa
$ tinymasker -t4 repDB.tmi contigs.fa > masks.gff3
$ tinymasker -t4 -g masks.gff3 contigs.fa > masked.fa    # unimplemented for now
```

Apply custom score profiles on building reference index:

```
$ tinymasker -t4 -c config.toml -d repDB.tmi repDB.fa
```

Example of custom profile:

```TOML
[[a]]
name = ""
comment = ""
match = 2
mismatch = 3
gap_open = 5
gap_extend = 1

[[b]]
name = ""
comment = ""
score_matrix = [
	[],
	[],
	[],
	[]
]
```



## Performance



## Algorithm

### Seeding with unidirected de Bruijn graph

### Filtering seeds

### Chaining seeds

### Filtering chains

### Finishing alignment



## License

Copyright (2019) Hajime Suzuki

MIT