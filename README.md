# BAMBAM - a rare k-mer read filtering tool
A tool for filtering BAM files based on the presence, position and structural compliance of rare k-mers.

## How to install
If you have not already, install Rust using this link:  
https://rust-lang.org/tools/install/

After cloning this repo, just run
```terminaloutput
$ cargo build --release
```
You will find the compiled binary inside ./target/release  

## Usage
```terminaloutput
$ bambam <INPUT_BAM_FILE> <OUTPUT_BAM_FILE> [options]
```

### Arguments
|Argument|Description|
|---|---|
|<INPUT_BAM_FILE>|The path to the input BAM file to be filtered.|
|<REF_FILE>|The file containing the reference genome (FASTA).|

### Options
| Flag | Long name | Description | Default |
|---|---|---|---|
|-o| --output <OUTPUT>|        The path where the filtered output BAM file will be written  | ./out.bam |
|| --min <MIN>|              Minimal k-mer count in the reference to be considered "rare" | 5 |  
|| --max <MAX> |            Maximal k-mer count in the reference to be considered "rare" | 10 |  
  |-i| --inf <MIN_COUNT> |       Minimal count of k-mers in a read for the read to be kept | 1 |
  |-p| --pct <PCT>  |            Minimum percentage of intact expected rare k-mers required to keep a read | 100 |
  |-l| --len <LEN>  |            Length of the unique k-mer (max: 32) | 31
  |-t| --tolerance <TOLERANCE> | Maximum allowed edit distance (or base tolerance if --dyn-tol is used) | 0 |
  |-d| --dyn-tol |                Enable dynamic, density-aware tolerance based on k-mer isolation | |
  |-b| --bed <BED>  |            Optional: Path to output a BED file of the rare k-mer coordinates | |
  |-h| --help    |               Print help  |
  |-V| --version  |              Print version|

