# Rare k-mer read filtering tool
A tool for filtering BAM files based on the presence, position and structural compliance of rare k-mers.
## Usage
```terminaloutput
$ bambam <INPUT_BAM_FILE> <OUTPUT_BAM_FILE> [options]
```

### Arguments
  <INPUT_BAM_FILE>  The path to the input BAM file to be filtered  
  <REF_FILE>        The file containing the reference genome (FASTA)

### Options
  -o, --output <OUTPUT>        The path where the filtered output BAM file will be written [default: ./out.bam]  
      --min <MIN>              Minimal k-mer count in the reference to be considered "rare" [default: 5]  
      --max <MAX>              Maximal k-mer count in the reference to be considered "rare" [default: 10]  
  -i, --inf <MIN_COUNT>        Minimal count of k-mers in a read for the read to be kept [default: 1]  
  -p, --pct <PCT>              Minimum percentage of intact expected rare k-mers required to keep a read [default: 100]  
  -l, --len <LEN>              Length of the unique k-mer [default: 31, max: 32]  
  -t, --tolerance <TOLERANCE>  Maximum allowed edit distance (or base tolerance if --dyn-tol is used) [default: 0]  
  -d, --dyn-tol                Enable dynamic, density-aware tolerance based on k-mer isolation  
  -b, --bed <BED>              Optional: Path to output a BED file of the rare k-mer coordinates  
  
  -h, --help                   Print help  
  -V, --version                Print version  
