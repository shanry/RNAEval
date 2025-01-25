# RNAEval
Evaluate RNA sequences across various metrics

## Dependency
```~/biology/LinearPartition/linearpartition```

## Online mode
```
echo "GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGAGGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA" | python ensemble.py -b 0 -d 0
```

## Process fasta
```
python ensemble.py -f data/test.fasta -b 0 -d 0
```

## Plot entropy
```
python plot.py -f data/test.fasta.csv -o entropy_test
```

## Plot base pair probabilities

### read sequences from standard input
```
echo -e "AUGAUCGAGGUGGUCUGCAAUGACCGCCUCGGUAAGAAGGUGAGGGUGAAGUG\nAUGAGCGAGGUCGUCUGCUAUGACCGCCUCGGUAAGAAGGUCCGGGUAAAGUG" | python bpp_plot.py -l
```

### read sequences by IDs in data/design_results.fasta
```python bpp_plot.py --name Q9BZL1```