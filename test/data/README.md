## Regression test datasets

By using genes in the *E. coli* K-12 proteome that were identified to have
their C-terminus in the cytoplasm yields a good set to calculate tAI for.

```bash
Rscript tAI_script.R test/Daley_gfp.fna > test/Daley_gfp.fna.tAI
Rscript tAI_script.R test/ecolik12.ffn > test/ecolik12.ffn.tAI
```
