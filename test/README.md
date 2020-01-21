README
======

Generate test/reference output files

```bash

codonw -nomenu -silent -nowarn -machine -all_indices -aro -hyd -dinuc input.fna ref/input.out ref/input.dinuc.blk
codonw -nomenu -silent -nowarn -machine -aau   input.fna /dev/null ref/input.aau.blk
codonw -nomenu -silent -nowarn -machine -raau  input.fna /dev/null ref/input.raau.blk
codonw -nomenu -silent -nowarn -machine -cu    input.fna /dev/null ref/input.cu.blk
codonw -nomenu -silent -nowarn -machine -cutab input.fna /dev/null ref/input.cutab.blk
codonw -nomenu -silent -nowarn -machine -cutot input.fna /dev/null ref/input.cutot.blk
codonw -nomenu -silent -nowarn -machine -rscu  input.fna /dev/null ref/input.rscu.blk
codonw -nomenu -silent -nowarn -machine -base  input.fna /dev/null ref/input.base.blk

```

