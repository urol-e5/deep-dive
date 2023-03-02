`deep-dive/A-Pver/data`

- [`Pver-karyotype.tab`](https://github.com/urol-e5/deep-dive/blob/main/A-Pver/data/Pver-karytotype.tab): Tab-delimited file containing scaffold name and scaffold length. Derived from NCBI _P.verrucosa_ genome GCA_014529365.1. [View notebook entry for details.](https://robertslab.github.io/sams-notebook/2023/02/15/Data-Wrangling-Create-P.verrucosa-GCA_014529365.1-Karyotype-File.html).


- [`metadata.RNAseq.csv`](https://github.com/urol-e5/deep-dive/blob/main/A-Pver/data/metadata.RNAseq.csv): Comma-separated file containing metadata for _P.verrucosa_ RNA-seq. Taken from [Danielle Becker's original repo](https://github.com/hputnam/Becker_E5/blob/master/RAnalysis/Data/RNA-seq/metadata.RNAseq.csv) on 20230224 by Sam White. Contains the following columns with headers:

  - `fragment.ID`:

  - `treatment`: `control` or `enriched`.

  - `block`:

  - `sample_id`: Unique sample id. `control` samples start with `C`. `enriched` samples start with `E`.