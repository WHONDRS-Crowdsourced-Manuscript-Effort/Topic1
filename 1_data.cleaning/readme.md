#-- Example folder for analysis folder

Create an analysis folder with the syntax "(number in sequence of analysis)_descriptive.word"
E.g. "1_data.cleaning", "2_mulvar.comparison"

For each analysis include:
- Scripts
- Input datasets

|Dataset                 |Replicate merging           |Used rarity cutoff         | Original number of mf/mz| Number of mf/mz kept| Percentage retained|
|:-----------------------|:---------------------------|:--------------------------|------------------------:|--------------------:|-------------------:|
|Molecular formulae (mf) |MF in at least 1 replicate  |More than 1 sample (rar1)  |                    37528|                23317|            62.13227|
|Molecular formulae (mf) |MF in at least 1 replicate  |More than 2 samples (rar2) |                    37528|                18738|            49.93072|
|Molecular formulae (mf) |MF in at least 2 replicates |More than 1 sample (rar1)  |                    37528|                10130|            26.99318|
|Molecular formulae (mf) |MF in at least 2 replicates |More than 2 samples (rar2) |                    37528|                 8663|            23.08410|
|Peaks (mz)              |MZ in at least 1 replicate  |More than 1 sample (rar1)  |                    91338|                45583|            49.90584|
|Peaks (mz)              |MZ in at least 1 replicate  |More than 2 samples (rar2) |                    91338|                35160|            38.49438|
|Peaks (mz)              |MZ in at least 2 replicates |More than 1 sample (rar1)  |                    91338|                16125|            17.65421|
|Peaks (mz)              |MZ in at least 2 replicates |More than 2 samples (rar2) |                    91338|                13253|            14.50984|
