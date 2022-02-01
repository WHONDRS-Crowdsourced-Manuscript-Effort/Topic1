#-- Example folder for analysis folder

Create an analysis folder with the syntax "(number in sequence of analysis)_descriptive.word"
E.g. "1_data.cleaning", "2_mulvar.comparison"

For each analysis include:
- Scripts
- Input datasets

|   |Dataset                 |Replicate merging           |Used rarity cutoff         | Original number of mf/mz| Number of mf/mz kept| Percentage retained|
|:--|:-----------------------|:---------------------------|:--------------------------|------------------------:|--------------------:|-------------------:|
|1  |Molecular formulae (mf) |MF in at least 1 replicate  |None                       |                    37528|                37528|               100.0|
|2  |Molecular formulae (mf) |MF in at least 1 replicate  |More than 1 sample (rar1)  |                    37528|                23317|                62.1|
|3  |Molecular formulae (mf) |MF in at least 1 replicate  |More than 2 samples (rar2) |                    37528|                18738|                49.9|
|21 |Molecular formulae (mf) |MF in at least 2 replicates |None                       |                    37528|                15089|                40.2|
|31 |Molecular formulae (mf) |MF in at least 2 replicates |More than 1 sample (rar1)  |                    37528|                10130|                27.0|
|4  |Molecular formulae (mf) |MF in at least 2 replicates |More than 2 samples (rar2) |                    37528|                 8663|                23.1|
|32 |Peaks (mz)              |MZ in at least 1 replicate  |None                       |                    91338|                91338|               100.0|
|5  |Peaks (mz)              |MZ in at least 1 replicate  |More than 1 sample (rar1)  |                    91338|                45583|                49.9|
|6  |Peaks (mz)              |MZ in at least 1 replicate  |More than 2 samples (rar2) |                    91338|                35160|                38.5|
|41 |Peaks (mz)              |MZ in at least 2 replicates |None                       |                    91338|                28528|                31.2|
|7  |Peaks (mz)              |MZ in at least 2 replicates |More than 1 sample (rar1)  |                    91338|                16125|                17.7|
|8  |Peaks (mz)              |MZ in at least 2 replicates |More than 2 samples (rar2) |                    91338|                13253|                14.5|
