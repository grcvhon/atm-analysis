### Determine samples of <i>Aipysurus apraefrontalis</i> and <i>A. foliosquama</i> with RADseq data

This step involves collating information regarding samples of <i>Aipysurus apraefrontalis</i> and <i>A. foliosquama</i> that have available RADseq data (i.e., stored in `PhoenixHPC:/uofaresstor/sanders_lab/`).<br>

The goal of this collation is to generate a sample sheet that has at least the following information: `order`,`dart_id`,`id_clean` (for an example, see https://github.com/a-lud/sea-snake-dart/blob/main/data/sample-sheets/240524-sample-linkage.csv).<br>

* `order` corresponds to the DaRT order number (`DNote##-####`)
* `dart_id` corresponds to the `.FASTQ.gz` prefix
* `id_clean` for example, <i>Hydrophis major</i> with KLS#1010 and FASTQ prefix 1234567: `HMA-KLS1010-1234567`

This format will improve efficiency when processing samples prior to any analyses and when using the workflow in genomics analyses.<br>

<i>NB: Scripts were not used entirely to collate information and some manual manipulation is required (e.g., via MS Excel).</i><br>

---

#### 1) Extract information from DaRTseq master spreadsheet
First, we refer to the file: `DARTseq_master.xlsx` (version as of 11 February 2025). We then filter for <i>A. apraefrontalis</i> and <i>A. foliosquama</i>. Take note of some of the comments as some samples may be contaminated or other outright low quality. Nonetheless, we take all rows that are either <i>A. apraefrontalis</i> or <i>A. foliosquama</i>.<br>

This step has been done manually and output is shown below:

|Source                  |SampleID      |Genus    |Species       |Location    |Latitude    |Longitude  |DaRT_set    |FASTQ.gz|Comments                             |Use|
|------------------------|--------------|---------|--------------|------------|------------|-----------|------------|--------|-------------------------------------|---|
|PreATM_sampling         |Aaprae 4.12.01|Aipysurus|apraefrontalis|Ashmore Reef|-12.24174549|123.04166  |DNote21-6332|2562202 |Coordinates approximate              |yes|
|PreATM_sampling         |KLS0834       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.166666  |114.2999988|DNote21-6332|2562130 |Coordinates approximate              |yes|
|PreATM_sampling         |SAM R68142    |Aipysurus|apraefrontalis|            |            |           |DNote21-6332|2571051 |Low quality DaRT                     |no |
|PreATM_sampling         |SS171013-03   |Aipysurus|apraefrontalis|Pilbara     |-19.6889305 |118.220874 |DNote21-6332|2562139 |                                     |yes|
|PreATM_sampling         |Afo1          |Aipysurus|foliosquama   |Ashmore Reef|-12.24174549|123.04166  |DNote21-6332|2562140 |Coordinates approximate              |yes|
|PreATM_sampling         |Afo8          |Aipysurus|foliosquama   |Ashmore Reef|-12.24174549|123.04166  |DNote21-6332|2562249 |Coordinates approximate              |yes|
|PreATM_sampling         |Afo8          |Aipysurus|foliosquama   |Ashmore Reef|-12.24174549|123.04166  |DNote21-6332|2571080 |Coordinates approximate              |yes|
|PreATM_sampling         |KLS1001       |Aipysurus|foliosquama   |            |            |           |DNote21-6332|2562209 |WA Coast apraefrontalis_contamination|no |
|PreATM_sampling         |KLS1001       |Aipysurus|foliosquama   |            |            |           |DNote21-6332|2584016 |WA Coast apraefrontalis_contamination|no |
|PreATM_sampling         |KLS1006       |Aipysurus|foliosquama   |            |            |           |DNote21-6332|2562240 |Failed DaRT                          |no |
|PreATM_sampling         |SS171014-02   |Aipysurus|foliosquama   |Pilbara     |-19.709453  |117.8305545|DNote21-6332|2562167 |                                     |yes|
|PreATM_sampling         |KLS1484       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.166666  |114.2999988|DNote23-8556|3517861 |Coordinates approximate              |yes|
|PreATM_sampling         |KLS1486       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.166666  |114.2999988|DNote23-8556|3517868 |Coordinates approximate              |yes|
|PreATM_sampling         |KLS1490       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.166666  |114.2999988|DNote23-8556|3517879 |Coordinates approximate              |yes|
|PreATM_sampling         |KLS1435       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.10533   |114.20917  |Dnote23-8773|3593375 |                                     |yes|
|PreATM_sampling         |KLS1436       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.120333  |114.321    |Dnote23-8773|3593362 |                                     |yes|
|PreATM_sampling         |KLS1454       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.1245    |114.134833 |Dnote23-8773|3593372 |                                     |yes|
|PreATM_sampling         |KLS1457       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.1245    |114.13483  |Dnote23-8773|3593394 |                                     |yes|
|PreATM_sampling         |KLS1459       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.134833  |114.2005   |Dnote23-8773|3593395 |                                     |yes|
|PreATM_sampling         |KLS1465       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.166666  |114.2999988|Dnote23-8773|3593393 |Coordinates approximate              |yes|
|PreATM_sampling         |KLS1468       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.090333  |114.246667 |Dnote23-8773|3593397 |                                     |yes|
|PreATM_sampling         |KLS1477       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.166666  |114.2999988|Dnote23-8773|3593356 |Coordinates approximate              |yes|
|PreATM_sampling         |KLS1509       |Aipysurus|apraefrontalis|Exmouth Gulf|-22.166666  |114.2999988|Dnote23-8773|3593337 |Coordinates approximate              |yes|
|PreATM_sampling         |KLS1202       |Aipysurus|foliosquama   |Pilbara     |-20.050212  |118.288886 |Dnote23-8773|3593377 |                                     |yes|
|DPIRD_Fisheries_June2024|KLS1696       |Aipysurus|foliosquama   |Shark Bay   |-25.24705   |113.4100415|DNote24-9763|4013436 |                                     |yes|
|DPIRD_Fisheries_June2024|KLS1700       |Aipysurus|foliosquama   |Shark Bay   |-25.623425  |113.1631085|DNote24-9763|4013440 |                                     |yes|
|DPIRD_Fisheries_June2024|KLS1701       |Aipysurus|foliosquama   |Shark Bay   |-25.623425  |113.1631085|DNote24-9763|4013441 |                                     |yes|
|DPIRD_Fisheries_June2024|KLS1702       |Aipysurus|foliosquama   |Shark Bay   |-25.5985585 |113.2126165|DNote24-9763|4013442 |                                     |yes|
|DPIRD_Fisheries_June2024|KLS1707       |Aipysurus|foliosquama   |Shark Bay   |-25.020975  |113.336883 |DNote24-9763|4013447 |                                     |yes|
|DPIRD_Fisheries_June2024|KLS1708       |Aipysurus|foliosquama   |Shark Bay   |-24.926183  |113.2711835|DNote24-9763|4013448 |                                     |yes|
|DPIRD_Fisheries_June2024|KLS1710       |Aipysurus|foliosquama   |Shark Bay   |-24.9617915 |113.19455  |DNote24-9763|4013450 |                                     |yes|

