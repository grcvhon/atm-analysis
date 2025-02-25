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
First, we refer to the file: `DARTseq_master.xlsx` (version as of 11 February 2025; file not stored in this repo). We then filter, in MS Excel, for <i>A. apraefrontalis</i> and <i>A. foliosquama</i>. Take note of some of the comments as some samples may be contaminated or other outright low quality. Nonetheless, we take all rows that are either <i>A. apraefrontalis</i> or <i>A. foliosquama</i>.<br>

This step has been done manually and output is shown below (first 10 entries):

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

<i>NB: For complete output, see: </i>`atm_genetic_dataset.csv`;<i> file not stored in this repo.</i>

#### 2) Use `awk` to filter and arrange subset data
In bash, we use the following commands to further manipulate our subset data:
* `awk -F, '{ if ($4 ~ /apraefrontalis/ && $11 ~ /yes/) { print $2, $3, $4, $9, $11 } }' atm_genetic_dataset.csv`
* `awk -F, '{ if ($4 ~ /foliosquama/ && $11 ~ /yes/) { print $2, $3, $4, $9, $11 } }' atm_genetic_dataset.csv`
These commands will print out rows under columns 2, 3, 4, 9, and 11 (SampleID, Genus, Species, FASTQ.gz prefix, Use) if column 4 (Species) contains apraefrontalis/foliosquama and if column 11 (Use) says "yes" as in good quality/usable sample.

