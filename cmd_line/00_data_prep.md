### Determine samples of <i>Aipysurus apraefrontalis</i> and <i>A. foliosquama</i> with RADseq data

This step involves collating information regarding samples of <i>Aipysurus apraefrontalis</i> and <i>A. foliosquama</i> that have available RADseq data (i.e., stored in `PhoenixHPC:/uofaresstor/sanders_lab/`).<br>

The goal of this collation is to generate a sample sheet that has at least the following information: `order`,`dart_id`,`id_clean` (for an example, see https://github.com/a-lud/sea-snake-dart/blob/main/data/sample-sheets/240524-sample-linkage.csv).<br>

* `order` corresponds to the DaRT order number (`DNote##-####`)
* `dart_id` corresponds to the `.FASTQ.gz` prefix
* `id_clean` for example, <i>Hydrophis major</i> with KLS 1010 and FASTQ prefix 1234567: `HMA-KLS1010-1234567` (no whitespaces)
* `barcode9l`
* `barcode`

This format will improve efficiency when processing samples prior to any analyses and when using the workflow in genomics analyses.<br>

<i>NB: Scripts were not used entirely to collate information and some manual manipulation is required (e.g., via MS Excel).</i><br>

---

#### 1) Extract information from DaRTseq master spreadsheet
First, we refer to the file: `DARTseq_master.xlsx` (version as of 11 February 2025; file not stored in this repo). We then filter, in MS Excel, for <i>A. apraefrontalis</i> and <i>A. foliosquama</i>. Take note of some of the comments as some samples may be contaminated or other outright low quality. Nonetheless, we take all rows that are either <i>A. apraefrontalis</i> or <i>A. foliosquama</i>.<br>

We also added columns to contain information on latitude and longitude (obtained, if present, from `The_One_Spreadsheet` and other field data sheets; files not stored in this repo), and if sample is usable (yes/no) based on information from `DARTseq_master.xlsx`.<br>

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

<br>

#### 2) Use command line to generate our sample sheet file (`atm_genetic_dataset.csv`)
In bash, we use the following commands to further manipulate our subset data:
* `awk -F, '{ if ($4 ~ /apraefrontalis/ && $11 ~ /yes/) { print $2, $3, $4, $9, $11 } }' atm_genetic_dataset.csv`
* `awk -F, '{ if ($4 ~ /foliosquama/ && $11 ~ /yes/) { print $2, $3, $4, $9, $11 } }' atm_genetic_dataset.csv`

These commands will print out row information for these columns: `SampleID`, `Genus`, `Species`, `FASTQ.gz prefix`, `Use`; <b>IF</b>: 
* column 4 (`Species`) contains `apraefrontalis`/`foliosquama`, and 
* column 11 (`Use`) says "yes" i.e., good quality/usable sample.
<br>

Doing this for <i>A. apraefrontalis</i>:

```bash
awk -F, '{ if ($4 ~ /apraefrontalis/ && $11 ~ /yes/) { print $2, $3, $4, $9, $11 } }' atm_genetic_dataset.csv
```

```
# output
Aaprae 4.12.01 Aipysurus apraefrontalis 2562202 yes
KLS0834 Aipysurus apraefrontalis 2562130 yes
SS171013-03 Aipysurus apraefrontalis 2562139 yes
KLS1484 Aipysurus apraefrontalis 3517861 yes
KLS1486 Aipysurus apraefrontalis 3517868 yes
KLS1490 Aipysurus apraefrontalis 3517879 yes
KLS1435 Aipysurus apraefrontalis 3593375 yes
KLS1436 Aipysurus apraefrontalis 3593362 yes
KLS1454 Aipysurus apraefrontalis 3593372 yes
KLS1457 Aipysurus apraefrontalis 3593394 yes
KLS1459 Aipysurus apraefrontalis 3593395 yes
KLS1465 Aipysurus apraefrontalis 3593393 yes
KLS1468 Aipysurus apraefrontalis 3593397 yes
KLS1477 Aipysurus apraefrontalis 3593356 yes
KLS1509 Aipysurus apraefrontalis 3593337 yes
```
<br>

Knowing that the command takes the samples of <i>A. apraefrontalis</i> with RADseq data (and their FASTQ.gz prefix) that we want to use (i.e., yes), we can expand the command to produce our sample sheet file.

```bash
echo "order","dart_id","id_clean" > aap-sample-sheet.csv
awk -F, '{ if ($4 ~ /apraefrontalis/ && $11 ~ /yes/) { gsub(/ /,"_"); print $8"," $9","toupper(substr($3,1,1))toupper(substr($4,1,2))"-"$2"-"$9 } }' atm_genetic_dataset.csv >> aap-sample-sheet.csv
```

Preview our `aap-sample-sheet.csv`:
|order       |dart_id|id_clean                  |
|------------|-------|--------------------------|
|DNote21-6332|2562202|AAP-Aaprae_4.12.01-2562202|
|DNote21-6332|2562130|AAP-KLS0834-2562130       |
|DNote21-6332|2562139|AAP-SS171013-03-2562139   |
|DNote23-8556|3517861|AAP-KLS1484-3517861       |
|DNote23-8556|3517868|AAP-KLS1486-3517868       |
|DNote23-8556|3517879|AAP-KLS1490-3517879       |
|Dnote23-8773|3593375|AAP-KLS1435-3593375       |
|Dnote23-8773|3593362|AAP-KLS1436-3593362       |
|Dnote23-8773|3593372|AAP-KLS1454-3593372       |
|Dnote23-8773|3593394|AAP-KLS1457-3593394       |
|Dnote23-8773|3593395|AAP-KLS1459-3593395       |
|Dnote23-8773|3593393|AAP-KLS1465-3593393       |
|Dnote23-8773|3593397|AAP-KLS1468-3593397       |
|Dnote23-8773|3593356|AAP-KLS1477-3593356       |
|Dnote23-8773|3593337|AAP-KLS1509-3593337       |
