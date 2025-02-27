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
First, we refer to the file: `DARTseq_master.xlsx` (version as of 11 February 2025; file not stored in this repo). We then filter, in MS Excel, for <i>Aipysurus apraefrontalis</i> and <i>Aipysurus foliosquama</i>. Take note of some of the comments as some samples may have been contaminated or of just low quality. Nonetheless, we take all rows that are either <i>A. apraefrontalis</i> or <i>A. foliosquama</i>.<br>

We also added columns to contain information on latitude and longitude (if present, obtained from `The_One_Spreadsheet` and other field data sheets; files not stored in this repo), and if sample is usable (yes/no) based on information from `DARTseq_master.xlsx`.<br>

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

#### 2) Use command line to initialise our sample sheet file
From our `atm_genetic_dataset.csv` file, we want to initialise the first 3 columns of our sample sheet in the desired format. Using the following command, let us extract the samples that have a "yes" (i.e., usable) in the `Use` column of our `atm_genetic_dataset.csv`.
<br>
```bash
awk -F, '{ if ( $11 ~ /yes/ ) { print $2, $3, $4, $9, $11 } }' atm_genetic_dataset.csv
```
This command goes: if column 11 (`Use`) is "yes", print out information for these columns: `SampleID`, `Genus`, `Species`, `FASTQ.gz`, `Use`<br>
```
# output
Aaprae 4.12.01 Aipysurus apraefrontalis 2562202 yes
KLS0834 Aipysurus apraefrontalis 2562130 yes
SS171013-03 Aipysurus apraefrontalis 2562139 yes
Afo1 Aipysurus foliosquama 2562140 yes
Afo8 Aipysurus foliosquama 2562249 yes
Afo8 Aipysurus foliosquama 2571080 yes
SS171014-02 Aipysurus foliosquama 2562167 yes
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
KLS1202 Aipysurus foliosquama 3593377 yes
KLS1696 Aipysurus foliosquama 4013436 yes
KLS1700 Aipysurus foliosquama 4013440 yes
KLS1701 Aipysurus foliosquama 4013441 yes
KLS1702 Aipysurus foliosquama 4013442 yes
KLS1707 Aipysurus foliosquama 4013447 yes
KLS1708 Aipysurus foliosquama 4013448 yes
KLS1710 Aipysurus foliosquama 4013450 yes
```
Knowing that the command takes the samples we want, we can expand the command to produce the first 3 columns of our sample sheet file in the desired format.

```bash
# generate headers
echo "order","dart_id","id_clean" > sample-sheet.csv

# append output of command below on to the `sample-sheet.csv` file
awk -F, '{ if ( $11 ~ /yes/ ) { gsub(/ /,"_"); print $8"," $9","toupper(substr($3,1,1))toupper(substr($4,1,2))"-"$2"-"$9 } }' atm_genetic_dataset.csv >> sample-sheet.csv
```

Preview our `sample-sheet.csv`:
|order       |dart_id|id_clean                  |
|------------|-------|--------------------------|
|DNote21-6332|2562202|AAP-Aaprae_4.12.01-2562202|
|DNote21-6332|2562130|AAP-KLS0834-2562130       |
|DNote21-6332|2562139|AAP-SS171013-03-2562139   |
|DNote21-6332|2562140|AFO-Afo1-2562140          |
|DNote21-6332|2562249|AFO-Afo8-2562249          |
|DNote21-6332|2571080|AFO-Afo8-2571080          |
|DNote21-6332|2562167|AFO-SS171014-02-2562167   |
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
|Dnote23-8773|3593377|AFO-KLS1202-3593377       |
|DNote24-9763|4013436|AFO-KLS1696-4013436       |
|DNote24-9763|4013440|AFO-KLS1700-4013440       |
|DNote24-9763|4013441|AFO-KLS1701-4013441       |
|DNote24-9763|4013442|AFO-KLS1702-4013442       |
|DNote24-9763|4013447|AFO-KLS1707-4013447       |
|DNote24-9763|4013448|AFO-KLS1708-4013448       |
|DNote24-9763|4013450|AFO-KLS1710-4013450       |

We will add the `barcode9l` and `barcode` columns in the next steps.<br>
<br>

#### 3) Extract `barcode9l`,`barcode` information from DaRTseq targets file
Each DaRTseq order comes with a `targets_*.csv` file which includes information about the specific run. This file contains `barcode9l` and `barcode` columns which are adapter sequences/cut sites we need to determine to get our raw sequences into usable form. We need to extract sample-specific information from these columns and then add them to our species 

for i in 2562202; do awk -F, '$1 ==col1 {print $2}' col1="$i" targets_HLCFMDRXY_1.csv; done
