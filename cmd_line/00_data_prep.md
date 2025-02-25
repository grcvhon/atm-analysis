### Determine samples of <i>Aipysurus apraefrontalis</i> and <i>A. foliosquama</i> with RADseq data

This step involves collating information regarding samples of <i>Aipysurus apraefrontalis</i> and <i>A. foliosquama</i> that have available RADseq data (i.e., available in Phoenix HPC:/uofaresstor/sanders_lab/).<br>

The goal of this collation is to generate a sample sheet that has at least the following information: `order`,`dart_id`,`id_clean` (for an example, see https://github.com/a-lud/sea-snake-dart/blob/main/data/sample-sheets/240524-sample-linkage.csv).<br>

* `order` corresponds to the DaRT order number (`DNote##-####`)
* `dart_id` corresponds to the `.FASTQ.gz` prefix
* `id_clean` for example, <i>Hydrophis major</i> with KLS#1010 and FASTQ prefix 1234567: `HMA-KLS1010-1234567`

This format will improve efficiency when processing samples prior to any analyses and when using the workflow in genomics analyses.<br>

<i>NB: Scripts were not used entirely to collate information and some manual manipulation is required (e.g., via MS Excel).</i><br>

---

#### 1) Extract information from DaRTseq master spreadsheet
First, we refer to the file: `DARTseq_master.xlsx` (version as of 11 February 2025)