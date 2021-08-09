# Dependencies

Most dependencies for the pipeline are met by the two conda_requirements
```conda_protmap.yml``` and ```conda_protmap_R.yml```.

The dependencies are tested on a ubuntu docker enviroment. To setup the docker
enviroment the script ```./setup_docker.sh``` is used. All further depencencies
are shown there. To meet the other depencies that are not covered by docker
adapte this script accordingly.

# Transcriptom

In the contribution we used all raw reads from the bio sample ```PRJNA655119```
for the expression analysis.
```https://www.ncbi.nlm.nih.gov/bioproject/?term=prjna655119```


# Genomes

For the publication the genomes for the following bacteria are downloaded.

This is handled by the script ```./build_db/download_all_genomes.sh```.

The short names are used throughout the scripts and should not be changed.
The full names are the following.

| short name | scientific name                                                         | taxid   | assembly id     |
| -          | -                                                                       | -       | -               |
| anaero     | Anaerostipes caccae DSM 14662                                           | 411490  | GCA_014131675.1 |
| bact       | Bacteroides thetaiotaomicron VPI5482                                    | 226186  | GCA_014131755.1 |
| bifi       | Bifidobacterium longum NCC2705                                          | 206672  | GCF_000007525.1 |
| blautia    | Blautia producta ATCC 27340 DSM 2950                                    | 1121114 | GCA_014131715.1 |
| clostri    | Clostridium butyricum DSM 10702                                         | 1316931 | GCA_014131795.1 |
| ecoli      | Escherichia coli str K12 substr MG1655                                  | 511145  | GCF_000005845.2 |
| ery        | Erysipelatoclostridium ramosum DSM 1402                                 | 445974  | GCA_014131695.1 |
| lacto      | Lactobacillus plantarum subsp plantarum ATCC 14917 JCM 1149 CGMCC 12437 | 525338  | GCA_014131735.1 |

# Structure of scripts

All major steps are covered in ther own sub directory.

 1. build_db:
   * Builds databases for comet, downloads genomes
 4. comet:
   * MS data is gathered, PSMs are generated.
 5. transcritptom
   * All scripts for the mapping of the transcriptomic reads are here.
 6. data_accumulation
   * Here most of the final analysis are done
 7. candidates
   * Here the candidate selection and evaluation is done.
 8. UCSC_track_tools
   * The UCSC track hub is generated here
 9. figure_plotting
   * All scripts for plots that where automatically generated from data are here.
 8. start_anno_html
   * The result for evidence of early annotation startsites are generated here.

# Parameters

The file parameters.json holds paramters for the script to run and must be changed for each system.

 1. session_id
   * trackhub session id
 2. hub_id
   * hub_id for UCSC genome browser
 3. data_dir
   * dir to store all data. 1.5 TB at least.
 4. tmp_dir
   * dir for temporary files
 5. publication_dir
   * dir to output figures and infos
 6. ms_dir
   * dir of ms data. Not needed if PRIDE is available.
 7. chrome_bin
   * dir to CORRECT (version) chrome binary .
 8. bin_path
   * dir where comet etc is expected.
 9. blastdb_dir
   * blast nt db.
