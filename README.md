# combineSomGerm
This program reports the overlapped results between different germline and somatic results

Syntax: ./combineSomGerm [options]

Options:
   To specify the input files
      -g [main germline result file]
      -r [other germline result file]
      -s [combine somatic result file]

   To specify the output files
      -o [output vcf file for the overlapping and consistent results]
      -e1 [output log file with overlapping but inconsistent results]
      -e2 [output log file with non-overlapping records]

   To apply the 'PASS' filter on the files
      -f [file contains the list of files to apply the filter]

   Other options:
      -t : Only keep the records with normal cells which are 0/0
      -c : Only keep the consistent records

Remark: 1. If there are more than one germline result files (except the main one)
           or more than one somatic result files, please use the options '-r' or '-s'
           more than once.
        2. Options '-g', '-o', and '-e' are compulsory.

Example:
   $./combineSomGerm -g patient16011.gatk.germline.vcf -r patient16011.Monovar.germline.vcf -s patient16011.gatk.somatic.vcf -s patient16011.somaticsniper.somatic.vcf -s patient16011.varscan.somatic.vcf -o patient16011.overlap.vcf -e1 patient16011.inconsist.log -e2 patient16011.err.log -f applyPassFilterFiles.txt
