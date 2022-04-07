# Software: combineSomGerm

This program reports the overlapped results between different germline and somatic results

## Installation

The software was written in C++, and it has been tested under linux and MacOS platform. You need
to have C++ compiler installed in the machine in order to compile the source codes. The compilation
steps are shown as follows:

```
tar -zxvf combineSomGerm-2.2.tar.gz
cd combineSomGerm-2.2
make
```

Then an executable file named *combineSomGerm* will appear

## Usage

Syntax:

```
./combineSomGerm [options]
```

Options:
   
  - To specify the input files

```
-g [main germline result file]
-r [other germline result file]
-s [combine somatic result file]
```

  - To specify the output files

```
-o [output vcf file for the overlapping and consistent results]
-e1 [output log file with overlapping but inconsistent results]
-e2 [output log file with non-overlapping records]
```

  - To apply the 'PASS' filter on the files

```
-f [file contains the list of files to apply the filter]
```

  - Other options:

```
-t : Only keep the records with normal cells which are 0/0
-c : Only keep the consistent records
```

Remark: 

1. If there are more than one germline result files (except the main one) or more than one somatic result files, please use the options '-r' or '-s' more than once.

2. Options '-g', '-o', and '-e' are compulsory.

## Example

Example files are enclosed in the folder.

Command:

```
./combineSomGerm -g patient16011.gatk.germline.vcf -r patient16011.Monovar.germline.vcf -s patient16011.gatk.somatic.vcf -s patient16011.somaticsniper.somatic.vcf -s patient16011.varscan.somatic.vcf -o patient16011.overlap.vcf -e1 patient16011.inconsist.log -e2 patient16011.err.log -f applyPassFilterFiles.txt
```

## Method

The followings are the methods implemented inside combineSomGerm

1. For each file, we ignore the positions without any non-empty tumor record.

2. For the main germline file, we ignore the position if all the non-empty tumor records on that position are the same as the normal record.

3. We use the main germline file as blueprint, then compare all the records with the other germline files and the somatic files. The positions with normal cells 0/1 will not be checked among the somatic records. The following lists out the detailed logic:

```
For each position on the main germline file
	if the normal cell record is 0/0
		if this position appears in all other germline files and somatic files (and all non-empty tumor records are consistent when option -c is used),
			keep the position
		else
			ignore the position
	    	end if
	else
		if this position appears in all other germline files (and all non-empty tumor records are consistent when option -c is used),
			keep the position
		else
			ignore the position
		end if
	end if
end for
```
