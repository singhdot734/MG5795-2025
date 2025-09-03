# Class 2 
## Opening up a terminal window
1. Log into OSC (ondemand.osc.edu)
2. Using the menu bar at the top, navigate to Files > /fs/scratch/PAS3124
3. Click the **>_Open in terminal** button at the top to open a new UNIX terminal window in a new tab.
4. Optional: choose a desired theme from the top right corner.

## Basic UNIX commands
### Where am I and What is in this directory?
1. `pwd` prints the current directory
2. `ls` lists contents of the current directory

### Creating new files and directories; basic navigation
1. To create a personal directory in scratch space of the classroom project (PAS3124), use `mkdir yourname`
2. The general structure of unix commands is: program options input output. In this case, program is "mkdir" and output is "unixmore", which is being run as follows (there are no options or input file being specified here):
`mkdir yourname`
3. To enter/change directory to "yourname", use `cd yourname`
4. Now try, `pwd`
5. Let's do the customary first command on UNIX: `echo "Hello world!"`
6. This will write **"Hello world!"** to the screen.
7. Output of any program (or command) can be directed to a file using `>`
8. "Hello world!" can be written to a file (instead of screen) by changing the command on line 13 to `echo "Hello world!" > myfile_1.txt`
9. Now do `ls` again. Now it will show "myfile_1.txt" in your personal directory. This is your first file in UNIX!
10. Try `ls -l`. "-l" is an option that will tell "ls" to list items with more details.
11. To move up a directory: `cd ..`
12. To move up two directories: `cd ../..` ...and so forth.
13. If you ever type just `cd`, this will take you to your home directory on OSC. To come back to your personal directory in the classroom project, you will have to type the **absolute path** to your directory: `cd /fs/scratch/PAS3124/yourname`

### Getting data and copying/moving files
1. To get SARS-CoV2 genome as a FASTA file, we will use a command called `curl`.
2. Copy this into your terminal `curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" -o NC_045512.2.fasta`
3. In the above command option `-o` tells what the output of `curl` should be.
4. This will download SARS-CoV2 genome sequence FASTA file into the directory you were in where you ran the code on line 23. To make sure you are in the correct directory, you can do `pwd` before proceeding with line 23.
5. Confirm with `ls` to make sure the downloaded file "NC_045512.2.fasta" is there.
6. Now we will move this FASTA file into a new directory.
7. In your personal directory, use `mkdir` to make a new directory called "sarscov2"
8. Do `ls` to confirm.
9. Change into the "sarscov2" directory using `cd`
10. `cp ../NC_045512.2.fasta .` will copy "NC_045512.2.fasta" file into this directory. `cp` is the program. `../NC_045512.2.fasta` gives the relative path to the file to be copied. `.` means copy to current directory.
11. Instead of copying, you can move the the FASTA file using this command `mv ../NC_045512.2.fasta .`. Now the FASTA file will no longer remain in the directory one level up.

# Class 3
## Opening up a terminal window
1. Log into OSC (ondemand.osc.edu)
2. Using the menu bar at the top, navigate to Files > /fs/scratch/PAS3124
3. Click the **>_Open in terminal** button at the top to open a new UNIX terminal window in a new tab.
4. Optional: choose a desired theme from the top right corner.
   
## More basic UNIX commands
1. Go to the "sarscov2" directory in "yourname" directory. You will use `cd`.
2. Tip: Start typing the directory name and then press tab. It will try to complete the name until it hits a unique character. This will save you lot of typing!

### Viewing files
1. `head` or `tail` commands can be used to view contents of a file. `head NC_045512.2.fasta` will show the first 10 lines of the FASTA file. `head -n 20 NC_045512.2.fasta` uses option `-n` with specified variable 20 and will show first 20 lines.
2. `tail NC_045512.2.fasta` will show the last 10 lines of the FASTA file.` How will you look at the last 25 lines of the FASTA file?
3. More ways to view files: `more` or `less` followed by `filename` will print it to the screen in scrollable format. Scroll up or down using arrows, or down by hitting space. Press `q` to quit more and return to command prompt
4. Files can also be viewed using the `cat` command. It will print it all to the screen.
5. `cat` can be used to concatenate files. It may show up later.

### Exercise 1:
Create a new text file that contains the first 20 lines of the SARS-CoV2 fasta file. Save this file as "NC_045512.2.head20.txt". Take a screenshot of your work to submit.

### Getting help
1. `man` followed by program name will open manual for that program to  view options available for that program/command. For example, `man head` will open manual for `head`.
2. What does the option `-n` mean?
3. Scroll up or down manual pages using arrows, or down by hitting space. Press q to quit and return to command prompt.

### Counting lines and characters
1. Can we get the size of the SARS-CoV2 genome by counting the number of nucleotide characters in the FASTA file? Yes we can!
2. Make sure you are in the `sarscov2` directory. Check if the FASTA file is there by using `ls`
3. Command `wc` gives word, line and character counts for a particular file. Type `wc NC_045512.2.fasta` and press enter
4. To get only line numbers, `wc` can be run with `-l` option. `wc -c` will give only the character count.
5. What is the character count? Compare it to the genome size on Genome Browser. Why is it different?
6. Is it the first line in the FASTA file that has the name and other info about the file? That can't be it.
7. We can actually investigate this if we learn couple of more tricks, i.e., commands.

### Combining two different commands
1. **`|`** or pipe can be used to daisy-chain two or more commands. It pipes the output of the first command as input into the second.
2. Try this.
3. `tail -23 NC_045512.2.fasta | head -1` will pipe the output of tail (last 23 lines of the FASTA file) into `head -1` to show the 23rd line from the bottom of the FASTA file.
4. Always think of more ways to use the pipe.

### Finding patterns using grep
1. **`grep`** or "global regular expression print", is a search tool in UNIX. It can be used to find patterns and either select them, show them or even ignore them.
2. `grep "A" NC_045512.2.fasta` will find all A's in the FASTA file.
3. Now if we "pipe" the output of the previous line into `wc -l`, we can count how many lines have A's. Type this `grep "A" NC_045512.2.fasta | wc -l`
4. How many lines have A's in the FASTA file?
5. `grep` with `-v` option will find a pattern and hide it. So `grep -v "A" NC_045512.2.fasta` will show lines that DO NOT have A's. To count these lines, you can pipe the output of `grep` to `wc -l`
6. `grep` with `-o` option will print only the matched parts, each on a separate line. If we do `grep -o "A" NC_045512.2.fasta` and pipe it into `wc -l` we can count the number of A's in the file!
7. `grep -o "[AT]" NC_045512.2.fasta` will count all A's or T's in the FASTA file.
8. Tip: You can press up arrow on command prompt to show previous commands you have typed. Keep hitting the up (and down) arrows as needed to bring up an old command and then you can edit it as necessary. Also saves lots of typing! And errors!!

### Exercise 2: 
1. Can you now tell the size of the SARS-CoV2 genome? Take a screenshot of your work to submit.
2. Why is the number still not exactly the same as on the Genome Browser?

### Getting human chromosome 22 sequence
1. Make a new directory called `chr22` in your personal directory 
2. Go to this new directory `cd chr22`
3. Get human chr22 fasta file using `curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz > chr22.fa.gz`
4. Note how `>` is used to save the curl output to a file.

### Data decompression
1. Three ways to decompress the `chr22.fa.gz` file (any one can be used):
   
```
      gzip -d chr22.fa.gz
      gunzip chr22.fa.gz
      gunzip -c chr22.fa.gz > chr22.fa
```
2. Third option is preferred as `-c` option saves the original file while the first two options do not.
    
### Data compression (skip if no time)
1. Files and directories can be zipped into `tar` (tape archive) format.
2. To do this, let's make new files out of `chr22.fa`, move them into a directory of their own.
3. Write out first 10 lines of the `chr22.fa` into a new file called `chr22.head.fa`. There are many ways to do this (try `cat chr22.fa | head -10 > chr22.head.fa`).
4. Now let's write the last 10 lines of the `chr22.fa` file into another file and call this `chr22.tail.fa`
5. Make a new directory. Let's call it `headtail`
6. Move `chr22.head.fa` and `chr22.tail.fa` to this directory: `mv chr22.*.fa headtail/`. Note the use of `*`, which is called a "wild-card" and it will match all files that contain rest of the filename characters. In this case both chr22.head.fa and chr22.tail.fa will be moved. Saves typing and errors!
7. Now archive this directory: `tar czfv chr22headtail.tar.gz headtail/*`

### Sizing up chr22
1. Check the sequence at the start and end of chromosome: `head` or `tail` for a defined set of lines; try `more` or `less` and use space to scroll through, ctrl+c to terminate.
2. How many nucleotides are in chr22? Per Genome Reference Consortium, it is **50,818,468**
3. Try `wc chr22.fa`. It gives line, word and character counts, respectively. The character count can give us chr length but it is a bigger number than expected. Why?
4. Use `grep -v ">"` to exclude the first line, `pipe` the output to `wc -c`.
5. Still a bigger number! Why? How about those hidden newline characters?
