# Class 2 
1. Log into OSC (ondemand.osc.edu)
2. Using the menu bar at the top, navigate to Files > /fs/scratch/PAS3124
3. Click the **>_Open in terminal** button at the top to open a new UNIX terminal window in a new tab.
## Basic UNIX commands
4. Where am I?
5. `pwd`
6. What is in this directory?
7. `ls`
8. Create a personal directory in scratch space of the classroom project (PAS3124) 
9. `mkdir yourname`
10. The general structure of unix commands is: program options input output. In this case, program is "mkdir" and output is "unixmore", which is being run as follows (there are no options or input file being specified here):
`mkdir yourname`
11. To enter/change directory to "yourname", use `cd yourname`
12. Now try, `pwd`
13. Let's do the customary first command on UNIX: `echo "Hello world!"`
14. This will write **"Hello world!"** to the screen.
15. Output of any program (or command) can be directed to a file using `>`
16. "Hello world!" can be written to a file (instead of screen) by changing the command on line 13 to `echo "Hello world!" > myfile_1.txt`
17. Now do `ls` again. Now it will show "myfile_1.txt" in your personal directory. This is your first file in UNIX!
18. Try `ls -l`. "-l" is an option that will tell "ls" to list items with more details.
19. To move up a directory: `cd ..`
20. To move up two directories: `cd ../..` ...and so forth.
21. If you ever type just `cd`, this will take you to your home directory on OSC. To come back to your personal directory in the classroom project, you will have to type the **absolute path** to your directory: `cd /fs/scratch/PAS3124/yourname`
22. To get SARS-CoV2 genome as a FASTA file, we will use a command called `curl`.
23. Copy this into your terminal `curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" -o NC_045512.2.fasta`
24. In the above command option `-o` tells what the output of `curl` should be.
25. This will download SARS-CoV2 genome sequence FASTA file into the directory you were in where you ran the code on line 23. To make sure you are in the correct directory, you can do `pwd` before proceeding with line 23.
26. Confirm with `ls` to make sure the downloaded file "NC_045512.2.fasta" is there.
27. Now we will move this FASTA file into a new directory.
28. In your personal directory, use `mkdir` to make a new directory called "sarscov2"
29. Do `ls` to confirm.
30. Change into the "sarscov2" directory using `cd`
31. `cp ../NC_045512.2.fasta .` will copy "NC_045512.2.fasta" file into this directory. `cp` is the program. `../NC_045512.2.fasta` gives the relative path to the file to be copied. `.` means copy to current directory.
32. Instead of copying, you can move the the FASTA file using this command `mv ../NC_045512.2.fasta .`. Now the FASTA file will no longer remain in the directory one level up.

# Class 3
## More basic UNIX commands
1. Go to the "sarscov2" directory in "yourname" directory. You will use `cd`. Tip: Start typing the directory name and then press tab. It will try to complete the name until it hits a unique character. This will save you lot of typing!
2. **Viewing files**
3. `head` or `tail` commands can be used to view contents of a file. `head NC_045512.2.fasta` will show the first 10 lines of the FASTA file. `head -n 20 NC_045512.2.fasta` uses option `-n` with specified variable 20 and will show first 20 lines.
4. `tail NC_045512.2.fasta` will show the last 10 lines of the FASTA file.` How will you look at the last 25 lines of the FASTA file?
5. More ways to view files: `more` or `less` followed by `filename` will print it to the screen in scrollable format. Scroll up or down using arrows, or down by hitting space. Press `q` to quit more and return to command prompt
6. Files can also be viewed using the `cat` command. It will print it all to the screen.
7. **Getting help**
8. `man` followed by program name will open manual for that program to  view options available for each program/command. For example, `man head` will open manual for `head`. Scroll up or down using arrows, or down by hitting space. Press q to quit and return to command prompt.
9. 

In this directory, we will get a file with yeast (S. cerevisiae) genome features and investigate some features of yeast genes and genome using basic unix commands. Originally, these files are available from Saccharomyces Genome Database (SGD; yeastgenome.org) but we will get the files that are hosted by biostar handbook (allegedly, links to files and datasets in SGD can change frequently).
Two commands to get data from http (or ftp) links are `curl` (step 8) or `wget` (step 9). The usage is like any other unix command: program flags input output. We will use `>` to direct the output to a file.
Getting data with `curl` with option `-s` (silent) that will not print progress meter or errors on screen: `curl -s http://data.biostarhandbook.com/data/SGD_features.tab > SGD_features.tab`
Getting data with wget (no flags available): `wget http://data.biostarhandbook.com/data/SGD_features.tab`. wget will automatically save the file, so you don't need `> filename`. To save a file with a different name using wget you need to specify the `-O filename` flag.
Check if the file is now in your directory; use `ls` or `ll`
Now get the README file that explains what is in the SGD_features.tab file: `wget http://data.biostarhandbook.com/data/SGD_features.README` (or use `curl`).
Use `ls` or `ll` to check directory contents again. The README file should also be there now.
View contents of the .tab file. You can use `more` or `less`. Press “space” to move forward, “b” to move backward, q or ESC to exit. Using `head` or `tail` can show 10 lines by default. Again, usage is same: program flags input output. In this case only program and input needs specified. Output will be on your screen. Can you make sense of what type of information is there in the file?
View the README file and read through the description. Now look at a few lines of the .tab file (use `head` or `tail`; use `-n` followed by a small number as a flag) and it may make more sense what all information is in the .tab file. 

