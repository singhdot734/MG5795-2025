## Using WGCNA and Cytoscape to explore regulatory relationships

 
Hello! Today's activity is going to be focused on coexpression analyses to find relationships between genes and TEs
potentially acting as their regulators. 

First step is to download and install the desktop version of Cytoscape, which can be done below. 
##### 1. Cytoscape install 
```
https://www.https://cytoscape.org/
```

You can do this either when I mention it in lecture or before class. 

Once installed, next thing we're going to do is process coexpression networks in R. Start by initializing an interactive
RStudio session through OSC.

##### 2. Launch RStudio
```
https://ondemand.osc.edu/pun/sys/dashboard/batch_connect/sys/bc_osc_rstudio_server/session_contexts/new
```

RStudio will automatically connect to your OSC account, so your directories will all be available. Navigate to the clone of 
the class git repo using the **Terminal** of the RStudio **Console Pane**, the bottom left window. This acts just like the 
command line. Change directories to the subdirectory for this activity, the same one you used to read this file in Github. 

Now, in **Files** located in the **Viewer Pane** (bottom right), click the .Rproj file. This will initialize the
enviornment and dependencies automatically. Now you can open the *Regulation.Rmd* file by clicking on it in **Files** and 
it will pop up in the **Editor Pane**. Execute everything there before continuing below. 


##### 3. Explore in Cytoscape 
Use FileZilla, Cyberduck, or any other sftp of your choice to move the .zip file to your local machine. Launch Cytoscape,
then import the .sif and nodes.csv files. I'll demonstrate this as well. 
