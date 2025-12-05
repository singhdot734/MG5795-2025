## Using WGCNA and Cytoscape to explore regulatory relationships

 
Hello! Today's activity is going to be focused on coexpression analyses to find relationships between genes and TEs
potentially acting as their regulators. 

First step is to download and install the desktop version of Cytoscape, which can be done below. 
##### 1. Cytoscape install 
```
www.cytoscape.org
```

You can do this either when I mention it in lecture or before class. 

Once installed, next thing we're going to do is process coexpression networks in R. Start by initializing an interactive
RStudio session through OSC.

##### 2. Launch RStudio
Login OSC Classroom with the following link:
https://ondemand.osc.edu/pun/sys/dashboard/batch_connect/sys/bc_osc_rstudio_server/session_contexts/new

Put "2" in the **Number of hours** to apply for an 2-hour session. After submitting the request, you should see the status "Queued". 
After the request is approved, the status changes to **Running**. Click the **Connect to RStudio Server** button to open RStudio on OSC.   

RStudio will automatically connect to your OSC account, so your directories will all be available. 
To make sure you are at your root directory, at your R **Console** tab of the RStudio **Console Pane** (the bottom left window), run this 
command to set your directory to your home directory at OSC. 
```R
setwd("~")
```

Then, navigate to the clone of the class git repo using the **Terminal** of the RStudio **Console Pane**. This acts just like the 
command line. Change directories to the subdirectory for this activity, the same one you used to read this file in Github. 

Update your course folder
```bash
cd ~/MG5645/name.#/MG5795-2025
git pull
```

Now, in **Files** located in the **Viewer Pane** (bottom right window), click the `27.TEs_and_regulation.Rproj` file. Select **General** 
if it pops up a "Project Options" window. This will initialize the enviornment and dependencies automatically. Now you can open the 
`Regulation.Rmd` file by clicking on it in **Files** and it will pop up in the **Editor Pane** (upper left window). Execute everything there 
before continuing below. 


##### 3. Explore in Cytoscape 
Now you prepared network file (.sif) and annotation file (.csv) for the second coexpression module (M2), which are stored in the `cyto_ME(M2)` folder.   

Use FileZilla or OSC OnDdmand to download the `ME(M2)_cytoscape_bundle.zip` file in the `cyto_ME(M2)` folder to your local machine. 
Unzip the file, then you see 'ME(M2)_edges_LTR_Gene.csv', 'ME(M2)_network.sif', 'ME(M2)_nodes.csv'.   

Launch Cytoscape, then import the *.sif and *nodes.csv files. I'll demonstrate this as well. 


