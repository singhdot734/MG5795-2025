# Installing R packages in OSC classroom R studio app

R packages have to be installed in the classroom directory and not the classroom directory in your home directory, which is the default.
Get paths to libraries where packages can be installed from R studio console: 
```
.libPaths()
[1] "/users/PAS1067/osu8197/osc_classes/MOLGEN_5795_OSU/R"                              
[2] "/fs/ess/PAS3124/MOLGEN_5795_OSU/Rpkgs"                                             
[3] "/apps/spack/0.21/ascend/linux-rhel9-zen2/r/gcc/12.3.0/4.4.0-o2566zh/rlib/R/library"
```
Above, it can be seen that the first path for package installation is in my home directory (not accessible to class).
To remove the first path:
```
.libPaths(.libPaths()[-1])
```
In the above command, -1 can be replaced with -2 or another number corresponding to the path to be removed.

Now do this again:
```
.libPaths()
[1] "/fs/ess/PAS3124/MOLGEN_5795_OSU/Rpkgs"                                             
[2] "/apps/spack/0.21/ascend/linux-rhel9-zen2/r/gcc/12.3.0/4.4.0-o2566zh/rlib/R/library"
```
Now the first path for package installation is the shared directory accessible to all students as well.

Install packages as follows (an example):
`install.packages("BiocManager", lib = "/fs/ess/PAS3124/MOLGEN_5795_OSU/Rpkgs")`
