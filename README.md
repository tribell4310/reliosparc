# reliosparc - interfacing between cryoSPARC and Relion

Tristan Bell

Chao Group

Massachusetts General Hospital and Harvard Medical School

## Transferring particle coordinates from cryosparc to relion 

The cs_to_stars.py script can take a cryosparc *.cs file, extract the particle coordinates, and format them into star files for easy import into relion.  **Please note** that this script only transfers particle coordinates.  If you would rather extract a particle stack in cryosparc and process it further in relion, the csparc2star script in Daniel Asarnow's PyEM repository works well.



**Dependencies**

Python 3.5+ with the Numpy package.  Numpy can be installed with

    python -m pip install numpy



**Accepted input file formats**

The script requires a cryosparc *.cs particle or particle passthrough file.  If one file gives an error, try using a different particle file from the same job.  Not all cryosparc files contain particle coordinate info.

The following cryosparc 3.x.x file types have been validated:

| Cryosparc Job Type | Preferred Filetype | Example Filename |
|--|--|--|
| 2D Classification | Particle Passthrough File | P45_J103_passthrough_particles.cs |
| Select Particles | Particle Passthrough File | P45_J104_passthrough_particles_selected.cs |
| Ab Initio | Particle Passthrough File | P45_J123_passthrough_particles_all_classes.cs |
| Refinement | Particle Passthrough File | P45_J124_passthrough_particles.cs |
| Particle Inspection | Extracted Particle File | extracted_particles.cs |
| Particle Extraction | Particle Passthrough File | P45_J136_passthrough_particles.cs |
| Blob Picker | Picked Particle File | picked_particles |
| Topaz Extract | Topaz Picked Particle File | topaz_picked_particles.cs |
| Remove Duplicates |Particles Kept or Particles Excluded File | particles_kept.cs |

The script also requires the corresponding "micrographs_ctf.star" file from your relion job.  If you have imported the same micrographs into relion as you did into cryosparc (then performed motion correction and Ctf correction in relion), this file can be found in the folder for your relion CtfFind job.

*Note:* For cryosparc 4.x.x, use a particle extraction job, and the extracted_particles.cs file in the job directory.



**Running the Script**

Run the script specifying your input cryosparc file and micrographs_ctf.star file:

    python cs_to_stars.py --cs path/to/your/cryosparc/file.cs --star /path/to/your/micrographs_ctf.star

A folder called Raw_data will be created in your working directory and populated with a star files (one for each micrograph).

**Options**

| Flag | Description |
|--|--|
| --flipx | inverts particle coordinates in x dimension |
| --flipy | inverts particle coordinates in y dimension |
| --no_nudge | overrides default coordinate identification behavior for cs files (*not recommended for most users*, only consider if prompted by the script at runtime)



**Using the star outputs**

 - The star outputs mimic those of a relion autopicking job.
 - In relion, run a small autopicking job, specifying all of your CTF-corrected micrographs as the template micrograph set.
 - After completing the job, navigate to the relion job folder for the picking job.  The relion job folder will contain a subdirectory called "Raw_data" that contains star files specifying the picked particle coordinates.
 - Delete the relion Raw_data folder and replace it with the Raw_data folder created by the python script.  (If you're transferring the files over a network, wrapping all the files up into a zip or tarball file will dramatically accelerate transfer.)
 - You're done!  Any further jobs in relion using this job will pull the coordinates defined in your cryosparc file.  Happy processing!



## Bugs & Troubleshooting

Feel free to reach out with any issues you encounter by opening a new issue item in the GitHub.


## Other Goodies

The converge.py and converge_continue.py scripts track movement of particles between classes during 2D or 3D classification jobs in relion.  converge.py is used for regular classification jobs, converge_continue.py is used for "continue" jobs.  The scripts produce png plots of class populations and a converge.log logfile describing overall volatility in the class assignments.  It can be executed during a run (to track convergence of classes) or after a run has finished.  Just place the script into your relion job folder and run as follows:

    python converge.py numberOfIterations &

... or ...

    python converge_continue.py numberOfInitialIterations numberOfTotalIterations &

