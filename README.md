# reliosparc - tools for interfacing between cryoSPARC and Relion

Tristan Bell

Chao Group

Mass General Hospital and Harvard Medical School

## Transferring particle coordinates from cryosparc to relion 

The cs_to_stars.py script can take a cryosparc *.cs file, extract the particle coordinates, and format them into star files for easy import into relion.  **Please note** that this script only transfers particle coordinates.  If you would rather extract a particle stack in cryosparc and process it further in relion, the csparc2star script in David Asarnow's PyEM repository works well.

**Dependencies**

 - Python 3.5+ with the Numpy package.  Numpy can be installed with

    python -m pip install numpy

 - Your micrographs must have an invariant initial portion and a variant final portion.  I.e., mcg1, mcg2, ... , mcg999 will work, but 1mcg, 2mcg, ... , 999mcg will not.

**Accepted input file formats**

The script requires a cryosparc *.cs particle or particle passthrough file.  If one file gives an error, try using a different particle file from the same job.  Not all cryosparc particle files contain the necessary coordinate info.

The following cryosparc file types have been validated:

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

**Running the Script**

Run the script specifying your input cryosparc file:

    python cs_to_stars.py path/to/your/cryosparc/file.cs

A folder called Raw_data will be created in your working directory that contains a series of star files (one for each micrograph).

**Using the star outputs**

 - The star outputs mimic those of a relion manual picking or autopicking job.
 - In relion, run a small manual picking or autopicking job, specifying all of your CTF-corrected micrographs as the template micrograph set.
 - After completing the job, navigate to the relion job folder for the picking job.  The relion job folder will contain a subdirectory called "Raw_data" that contains star files specifying the picked particle coordinates.
 - Delete the relion Raw_data folder and replace it with the Raw_data folder created by the python script.
 - You're done!  Any further jobs in relion using this job will pull the coordinates defined in your cryosparc file.  Happy processing!

## Bugs & Troubleshooting

Feel free to reach out with any issues you encounter by opening a new issue item in the GitHub.  I will try to respond promptly to any issues.

