## Instructions
---

1. compile gambit as normal

2. copy the gambit executable into files/
    (overwrite the existing one)

3. setup the options below to you liking

4. set any other options in the yaml files
    found in files/yaml_files/ but don't touch
    any line containing '~'

5. setup your desired plots in files/plots.pip
    but don't touch any line containing '~'

6. modify files/job.sh so that it works on
    your HPC but don't touch any line containing '~'

7. run this python script

8. run the file: 'gen_path/runScans.sh'

9. run the file: 'gen_path/runPippi.sh'


a GAMBIT is generated for all possible combinations below and a script is generated for running all simultaneously
WARNING: be careful not to generate too many, otherwise you will run out of storage