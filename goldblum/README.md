Goldblum
========

Goldblum is a Matlab program that collects FlyDisco experiments off of
a set of rig computers and analyzes them.  It is intended to be run
automatically on a schedule, typically once a day, at night after data
collection for the day has concluded.

Each lab that uses Goldblum runs a single Goldblum 'instance'.  Each
instance runs as a particular user, typically a special user created
to own data that is collectively owned by a particular lab.  For
instance, the Branson Lab Goldblum instance runs as the user
`bransonlab`.

Currently, Goldblum runs are launched from a cron job on
`submit.int.janelia.org`.  There is one cron job per lab, i.e. one
cron job per Goldblum instance.

When a Goldblum instance runs, it searches a particular folder on each
specified rig computer for `experiment` folders.  These are folders
containing two or more of the files `movie.ufmf`, `metadata.xml`, and
`ABORTED`.  (The first two file names are case-insensitive.  The last
one is case-sensitive.)  Any such folders are copied to a specified
folder on `dm11`, and then they are subsequently analyzed on the
Janelia cluster using analysis settings determined by the contents of
the `metadata.xml` folder.

Whether a folder is considered an experiment folder is decided by a point system.
Each of the following is worth a point:

1.  If the folder contains a file named `movie.ufmf` or `movie.avi` (case-insensitive).
2.  If the folder contains a file named `metadata.xml` (case-insensitive).
3.  If the folder contains a file names `ABORTED` (case-sensitive).

If a folder scores two points or more, it is considered an experiment folder.

Goldblum is part of the FlyDiscoAnalysis project, the source code for
which lives at
```
https://github.com/JaneliaSciComp/FlyDiscoAnalysis
```
Goldblum is implemented by the file
```
goldblum/goldblum.m
```
within the source repository.  This is the function that is run by the
cron job. See the documentation within that file for more details.

The core analysis pipeline is implemented by the file
```
FlyDiscoPipeline.m
```
within the source repository.  See the documentation within that file
for more details.

Folders of Goldblum
-------------------

Each Goldblum instance copies experiments from a set of rig computers to a
destination folder, typically named `flydisco_data` somewhere under
`/groups/<lab_head_name>` on the dm11 file system.  Within this folder, Goldblum
creates a folder named `goldblum-logs`.  For each run of Goldblum, a file is
created which records the Goldbum log for that run.  These files have names like
`goldblum-2021-01-09.log`.  These files will contain error messages if something
goes wrong with a Goldblum run.  It's also a good idea to check these files occasionally to 
make sure Goldblum is running correctly.

Another folder worth knowing about is the `to-process` folder, also within 
the `flydisco_data` folder.  After an experiment folder is transferred from 
one of the rig computers, a symbolic link is created within the `to-process` 
folder.  This link points to the just-transfered experiment folder.  When Goldblum
finishes transfering experiments from the rig computers, it uses the contents of
the `to-process` folder to determine what experiment folders need to be
analyzed.  Note that if you manually create a symbolic link to an experiment folder 
in the `to-process` folder, that experiment will get analyzed just as if it had been 
transfered by Goldblum.  This is often useful for reanalyzing experiment folders when 
something fixable goes wrong with the original attempt at analysis.


How to set up Goldblum for a new lab
------------------------------------

Suppose a new lab head joined Janelia, named Gina Davis.  To set up
Goldblum for the Davis Lab, one would:

1.  On each rig machine you want to transfer data off of, install Cygwin and the
    Cygwin ssh server. We currently use Cygwin 1.7.1, which is super old, but a
    more recent version of Cygwin should also work.

2.  Have SciComp Systems create a `davislab` user.

3.  Have SciComp Systems add your personal public RSA key to the list of
    authorized keys for the account.  You may need to create these if you don't
    have one already.  This will enable you to login as `davislab` on
    `submit.int.janelia.org` without entering a password. For instance, when
    logged in as yourself, you should be allowed to run `ssh davislab@submit.int.janelia.org`.

4.  Create an RSA keypair for the `davislab` account using `ssh-keygen`.  Leave
    the passphrase empty.

5.  On your local macOS or Linux machine, which is logged into your normal account, start a terminal and do:
    ```
    git clone https://github.com/JaneliaSciComp/FlyDiscoAnalysis
    cd FlyDiscoAnalysis
    git submodule update
    ```
    This will create a folder named `FlyDiscoAnalysis` with the FlyDiscoAnalysis source code in it,
    including Goldblum. 

6.  In your local repo, there is a file named `goldblum/branson_configuration.m`.  This file looks
    something like this:
    ```
    function result = branson_configuration()
        result = struct() ;
        result.lab_head_last_name = 'branson' ;
        result.host_name_from_rig_index = { 'arrowroot.hhmi.org', 'beet.hhmi.org', 'carrot.hhmi.org', 'daikon.hhmi.org' } ;
        result.rig_user_name_from_rig_index = repmat({'bransonk'}, [1 4]) ;
        result.data_folder_path_from_rig_index = repmat({'/cygdrive/e/flydisco_data'}, [1 4]) ;
        result.destination_folder = '/groups/branson/bransonlab/flydisco_data' ;
        this_folder_path = fileparts(mfilename('fullpath')) ;
        flydisco_analysis_path = fileparts(this_folder_path) ;
        result.settings_folder_path = fullfile(flydisco_analysis_path, 'settings') ;
        result.does_use_per_user_folders = false ;
    end
    ```
    Make a copy of this file (in the `goldblum/` folder) and rename in `davis_configuration.m`.
    Modify the file contents to look something like this, as appropriate:
    ```
    function result = davis_configuration()
        result = struct() ;
        result.lab_head_last_name = 'davis' ;
        result.host_name_from_rig_index = { 'davis-rig-1.hhmi.org', 'davis-rig-2.hhmi.org' } ;
        result.rig_user_name_from_rig_index = repmat({'labadmin'}, [1 2]) ;
        result.data_folder_path_from_rig_index = repmat({'/cygdrive/e/flydisco_data'}, [1 4]) ;
        result.destination_folder = '/groups/davis/davislab/flydisco_data' ;
        this_folder_path = fileparts(mfilename('fullpath')) ;
        flydisco_analysis_path = fileparts(this_folder_path) ;
        result.settings_folder_path = fullfile(flydisco_analysis_path, 'settings') ;
        result.does_use_per_user_folders = false ;
    end
    ```    
    The `lab_head_last_name` field is used for two things: it determines what
    lab is billed for the analysis job on the cluster, and it determines what
    folder on each rig is searched for experiments (see the description of
    `data_folder_path_from_rig_index` below).

    The `host_name_from_rig_index` field determines the host name of
    each rig that will have experiments transfered off of it and analyzed.

    The `rig_user_name_from_rig_index` gives, for each rig machine,
    the username that Goldblum will use to connect.  In this example,
    Goldblum runs as user `davislab` on a cluster node.  Each rig user
    account must have `davislab`s public RSA key added to its
    `~/.ssh/authorized_keys` file, so that Goldblum, running as the
    `davislab` user, can connect without using a password.  Adding a
    key to the `authorized_keys` file is usually done using the
    `ssh-copy-id` command.

    The `data_folder_path_from_rig_index` specifies, for each rig machine, the
    folder where Goldblum will look for experiment folders.  Within this folder,
    Goldblum expects to find a folder with the same name as the
    `lab_head_last_name`, and this `lab_head_last_name` folder is searched for
    experiment folders.

    The `settings_folder_path` field specifies where Goldblum will
    look for analysis settings folders.  Leave this alone unless you
    know what you're doing.

    The `does_use_per_user_folders` field specifies whether the folder
    specified by `data_folder_path_from_rig_index` on each rig will
    contain per-user folders that need to be kept separate when
    transfered to `dm11`.  For instance, if this field is `true`, and
    each rig folder contains the subfolders `alice/` and `bob/`, then
    the `destination_folder` on `dm11` will also contain subfolders
    `alice/` and `bob/`, with experiments from the matching subfolder
    on the rigs.  If `does_use_per_user_folders` is `false`, then
    experiments from the `alice/' and 'bob/' subfolders would all be
    mixed together in the `destination_folder` on `dm11`.

7.  On your local machine, start Matlab and change the working directory to the `FlyDiscoAnalysis`
    folder you created above.
    
8.  At the Matlab command line, do:
    ```
    modpath
    ```
    This will update your Matlab path so that you can access all the code in FlyDiscoAnalysis.
    (The path is only update for the current session.  Do `savepath` if you want to make this
    change permanent.)

9.  At the Matlab command line, do:
    ```
    edit push_goldblum_into_production.m
    ```
    
10. In the Matlab editor, you will see the function looks something like this:
    ```
    function push_goldblum_into_production()
        % Determine the FlyDiscoAnalysis folder path
        goldblum_folder_path = fileparts(mfilename('fullpath')) ;
        fda_folder_path = fileparts(goldblum_folder_path) ;

        % Make sure there are no uncommitted changes
        error_if_uncommited_changes(fda_folder_path) ;

        % Do Branson Lab instance
        copy_to_single_user_account('bransonlab', fda_folder_path) ;

        % Do Rubin Lab instance
        copy_to_single_user_account('rubinlab', fda_folder_path) ;

        % If get here, everything went well
        fprintf('Successfully copied %s into all the *lab user accounts\n', fda_folder_path) ;
    end
    ```
    If you are an administrator, and want to update many lab's code, you want to have many labs
    listed in this file. Add a line for the Davis Lab.  Afterwards, the text should look
    something like this:
    ```
    function push_goldblum_into_production()
        % Determine the FlyDiscoAnalysis folder path
        goldblum_folder_path = fileparts(mfilename('fullpath')) ;
        fda_folder_path = fileparts(goldblum_folder_path) ;

        % Make sure there are no uncommitted changes
        error_if_uncommited_changes(fda_folder_path) ;

        % Do Branson Lab instance
        copy_to_single_user_account('bransonlab', fda_folder_path) ;

        % Do Rubin Lab instance
        copy_to_single_user_account('rubinlab', fda_folder_path) ;

        % Do Davis Lab instance
        copy_to_single_user_account('davislab', fda_folder_path) ;

        % If get here, everything went well
        fprintf('Successfully copied %s into all the *lab user accounts\n', fda_folder_path) ;
    end
    ```
    If you just want to update code for your own lab, make a copy of this function with only
    you lab in it, and run. The text should look something like this:
    ```
    function push_goldblum_into_production_davis()
        % Determine the FlyDiscoAnalysis folder path
        goldblum_folder_path = fileparts(mfilename('fullpath')) ;
        fda_folder_path = fileparts(goldblum_folder_path) ;

        % Make sure there are no uncommitted changes
        error_if_uncommited_changes(fda_folder_path) ;

        % Do Davis Lab instance
        copy_to_single_user_account('davislab', fda_folder_path) ;

        % If get here, everything went well
        fprintf('Successfully copied %s into all the *lab user accounts\n', fda_folder_path) ;
    end
    ```

11. At the Matlab command line, do:
    ```
    push_goldblum_into_production
    ```
    Behind the scenes, this will use `scp` to copy your `FlyDiscoAnalysis` folder to
    each of the `bransonlab`, `rubinlab`, and `davislab` accounts.

12. Back at the macOS/Linux terminal (logged in as your normal user account), do:
    ```
    ssh davislab@submit.int.janelia.org
    ```

13. You should see the `FlyDiscoAnalysis` folder in your home folder. Do:
    ```
    cd FlyDiscoAnalysis
    ```
    to `cd` into it.

14. At the terminal, do:
    ```
    /misc/local/matlab-2019a/bin/matlab -singleCompThread -nodisplay -batch 'modpath; turn_on_goldblum()'
    ```
    This will add a line to the `davislab` user crontab to launch Goldblum at
    10&nbsp;PM every night.  To instead run Goldblum at, say, 11:45&nbsp;PM, do:
    ```
    /misc/local/matlab-2019a/bin/matlab -singleCompThread -nodisplay -batch 'modpath; turn_on_goldblum(23,45)'
    ```    
    (You may get a message like:
    ```
    There was a problem during execution of the set_parpool_job_storage_location() function:
    Error using mkdir
    Permission denied
    ```
    but this is normal, and you can ignore it.)

15. To see if that worked, do:
    ```
    crontab -l
    ```
    You should see something like this:
    ```
    MAILTO=""

    00 22 * * *     flock --nonblock '/groups/davis/davislab/flydisco_data' --command ''\''/groups/davis/home/davislab/FlyDiscoAnalysis/goldblum/goldblum_launcher.sh'\'' '\''/groups/davis/home/davislab/.bash_profile'\'' '\''/groups/davis/home/davislab/FlyDiscoAnalysis'\'' davis '\''/groups/davis/davislab/flydisco_data/goldblum-logs'\'''   #GOLDBLUM
    ```
    Note that the cron job launches Goldblum using the `goldblum/goldblum_launcher.sh` Bash script,
    which then launches Matlab and runs then `goldblum()` function within Matlab.

16. Make sure the folder for the Goldblum logs exists by doing:
    ```
    mkdir -p /groups/davis/davislab/flydisco_data/goldblum-logs
    ```
    (You *have* to do this as the `davislab` user.)

17. If you want to turn Goldblum off, do this:
    ```
    /misc/local/matlab-2019a/bin/matlab -singleCompThread -nodisplay -batch 'modpath; turn_off_goldblum()'
    ```

    (Again, you may get a message like:
    ```
    There was a problem during execution of the set_parpool_job_storage_location() function:
    Error using mkdir
    Permission denied
    ```
    but this is normal, and you can ignore it.)
