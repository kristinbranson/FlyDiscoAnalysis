#!/usr/bin/python3

import os
import tempfile
import stat

from tpt.utilities import *



def main() :
    url = 'https://github.com/kristinbranson/FlyDiscoAnalysis'
    username_from_user_index = [ 'bransonlab', 'rubinlab', 'projtechreslab', 'geniegeneric', 'otopalikrobot' ]
    clone_and_copy_github_repository_into_user_home_folders(url, username_from_user_index, branch_name='deploy')



if __name__ == "__main__":
    main()
