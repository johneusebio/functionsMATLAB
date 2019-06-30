# MATLAB functions and stuff

This is a fork of EEG functions originally written by Hause Lin (see https://github.com/hauselin/functionsMATLAB). Please ensure that both (EEGLAB)[https://sccn.ucsd.edu/eeglab/index.php] and (ERPLAB)[https://erpinfo.org/erplab] are installed, as they are necessary for the functions provided.

This fork is designed to be more use-accessible and customizable, allowing for easy preprocessing and analysis of EEG/ERP data for a wide variety of study designes with minimal changes to the underlying code.

Most of the necessary functionality can be accessed using the `run_scripts.m` file located within the EEG directory.  The user need simply change the following fields:

- data.path: Enter in the parent directory path for the subject EEG data to be preprocessed and analyzed.
- data.ids: A character array of subject IDs to be used 
- data.exclude: A character array of subject ids to be excluded from the analysis. This is included for easy batching (e.g., subject ids = 1:100, excluding 35).
- data.conds: A structure containing condition information to be used in the analysis. Please ensure that the order in which condition information is provided is consistent (e.g., GO is always first, NOGO is always second). Multiple files per condition per subject (i.e., more than one run) are allowed, but they will be analyzed one file at a time, independent of the others.
- data.conds.name: A character array of condition names you will be using. These names must be entered in the same way they appear in the filenames of the EEG data to be analyzed, as substring comparison will be used to fetch the proper files.
- data.conds.binlist = A character array of the filepath of a text file containing the BINLIST information to be used. 
- data.conds.contrasts = A character array of file paths leading to the a text file containing contrast information for a given condition. Each condition should have its own unique contrast file. If no contrast is necessary, simply leave the corresponding entry within the array blank (i.e., `[]`). Contrast data must be entered in the following format:

    CONTRAST NAME 1: 1 - 2
    CONTRAST NAME 2: 2 - 1

    Each line is a seperate contrast. Numbers correspond to the binlist provided for the given condition. Only simple addition and subtraction of bins is supported. Spaces must seperate each bin number and mathematical operator, as shown above. Contrast names cannot contain colons, as that is treated as a special character.
- data.conds.preproc_cfg: A character array of filepaths pointing to a plain text file containing preprocessing configuration parameters. If not provided (i.e., field is left blank, `[]`), then the default values will be used. See `EEG/s1_preprocessEEG.m` for details. Each line corresponds is treated as a preprocessing parameter to over-ride, and must be entered in using the following format:

    PARAMETER_NAME1 = value1
    PARAMETER_NAME2 = value2
