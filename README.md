# FLOYDS_pipeline

This is an unofficial data reduction pipeline for the [LCO FLOYDS spectrograph](https://lco.global/observatory/instruments/floyds/). 

See the `floyds_default.yaml` to see the example config file for reducing AT2022DML.

The output folder is defined in `output_folder_path`, it is a relative path to where the config yaml file is.

The default setting output all the possible intermediate figures and intemediate data products.

## Running the script

`python reduce_floyds_data.py name_of_the_config_file.yaml`

## If you wish to use the script anywhere

Create an alias for the `reduce_floyds_data.py` file. For example in bash, add this to your `.bashrc` or `.bash_profile`

`alias reduce_floyds="/path/to/where/you/have/the/script/reduce_floyds_data.py"`
