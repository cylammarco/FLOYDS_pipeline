# FLOYDS_pipeline

This is an unofficial data reduction pipeline for the [LCO FLOYDS spectrograph](https://lco.global/observatory/instruments/floyds/). 


## Fully automated data reduction

The `floyds_pipeline.py` makes use of the `tns_api_search.py` to resolve name into coordinates, then query the LCO archive with `query_lco_archive.py` and download with `download_floyds_data.py`. After downloading all the relevant science and standard data within the defined date range, it reduce the data by calling the script `reduce_floyds_data.py`.

At the command line, these are the options, brackets show the default values.

```
--target_name:  The target name to be queried on the TNS. (None)
--ra:           Right Ascension in decimal. Only used if target_name is None. (None)
--dec:          Declination in decimal. Only used if target_name is None. (None)
--directory:    Path to store the raw and reduced data products. (None)
--login:        Path to the login details. (None)
--lco_token:    LCO token. Only used if --login is None. (None)
--tns_bot_id:   TNS Bot ID. Only used if --login is None. Not used if ra and dec are provided. (None)
--tns_bot_name: TNS Bot name. Only used if --login is None. Not used if ra and dec are provided. (None)
--tns_token:    TNS token. Only used if --login is None. Not used if ra and dec are provided. (None)
--date_start:   The date of the beginning of the night of the observation. ("1900-01-01")
--date_end:     The date of the beginning of the night of the observation. ("2100-12-31")
```

## Reduction

The `floyds_template.yaml` is the template for creating a config file for use with the `reduce_floyds_data.py`. It is automated generated from the `floyds_pipeline.py` if you let it run fully automatically.

The default setting output all the possible intermediate figures and most of data products.

## Running the script

To run a reduction with a custom-built configuration script:

`python reduce_floyds_data.py name_of_the_config_file.yaml`

To run a fully automated reduction with a target name and login details stored in the file login_details.yaml:

`python floyds_pipeline.py --directory="2022juw" --target_name="2022juw" --login="login_details.yaml"`

or

`python floyds_pipeline.py --directory="2022juw" --login="login_details.yaml"`

where it will use the directory name as the target's name.


For the same target, if you have the position, we can skip the query through TNS and query directly with the LCO archive, this way the TNS credentials will not be needed:

`python .\floyds_pipeline.py --directory="2022juw" --ra=127.631104 --dec=18.203840 --login="login_details.yaml"`

## If you wish to use the script anywhere

Create an alias for the `reduce_floyds_data.py` file. For example in bash, add this to your `.bashrc` or `.bash_profile`

`alias reduce_floyds="/path/to/where/you/have/the/script/reduce_floyds_data.py"`
`alias floyds_pipeline="/path/to/where/you/have/the/script/floyds_pipeline.py"`

## More to do 

1. More logging
2. Summary pdf
