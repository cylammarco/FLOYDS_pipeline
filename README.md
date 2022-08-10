# FLOYDS_pipeline

This is an unofficial data reduction pipeline for the [LCO FLOYDS spectrograph](https://lco.global/observatory/instruments/floyds/).

## Fully automated data reduction

To reduce a FLOYDS spectrum already exist in the LCO archive from scratch, we only require user to provide the target name or the coordinates of the targets (and the TNS + LCO API tokens). The main script `floyds_pipeline.py` handles the search, download, paths, and reduction. It works by:

1. Makeing use of the `tns_api_search.py` to resolve name into coordinates (unless coordinates are not provided). It will ask for the API key, BOT ID, BOT name if TNS is needed.
2. Query the LCO archive with `query_lco_archive.py` with the coordinates and find the calibration files at the same time. Optionally the date range of the observation can be provided. It will ask for the LCO API token if needed, put "public" if you want to query public data that doesn't require a token. Or you already have the data in place.
3. Query for the standard observation taken at a closest time from the observation. (The time of the standard observation taken is not limited by the time set to query the science target.)
4. Download the data with `download_floyds_data.py` and sort out the file relations at the same time.
5. yaml configuration files are generated from updating the `floyds_template.yaml` with the target frames.
6. Runs the `reduce_floyds_data.py` to reduce data with ASPIRED.
7. Wait for a few minutes.

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

## Running the reduction manually or automatically

The `floyds_template.yaml` is the template for creating a config file for use with the `reduce_floyds_data.py`. It is automated generated from the `floyds_pipeline.py` if you let it run fully automatically.

The default setting output all the possible intermediate figures and most of data products.

If you wish to reduce the data much finer control, you can manually modify the yaml file and run:

`python reduce_floyds_data.py name_of_the_config_file.yaml`

To run a fully automated reduction with a target name and login details stored in the file login_details.yaml:

`python floyds_pipeline.py --directory="2022juw" --target_name="2022juw" --login="login_details.yaml"`

or

`python floyds_pipeline.py --directory="2022juw" --login="login_details.yaml"`

where it will use the directory name as the target's name.

If you wish to enter your keys explicitly, you can do

`python floyds_pipeline.py --directory="2022juw" --lco-token="xxxxxxxxxxxxxxxx" --tns_bot_id="xxxxxxxxxxxxxxxx" --tns_bot_name="xxxxxxxxxxxxxxxx" --tns_token="xxxxxxxxxxxxxxxx"`

For the same target, if you have the position, we can skip the query through TNS and query directly with the LCO archive, this way the TNS credentials will not be needed:

`python .\floyds_pipeline.py --directory="2022juw" --ra=127.631104 --dec=18.203840 --login="login_details.yaml"`

## If you wish to use the script anywhere

Create an alias for the `reduce_floyds_data.py` file. For example in bash, add this to your `.bashrc` or `.bash_profile`

`alias reduce_floyds="python /path/to/where/you/have/the/script/reduce_floyds_data.py"`

`alias run_floyds_pipeline="python /path/to/where/you/have/the/script/floyds_pipeline.py"`

## More to do

1. More logging
2. Summary pdf
