Plotting multiple hazard curves:

* First set PYTHONPATH in environmental variables (one time only), e.g.:

	export PYTHONPATH=~/2015_gsc_nshm/pycodes/tools/oq_tools:$PYTHONPATH
	
* To run script from any folder, set alias (Optional - one time only), e.g.:

	alias hcp='python ~/2015_gsc_nshm/pycodes/post_processing/cmp_multi_oq_haz_curves.py'

* Make configuration file.  See link below for example:

	https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/zones/2012_mw_ge_4.0/NSHA13_Background/output/cmp_bst_mxv_bvv_col_curves.conf
	
	
* Set "siteListFile" to:
	
	~/2015_gsc_nshm/shared/swcan_sites_cities.csv (for example)

* "relativePaths" are folder paths relative to the config file

* To run:

	> python ~/2015_gsc_nshm/pycodes/post_processing/cmp_multi_oq_haz_curves.py <path/to/config/file>
	
	or, if alias is set:
	
	> hcp <path/to/config/file>