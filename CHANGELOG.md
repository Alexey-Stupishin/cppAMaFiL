## Update History

* 16 June 2020 - 1<sup>st</sup> release (v 2.1.20.428), [doi:10.5281/zenodo.3896223](https://zenodo.org/record/3896223#.Y13LRHZBxJQ)
* 04 October 2020 (v 2.1.20.1004, rev.363): _weight_bound_size_ key added to _gx_box_make_nlfff_wwas_field.pro_, see comment inside
* 25 January 2021 (v 2.2.21.125, rev.384): bug fixed (_extra parameters could prevent NLFFF in some cases)
* 18 February 2021 (v 2.3.21.217, rev.392):
	* bug fixed (it was crash when calculated with seeds)
	* improvement (all lines calculated, even if seed below chromo_level); but only part above chromo_level is stored
	* comment in _gx_box_calculate_lines.pro_ utility corrected, so that is less confusing
* 16 October 2022 (v 2.4.22.1016, rev.625):
	* small bug for short low loops fixed
	* small lines wrapper correction
	* Linux version started (at separate repository)
* 30 October 2022 (v 3.4.22.1025, rev.626): 
	* Major structure changes 
	* Significant code reorganization for multiplatforming, change implementation to C++11 standard
	* Linux version built
* 23 August 2023
	* Python wrapper
* 30 August 2023 (v 3.4.23.745, rev.725): 
	* Small polishing, version changed
* 03 September 2023
	* Python wrapper update
* 04 December 2023 (v 3.4.23.1203, rev. 797) - One-footpoint-line status bit added
* 01 June 2024 (v 4.0.24.601, rev. 837) - Major update
* 11 June 2024 (v 4.0.24.611, rev. 839) - Polishing
* 01 June 2026 (v 4.4.24.601, rev. 61)
	* Refactoring, reorganization, polishing
	* Non-physical lines don't calculated (time reducing)
	* Calculation time for lines stored (for debugging)
	* Memory requirement  estimation implemented (not fully tested)
	* Disambiguation functionality implemented (not fully tested)
* 03 June 2026 IDL and Python wrappers moved to separate repositories
