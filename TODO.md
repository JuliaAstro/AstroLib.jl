TODO List
=========

Here is the list of procedures and function provided by IDL AstroLib.  Not all
utilities listed under "Missing" need to be translated, they can be of no or
little utility or deprecated in original AstroLib.  Many functions really don’t
fit into this package (like mathematical ones or plotting-related procedures)
and are kept in the "Missing" section only for reference.

Our priority is to port the "Astronomical utilities", and we should avoid
overlap with other JuliaAstro projects (e.g., don't port cosmology or
FITS-related functions, they're already part of `Cosmology.jl` and `FITSIO.jl`).

You can also help reporting functions in IDL AstroLib already present elsewhere
in Julia ecosystem (like in Julia standar library or third-party packages), so
that we don’t duplicate our efforts.

Deprecated or Not Needed
------------------------

* `findpro`

Already Present in Julia
------------------------

* `asinh`
* `cirrange`.  It is equivalent to `mod(x, 360)`,
  or to `mod2pi(x)` for the `[0, 2pi)` range.
* `minmax`.  It is called `extrema` in Julia.
* `permute`.  It is called `randperm` in Julia.
* `to_hex`.  It is called `hex` in Julia.

Present in Other Libraries
--------------------------

* `aper`, see <https://github.com/kbarbary/AperturePhotometry.jl>
* `cosmo_param`, see `Cosmology` package
  (<https://github.com/JuliaAstro/Cosmology.jl>)
* `galage`, see `Cosmology` package
  (<https://github.com/JuliaAstro/Cosmology.jl>)
* `glactc_pm`, see `SkyCoords` package
  (<https://github.com/kbarbary/SkyCoords.jl>)
* `glactc`, see `SkyCoords` package (<https://github.com/kbarbary/SkyCoords.jl>)
* `jplephinterp`, see `JPLEphemeris.jl` package
  (<https://github.com/helgee/JPLEphemeris.jl>)
* `jplephread`, see `JPLEphemeris.jl` package
  (<https://github.com/helgee/JPLEphemeris.jl>)
* `jplephtest`, see `JPLEphemeris.jl` package
  (<https://github.com/helgee/JPLEphemeris.jl>)
* `lumdist`, see `Cosmology` package
  (<https://github.com/JuliaAstro/Cosmology.jl>)
* `readcol`, use `DelimitedFiles.readdlm`

Missing in AstroLib.jl
----------------------

### Astronomical utilities

* `aitoff_grid`
* `arcbar`
* `arrows`
* `astdisp`
* `astro`
* `ccm_unred` (should go to [DustExtinction.jl](DustExtinction.jl))
* `date`
* `date_conv`
* `eqpole_grid`
* `fm_unred`
* `gal_flat`
* `get_coords`
* `imcontour`
* `qdcb_grid`
* `tdb2tdt`
* `ticlabels`
* `zang`

### DAOPHOT-Type Photometry Procedures

* `cntrd`
* `dao_value`
* `daoerf`
* `find`
* `gcntrd`
* `getpsf`
* `group`
* `mmm`
* `nstar`
* `pixwt`
* `pkfit`
* `rdpsf`
* `sky`
* `srcor`
* `substar`
* `t_aper`
* `t_find`
* `t_getpsf`
* `t_group`
* `t_nstar`
* `t_substar`

### Database Procedures

* `db_ent2ext`
* `db_ent2host`
* `db_info`
* `db_item`
* `db_item_info`
* `db_or`
* `db_titles`
* `dbbuild`
* `dbcircle`
* `dbclose`
* `dbcompare`
* `dbcreate`
* `dbdelete`
* `dbedit`
* `dbedit_basic`
* `dbext`
* `dbext_dbf`
* `dbext_ind`
* `dbfind`
* `dbfind_entry`
* `dbfind_sort`
* `dbfparse`
* `dbget`
* `dbhelp`
* `dbindex`
* `dbindex_blk`
* `dbmatch`
* `dbopen`
* `dbprint`
* `dbput`
* `dbrd`
* `dbsearch`
* `dbsort`
* `dbtarget`
* `dbtitle`
* `dbupdate`
* `dbval`
* `dbwrt`
* `dbxput`
* `dbxval`
* `imdbase`

### Disk I/O (e.g. IRAF files)

* `irafdir`
* `irafrd`
* `irafwrt`
* `read_fmr`
* `wfpc2_read`

### FITS Header Astrometry

* `ad2xy`
* `add_distort`
* `adxy`
* `cons_dec`
* `cons_ra`
* `extast`
* `fits_cd_fix`
* `get_equinox`
* `getrot`
* `gsss_stdast`
* `gsssadxy`
* `gsssextast`
* `gsssxyad`
* `hastrom`
* `hboxave`
* `hcongrid`
* `heuler`
* `hextract`
* `hprecess`
* `hrebin`
* `hreverse`
* `hrot`
* `hrotate`
* `make_astr`
* `putast`
* `sip_eval`
* `solve_astro`
* `starast`
* `tnx_eval`
* `tpv_eval`
* `update_distort`
* `wcs_check_ctype`
* `wcs_demo`
* `wcs_getpole`
* `wcs_rotate`
* `wcssph2xy`
* `wcsxy2sph`
* `wfpc2_metric`
* `xy2ad`
* `xyad`
* `xyxy`

### STSDAS Image manipulation

* `extgrp`
* `st_diskread`
* `sxginfo`
* `sxgpar`
* `sxgread`
* `sxhcopy`
* `sxhmake`
* `sxhread`
* `sxhwrite`
* `sxmake`
* `sxopen`
* `sxread`
* `sxwrite`

### FITS ASCII & Binary Table I/O

* `ftab_delrow`
* `ftab_ext`
* `ftab_help`
* `ftab_print`
* `ftaddcol`
* `ftcreate`
* `ftdelcol`
* `ftdelrow`
* `ftget`
* `fthelp`
* `fthmod`
* `ftinfo`
* `ftkeeprow`
* `ftprint`
* `ftput`
* `ftsize`
* `ftsort`
* `tbdelcol`
* `tbdelrow`
* `tbget`
* `tbhelp`
* `tbinfo`
* `tbprint`
* `tbsize`

### FITS Binary Table Extensions I/O

* `fxaddpar`
* `fxbaddcol`
* `fxbclose`
* `fxbcolnum`
* `fxbcreate`
* `fxbdimen`
* `fxbfind`
* `fxbfindlun`
* `fxbfinish`
* `fxbgrow`
* `fxbheader`
* `fxbhelp`
* `fxbhmake`
* `fxbintable`
* `fxbisopen`
* `fxbopen`
* `fxbparse`
* `fxbread`
* `fxbreadm`
* `fxbstate`
* `fxbtdim`
* `fxbtform`
* `fxbwrite`
* `fxbwritm`
* `fxfindend`
* `fxhclean`
* `fxhmake`
* `fxhmodify`
* `fxhread`
* `fxpar`
* `fxparpos`
* `fxread`
* `fxwrite`

### FITS I/O

Use FITSIO.jl or FITSFiles.jl

* `check_fits`
* `fits_add_checksum`
* `fits_ascii_encode`
* `fits_close`
* `fits_help`
* `fits_info`
* `fits_open`
* `fits_read`
* `fits_test_checksum`
* `fits_write`
* `fitsdir`
* `fitsrgb_to_tiff`
* `fxmove`
* `fxposit`
* `headfits`
* `mkhdr`
* `modfits`
* `mrd_hread`
* `mrdfits`
* `mwrfits`
* `rdfits_struct`
* `readfits`
* `sxaddhist`
* `sxaddpar`
* `sxdelpar`
* `sxpar`
* `writefits`

### Image Manipulation

* `boxave`
* `convolve`
* `correl_images`
* `correl_optimize`
* `corrmat_analyze`
* `cr_reject`
* `dist_circle`
* `dist_ellipse`
* `filter_image`
* `frebin`
* `imlist`
* `max_entropy`
* `max_likelihood`
* `medarr`
* `positivity`
* `psf_gaussian`
* `rinter`
* `sigma_filter`
* `skyadj_cube`
* `xmedsky`

### Math and Statistics

* `avg` - Use `Statistics.mean`
* `cic`
* `cspline`
* `factor`
* `fitexy`
* `flegendre`
* `gaussian` - Use `Distributions.Normal`
* `hermite`
* `ksone`
* `kstwo`
* `kuiperone`
* `kuipertwo`
* `linmix_err`
* `linterp`
* `meanclip`
* `minf_bracket`
* `minf_conj_grad`
* `minf_parabol_d`
* `minf_parabolic`
* `mlinmix_err`
* `mrandomn`
* `multinom`
* `ngp`
* `pca`
* `pent`
* `poidev`
* `polint`
* `poly_smooth`
* `polyleg`
* `prob_ks`
* `prob_kuiper`
* `qsimp`
* `qtrap` - Use Trapz.jl
* `quadterp`
* `randomchi` - Use Distributions.jl
* `randomdir` - Use Distributions.jl
* `randomgam` - Use Distributions.jl
* `randomp`
* `randomwish`
* `safe_correlate`
* `sixlin`
* `tabinv`
* `transform_coeff`
* `trapzd`
* `tsc`
* `tsum`
* `zbrent`

### Plotting Procedures

Use any of Plots.jl, Makie.jl, or UnicodePlots.jl

* `al_legend`
* `al_legendtest`
* `cleanplot`
* `lineid_plot`
* `multiplot`
* `oploterror`
* `partvelvec`
* `ploterror`
* `plothist`
* `plotsym`
* `rdplot`
* `sunsymbol`
* `vsym`

### IDL Structure procedures

* `compare_struct`
* `copy_struct`
* `copy_struct_inx`
* `create_struct`
* `mrd_struct`
* `print_struct`
* `tag_exist`
* `where_tag`

### Robust Statistics

* `autohist`
* `biweight_mean`
* `histogauss`
* `medsmooth`
* `resistant_mean`
* `rob_checkfit`
* `robust_linefit`
* `robust_poly_fit`
* `robust_sigma`

### Web Socket Procedures

* `query_irsa_cat`
* `querydss`
* `querygsc`
* `querysimbad`
* `queryvizier`
* `read_ipac_table`
* `read_ipac_var`
* `webget`
* `write_ipac_table`

### TV Display Procedures

* `blink`
* `curs`
* `curval`
* `pixcolor`
* `sigrange`
* `tvbox`
* `tvcircle`
* `tvellipse`
* `tvlaser`
* `tvlist`
* `unzoom_xy`
* `zoom_xy`

### Miscellaneous (Non-Astronomy) Procedures

* `blkshift`
* `boost_array`
* `break_path`
* `bsort`
* `checksum32`
* `concat_dir`
* `delvarx`
* `detabify`
* `expand_tilde`
* `f_format`
* `fdecomp`
* `file_launch`
* `find_all_dir`
* `find_with_def`
* `findpro`
* `forprint`
* `get_pipe_filesize`
* `getopt`
* `getpro`
* `gettok`
* `getwrd`
* `hgrep`
* `host_to_ieee`
* `hprint`
* `ieee_to_host`
* `is_ieee_big`
* `isarray`
* `list_with_path`
* `make_2d`
* `match`
* `match2`
* `minmax`
* `mrd_skip`
* `n_bytes`
* `nint`
* `nulltrim`
* `one_arrow`
* `one_ray`
* `qget_string`
* `rdfloat`
* `read_key`
* `readcol`
* `readfmt`
* `rem_dup`
* `remchar`
* `remove`
* `repchr`
* `repstr`
* `select_w`
* `spec_dir`
* `store_array`
* `str_index`
* `strcompress2`
* `strn`
* `strnumber`
* `textclose`
* `textopen`
* `to_hex`
* `valid_num`
* `vect`
* `wherenan`
* `xdispstr`
* `zparcheck`
