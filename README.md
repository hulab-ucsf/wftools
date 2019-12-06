# wftools
Waveform extraction/processing tools for patient monitoring archive files.

The project consists of the following executables:

- `wfconvert` process patient monitoring archive files, and convert the output to binary format for furtner processing.
- `wfshow` show waveform for quick inspection
- `wfedit` manipulate waveform files (example: select specific channels and combining multiple files, etc.) (coming soon)

It takes configuration in YAML format.  Please check the YAML file in the repository for example.

The output files in .adibin format can be displayed by LabChart's free viewer:

https://www.adinstruments.com/products/labchart-reader

This version has been tested under Python 3.7 on Windows.

## To get help
```
> wfconvert -h
```

## Example: wfconvert
```
>wfconvert -d D:\data_extraction\test -o D:\data_extraction\test_output
wfconvert v0.5 - Copyright(c) HuLab@UCSF 2019
reading config file: wfconvert_config.yaml
Start processing at: 2019-03-13 10:34:04
Processing XML file: 11-30-08-000Z.xml
Processing XML file: 12-00-08-000Z.xml
Processing XML file: 12-30-08-000Z.xml
...
...
```
## Example: wfshow
```
wfshow -f D:\test\test1.adibin -s 45000 -n 100 --size=1280x720
```

## Sample wfconvert_config.yaml file
```
# output_fn_pattern, supported names are: 
# starttime, endtime, exetime, id1, id2, id3, id4, id5
# example:
# output_fn_pattern: "{id1}_{starttime}_{endtime}_{exetime}_{id2}_{id3}"
# output_fn_pattern: "{id1}_{starttime}_{endtime}_{exetime}"
output_fn_time_format_list:
  - key: "starttime"
    value: "%Y-%m-%d"
  - key: "endtime"
    value: "%Y%m%d%H%M%S"
  - key: "exetime"
    value: "%Y%m%d%H%M%S"
output_fn_pattern: "{starttime}-{id1}"
output_fn_ext: "adibin"
sampling_rate: 300
# if true, it would not fill in gaps or fix overlap issue if source file segment has gap
ignore_gap: False
# if true, it would not fill in gaps or fix overlap issue between segments (i.e. between source xml files)
ignore_gap_between_segs: False
# if true, it would print out warnings when encountering gaps or overlaps (only when ignore_gap or ignore_gap_between_segs are False)
warning_on_gaps: True
# target pattern for channels to be extracted.  Comment out
# the following lines will extract all channels
channel_pattern_list:
  - "ECG.*"
  - "PLETH"
channel_info_list:
  - label: "ECG.*"
    renameTo: "ECG"
    uom: "mV"
    rangeLow: 0.0
    rangeHigh: 1.0
    offset: 0.0
    scale: 0.1
  - label: "PLETH"
    uom: "mV"
    rangeLow: 0.0
    rangeHigh: 1.0
    offset: 0.0
    scale: 0.1
  - label: "SPO2"
    uom: "%"
    rangeLow: 0.0
    rangeHigh: 100.0
    offset: 0.0
    scale: 0.1
```

## Release Notes
v0.66
- Support file path with spaces in BedMaster STP filename

v0.5
- Fixed bug of not using key "exe_param" in config file
- Prepare for open source release
- Updated Makefile
- Updated project structure, tests, etc.
- Renamed project to "wftools"
- Renamed executable to "wfconvert"
- Added executable "wfshow"

v0.4
- Support for BedMaster STP file format (based on STPToolkit v8.2 or later)
- Support for multiple file formats by using attribute ("converter_type", either "bernoulli" or "bedmaster") in YAML config file
- Support for converter specific parameters ("converter_options") in YAML config file
- Support all features, like select channels, up/down sampling, for all converter types
- Fixed renameTo bug if "{endtime}" is used in output filename pattern

v0.3
- Added support for customizable timestamp format in YAML file (option: output_fn_time_format_list)
- Added "renameTo" in channel settings in channel_info_list in YAML file (option: "renameTo")
- Use regular expression to match channel label in channel_info_list in YAML file
- Channel label matching are now case-insensitive

v0.2
- Handle gaps and overlap
- Use "warning_on_gaps" option to show gaps or overlap in XML

v0.1
- Initial release
- Use YAML file for configuration
- Allow the use of command line options to override configuration in YAML file
