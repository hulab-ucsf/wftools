"""
MIT License

Copyright (c) 2019 UCSF Hu Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os
import sys
import re
import yaml
import argparse
from datetime import datetime
from datetime import timedelta
from pathlib import Path
from myutil import dtFormat
from myutil import dtTimestampFormat
from myutil import elapsedFormat
from myutil import parsetime
from xmlconvert import Xml2BinState
from xmlconvert import XmlConverterForGE
from xmlconvert import XmlConverterForBedMaster
from typing import Dict

g_version = "0.63"
g_exename = "wfconvert"
default_config_fn = "{0}_config.yaml".format(g_exename)
default_sampling_rate = 240
default_fn_ext = "adibin"


def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input File")
    parser.add_argument("-d", "--dir", help="Input Directory")
    parser.add_argument("-o", "--output_dir", help="Output Directory", required=True)
    parser.add_argument("-c", "--config_file", help="configuration file")
    parser.add_argument("-s", "--sampling_rate", help="Target sampling rate")
    parser.add_argument("-p", "--channel_patterns", help="comma separated regex pattern for channels")
    parser.add_argument("--stime", help="specify start time, format: \"1/1/2019 8:00:00 AM\"")
    parser.add_argument("--etime", help="specify start time, format: \"1/1/2019 12:00:00 PM\"")
    # use store_const, instead of store_true, so that default value is None (instead of False)
    parser.add_argument("--ignore_gap", help="ignore gap or overlap within a source file", action="store_const", const=False)
    parser.add_argument("--ignore_gap_between_segs", help="ignore gap or overlap between segments (or xml files)", action="store_const", const=False)
    parser.add_argument("--warning_on_gaps", help="show warning when encountering gaps or overlaps (if gaps are not ignored)", action="store_const", const=False)
    parser.add_argument("--output_fn_pattern", help="output filename pattern")
    parser.add_argument("--output_fn_ext", help="output file extention, e.g. adibin, bin")
    parser.add_argument("--id1", help="id1 tag in output file")
    parser.add_argument("--id2", help="id2 tag in output file")
    parser.add_argument("--id3", help="id3 tag in output file")
    parser.add_argument("--id4", help="id4 tag in output file")
    parser.add_argument("--id5", help="id5 tag in output file")
    return parser.parse_args()


def printOptions():
    print("Options used:")
    if (g_file is not None) and len(g_file) > 0:
        print("\tinput file: {0}".format(g_file))
    if (g_dir is not None) and len(g_dir) > 0:
        print("\tinput dir: {0}".format(g_dir))
    print("\toutput dir: {0}".format(g_output_dir))
    print("\tsampling rate: {0}".format(g_sampling_rate))
    print("\tignore gap: {0}".format(g_ignore_gap))
    print("\tignore_gap_between_segs: {0}".format(g_ignore_gap_between_segs))
    print("\twarning_on_gaps: {0}".format(g_warning_on_gaps))
    if (g_channel_patterns is not None) and len(g_channel_patterns) > 0:
        print("\tchannel_patterns: {0}".format(",".join(g_channel_patterns)))
    else:
        print("\tchannel_patterns: *")


def getTempDir():
    exeDir = Path(sys.executable).parent
    defaultTempDir = Path(exeDir).joinpath("temp_dir")
    tempDir = g_converter_options.get("temp_dir", defaultTempDir)
    return tempDir


def execute_ext(ext_exe: str, srcFn: str, cmdParam: str, xmlOutputFn: str, startSegment: int, endSegment: int, stime: datetime, etime: datetime):
    tempDir = getTempDir()
    if not os.path.exists(tempDir):
        os.mkdir(tempDir)
    xmlOutputFullFn = Path(tempDir).joinpath(xmlOutputFn)
    cmd = "{0} {1} -o {2} {3}".format(ext_exe, srcFn, xmlOutputFullFn, cmdParam)
    if (stime is not None) and (etime is not None):
        cmd = cmd + " -stime \"{0}\"".format(stime.strftime("%m/%d/%Y %I:%M:%S %p"))
        cmd = cmd + " -etime \"{0}\"".format(etime.strftime("%m/%d/%Y %I:%M:%S %p"))
    else:
        if startSegment >= 0:
            cmd = cmd + " -s {0}".format(int(startSegment))
        if endSegment >= 0:
            cmd = cmd + " -e {0}".format(int(endSegment))
    print("execute: {0}".format(cmd))
    os.system(cmd)
    return xmlOutputFullFn


def runApp(flow: str, srcFile: str, srcDir: str, dstDir: str):
    timestampTm = datetime.now()
    tagsDict = {"id1": g_id1, "id2": g_id2, "id3": g_id3, "id4": g_id4, "id5": g_id5}
    xmlconverter = None

    if flow == "file":
        xml2BinState = Xml2BinState()
        xml2BinState.setTimestampTm(timestampTm)
        if g_converter_type == "bernoulli":
            xmlconverter = XmlConverterForGE(dstDir, g_output_fn_pattern, g_output_fn_ext,
                                             int(g_sampling_rate), g_channel_pattern_list, g_channel_info_list,
                                             g_ignore_gap, g_ignore_gap_between_segs, g_warning_on_gaps,
                                             g_output_fn_time_format_dict)
        elif g_converter_type == "bedmaster":
            xmlconverter = XmlConverterForBedMaster(dstDir, g_output_fn_pattern, g_output_fn_ext,
                                                    int(g_sampling_rate), g_channel_pattern_list, g_channel_info_list,
                                                    g_ignore_gap, g_ignore_gap_between_segs, g_warning_on_gaps,
                                                    g_output_fn_time_format_dict)
        xmlconverter.convert(srcFile, tagsDict, xml2BinState, print_processing_fn=True)
        xmlconverter.renameChannels(print_rename_details=True)
        return 0
    elif flow == "dir":
        numFilesProcessed = 0
        xml2BinState = Xml2BinState()
        xml2BinState.setTimestampTm(timestampTm)
        if g_converter_type == "bernoulli":
            xmlconverter = XmlConverterForGE(dstDir, g_output_fn_pattern, g_output_fn_ext,
                                             int(g_sampling_rate), g_channel_pattern_list, g_channel_info_list,
                                             g_ignore_gap, g_ignore_gap_between_segs, g_warning_on_gaps,
                                             g_output_fn_time_format_dict)
        elif g_converter_type == "bedmaster":
            xmlconverter = XmlConverterForBedMaster(dstDir, g_output_fn_pattern, g_output_fn_ext,
                                                    int(g_sampling_rate), g_channel_pattern_list, g_channel_info_list,
                                                    g_ignore_gap, g_ignore_gap_between_segs, g_warning_on_gaps,
                                                    g_output_fn_time_format_dict)
        for file in sorted(os.listdir(srcDir)):
            # debug
            fp = os.path.join(srcDir, file)
            ext = os.path.splitext(fp)[1] if len(os.path.splitext(fp)) == 2 else ""
            if ext == ".xml":
                xmlconverter.clearState()
                xmlconverter.convert(fp, tagsDict, xml2BinState, print_processing_fn=True)
                numFilesProcessed += 1
        xmlconverter.renameChannels(print_rename_details=True)
        print("Number of XML files processed = {0}".format(numFilesProcessed))
        return 0
    elif flow == "stp":
        # need to add date-range option for processing
        ext_exe = str(g_converter_options.get("stptools_exe", ""))
        ext_param = str(g_converter_options.get("stptools_param", ""))
        keep_temp_file = bool(g_converter_options.get("keep_temp_files", False))
        numSegmentsPerBatch = int(g_converter_options.get("num_segments_per_batch", 500))
        numSegmentsProcessed = 0
        xml2BinState = Xml2BinState()
        xml2BinState.setTimestampTm(timestampTm)
        if g_converter_type == "bedmaster":
            xmlconverter = XmlConverterForBedMaster(dstDir, g_output_fn_pattern, g_output_fn_ext,
                                                    int(g_sampling_rate), g_channel_pattern_list, g_channel_info_list,
                                                    g_ignore_gap, g_ignore_gap_between_segs, g_warning_on_gaps,
                                                    g_output_fn_time_format_dict)
        startSegment = 0
        endSegment = startSegment + numSegmentsPerBatch - 1
        if (ext_exe is None) or (len(ext_exe) == 0):
            print("ERROR: stptools_exe option not set!")
            return 1
        hasSegments = True
        useTime = (g_stime is not None) and (g_etime is not None)
        time_step = timedelta(days=1)
        current_stime = g_stime
        current_etime = g_stime
        if useTime:
            current_etime = g_stime + time_step
            if current_etime > g_etime:
                current_etime = g_etime
        while ((not useTime) and hasSegments) or (useTime and (current_stime < g_etime) and (current_stime < current_etime)):
            basefn = Path(srcFile).stem
            xmlOutputFn = "{0}_{1}_{2}_wf.xml".format(basefn, dtTimestampFormat(timestampTm), startSegment)
            xmlOutputFullFn = execute_ext(ext_exe, srcFile, ext_param, xmlOutputFn, startSegment, endSegment, current_stime, current_etime)
            xmlProcesingStarttime = datetime.now()
            if Path(xmlOutputFullFn).exists():
                if not useTime:
                    print("Processing XML from segment {0} to segment {1}...".format(startSegment, endSegment))
                else:
                    print("Processing XML from {0} to {1}...".format(current_stime.strftime("%m/%d/%Y %I:%M:%S %p"), current_etime.strftime("%m/%d/%Y %I:%M:%S %p")))
                xmlconverter.clearState()
                numSamplesWritten = xmlconverter.convert(xmlOutputFullFn, tagsDict, xml2BinState, print_processing_fn=True)
                xmlProcessingEndtime = datetime.now()
                xmlProcessingElapsedtime = xmlProcessingEndtime - xmlProcesingStarttime
                if numSamplesWritten == 0:
                    print("No more data segments for processing.")
                print("Processing XML takes {0}".format(elapsedFormat(xmlProcessingElapsedtime.total_seconds(), totalSecondsOnly=True)))
                if not keep_temp_file:
                    os.remove(xmlOutputFullFn)
            startSegment = endSegment + 1
            endSegment = startSegment + numSegmentsPerBatch - 1
            numSegmentsProcessed += numSegmentsPerBatch
            hasSegments = (numSamplesWritten > 0)
            if useTime:
                current_stime = current_etime
                current_etime = current_stime + time_step
                if current_etime > g_etime:
                    current_etime = g_etime
        xmlconverter.renameChannels(print_rename_details=True)
        print("Done")


print("{0} v{1} - Copyright(c) HuLab@UCSF 2019".format(g_exename, g_version))
args = getArgs()
g_converter_type = ""
g_converter_options = {}
g_config_file = default_config_fn
exeDir = Path(sys.executable).parent
g_config_file = Path(exeDir).joinpath(default_config_fn)
if not Path(g_config_file).exists():
    g_config_file = Path.cwd().joinpath(default_config_fn)
if not Path(g_config_file).exists():
    g_config_file = Path.cwd().joinpath("src", default_config_fn)
g_config_file = str(g_config_file)
if args.config_file is not None:
    g_config_file = args.config_file
g_file = args.file
g_dir = args.dir
g_output_dir = args.output_dir
g_output_fn_pattern = None
g_output_fn_ext = default_fn_ext
g_sampling_rate = default_sampling_rate
g_channel_patterns = None
g_channel_pattern_list = None
g_channel_info_list = None
g_ignore_gap = False
g_ignore_gap_between_segs = False
g_warning_on_gaps = False
g_stime = None
g_etime = None
if args.stime is not None:
    g_stime = parsetime(args.stime)
if args.etime is not None:
    g_etime = parsetime(args.etime)
g_id1 = args.id1 if args.id1 is not None else ""
g_id2 = args.id2 if args.id2 is not None else ""
g_id3 = args.id3 if args.id3 is not None else ""
g_id4 = args.id4 if args.id4 is not None else ""
g_id5 = args.id5 if args.id5 is not None else ""
configData = None
if (g_config_file is not None) and (len(g_config_file) > 0) and os.path.exists(g_config_file):
    print("Reading config file: {0}".format(g_config_file))
    with open(g_config_file, 'r') as stream:
        configData = yaml.load(stream)
if configData is not None:
    if configData.get("converter_type") is not None:
        g_converter_type = str(configData.get("converter_type")).lower()
    if configData.get("converter_options") is not None:
        g_converter_options = {}
        if isinstance(configData.get("converter_options"), list):
            for c in configData.get("converter_options"):
                g_converter_options[c.get("key")] = c.get("value")
    if configData.get("output_fn_time_format_list") is not None:
        if isinstance(configData.get("output_fn_time_format_list"), list):
            g_output_fn_time_format_dict = {}
            for t in configData.get("output_fn_time_format_list"):
                g_output_fn_time_format_dict[t["key"].lower()] = t["value"]
    if configData.get("output_fn_pattern") is not None:
        g_output_fn_pattern = configData.get("output_fn_pattern")
    if configData.get("output_fn_ext") is not None:
        g_output_fn_ext = configData.get("output_fn_ext")
    if configData.get("sampling_rate") is not None:
        g_sampling_rate = float(configData.get("sampling_rate"))
    if configData.get("ignore_gap") is not None:
        g_ignore_gap = bool(configData.get("ignore_gap"))
    if configData.get("ignore_gap_between_segs") is not None:
        g_ignore_gap_between_segs = bool(configData.get("ignore_gap_between_segs"))
    if configData.get("warning_on_gaps") is not None:
        g_warning_on_gaps = bool(configData.get("warning_on_gaps"))
    if configData.get("channel_pattern_list") is not None:
        g_channel_patterns = configData.get("channel_pattern_list")
    if configData.get("channel_info_list") is not None:
        if isinstance(configData.get("channel_info_list"), list):
            g_channel_info_list = []
            for c in configData.get("channel_info_list"):
                g_channel_info_list.append(c)
            for cinfo in g_channel_info_list:
                if len(cinfo.get("label", "")) > 0:
                    cinfo["labelPattern"] = re.compile(cinfo.get("label", ""), flags=re.IGNORECASE)
if args.output_fn_pattern is not None:
    g_output_fn_pattern = args.output_fn_pattern
if args.output_fn_ext is not None:
    output_fn_ext = args.output_fn_ext
if args.sampling_rate is not None:
    g_sampling_rate = float(args.sampling_rate)
if args.channel_patterns is not None:
    g_channel_patterns = args.channel_patterns.split(",")
if (g_channel_patterns is not None) and len(g_channel_patterns) > 0:
    g_channel_pattern_list = []
    for p in g_channel_patterns:
        g_channel_pattern_list.append(re.compile(p, flags=re.IGNORECASE))
if args.ignore_gap is not None:
    g_ignore_gap = bool(args.ignore_gap)
if args.ignore_gap_between_segs is not None:
    g_ignore_gap_between_segs = bool(args.ignore_gap_between_segs)
if args.warning_on_gaps is not None:
    g_warning_on_gaps = bool(args.warning_on_gaps)

flow = "unknown"
if g_file is None and g_dir is not None:
    flow = "dir"
elif g_file is not None and g_dir is None:
    flow = "file"

if g_converter_type == "bedmaster" and flow == "file":
    ext = os.path.splitext(g_file)[1] if len(os.path.splitext(g_file)) == 2 else ""
    if ext.lower() == ".stp":
        flow = "stp"

valid = True
if (not g_file) and (not g_dir):
    print("You must specify -f or -d option!!")
    valid = False
if valid and (not g_output_dir):
    print("You must specify -o option!!")
    valid = False
if valid and (not (os.path.exists(g_output_dir))):
    print("Output directory is not accessible!!")
    valid = False


if valid:
    starttime = datetime.now()
    print("Start processing at: {0}".format(dtFormat(starttime)))
    printOptions()
    result = runApp(flow, g_file, g_dir, g_output_dir)
    if result is not None and result != 0:
        print("Error during processing!")
    endtime = datetime.now()
    print("Finished processing at: {0}".format(dtFormat(endtime)))
    elapsedtime = endtime - starttime
    print("Total elapsed time: {0}".format(
        elapsedFormat(elapsedtime.total_seconds())))
