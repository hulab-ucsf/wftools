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
from pathlib import Path
import datetime
from dateutil import parser
import random
from .xml2bin_state import Xml2BinState
from .xmlconverter import XmlConverterError
from .fixsampling import fixsamplingarr
from myutil import parsetime
from myutil import dtTimestampFormat
from myutil import getOutputFilename
from binfilepy import BinFile
from binfilepy import CFWBINARY
from binfilepy import CFWBCHANNEL
from binfilepy import constant
from vitalfilepy import VitalFile
from vitalfilepy import VITALBINARY
import xml.etree.ElementTree as ET
import base64
from array import array
import numpy as np
import math
from typing import List
from typing import Dict
from typing import Any

DEFAULT_VS_LIMIT_LOW = -999999
DEFAULT_VS_LIMIT_HIGH = 999999


class XmlConverterForBedMaster:
    """
    This class provides XML conversion to Bin file
    """
    outputDir = ""
    outputFnExt = ""
    outputFnPattern = ""
    defaultSamplesPerSec = 0
    channelPatternList = None
    channelInfoList = None
    ignoreGap = False
    ignoreGapBetweenSegs = False
    warningOnGaps = False
    outputFnTimeFormatDict = None
    header = None
    headerStartDt = None
    channels = []
    name2Channel = {}
    outputFileSet = set()
    outputFileList = []

    def __init__(self, outputDir: str = "", outputFnPattern: str = "", outputFnExt: str = "", defaultSamplesPerSec: int = 0,
                 channelPatternList: List = None, channelInfoList: List = None,
                 ignoreGap: bool = False, ignoreGapBetweenSegs: bool = False, warningOnGaps: bool = False,
                 outputFnTimeFormatDict: Dict = None):
        self.outputDir = outputDir
        self.outputFnPattern = outputFnPattern
        self.outputFnExt = outputFnExt
        self.defaultSamplesPerSec = defaultSamplesPerSec
        self.channelPatternList = channelPatternList
        self.channelInfoList = channelInfoList
        self.ignoreGap = ignoreGap
        self.ignoreGapBetweenSegs = ignoreGapBetweenSegs
        self.warningOnGaps = warningOnGaps
        self.outputFnTimeFormatDict = outputFnTimeFormatDict

    def clearState(self):
        self.header = None
        self.headerStartDt = None
        self.channels = []
        self.name2Channel = {}

    def clearOutputFileList(self):
        self.outputFileSet = set()
        self.outputFileList = []

    def setDefaultSamplesPerSec(self, defaultSamplesPerSec: int):
        self.defaultSamplesPerSec = defaultSamplesPerSec

    def setChannelPatternList(self, channelPatternList: List):
        self.channelPatternList = channelPatternList

    def inChannelPatternList(self, label: str):
        if (self.channelPatternList is None) or (len(self.channelPatternList) == 0):
            return True
        for p in self.channelPatternList:
            if p.match(label) is not None:
                return True
        return False

    def getChannelInfo(self, label: str):
        if self.channelInfoList is not None:
            for cSettingInfo in self.channelInfoList:
                labelPattern = cSettingInfo.get("labelPattern")
                if labelPattern.match(label) is not None:
                    return cSettingInfo
        return None

    def renameOutputFnWithEndtime(self, numSamples: int, tagsDict: Dict, x: Xml2BinState, filename: str):
        # rename the file that we just closed (if filename pattern has {endtime})
        if "{endtime}" in self.outputFnPattern:
            fmt = self.outputFnTimeFormatDict.get("starttime", None) if (self.outputFnTimeFormatDict is not None) else None
            tagsDict["starttime"] = dtTimestampFormat(self.headerStartDt, fmt)
            fmt = self.outputFnTimeFormatDict.get("exetime", None) if (self.outputFnTimeFormatDict is not None) else None
            tagsDict["exetime"] = dtTimestampFormat(x.timestampTm, fmt)
            endDt = self.headerStartDt + datetime.timedelta(seconds=int(numSamples / self.defaultSamplesPerSec))
            fmt = self.outputFnTimeFormatDict.get("endtime", None) if (self.outputFnTimeFormatDict is not None) else None
            tagsDict["endtime"] = dtTimestampFormat(endDt, fmt)
            if "tempendtime" in filename:
                original_filename = filename
                filename = getOutputFilename(self.outputDir, self.outputFnPattern, tagsDict, self.outputFnExt)
                if original_filename != filename:
                    os.rename(original_filename, filename)
                x.lastBinFilename = filename
            # end-if
        # end-if

    # return total number of samples written
    def convert(self, xmlFile: str, tagsDict: Dict, x: Xml2BinState, print_processing_fn: bool = False):
        """
        convert to BIN file from XML
        """
        binFileOut = None
        filename = ""
        chanLabel = []
        tempChanInfo = []
        tempChanLabel = []
        tempChanLabel2Index = {}
        startWaveTm = datetime.datetime.min
        numSamples = 0
        totalNumSamplesWritten = 0
        firstBinFile = True
        firstMeasurement = True

        # array of parName, startTm, vitalFileOut, filename
        startVitalTm = datetime.datetime.min
        vitalFileInfoArr = []
        vitalParName2Info = {}

        # progress
        if (x.lastVitalFileInfoArr is not None) and (len(x.lastVitalFileInfoArr) > 0):
            for vinfo in x.lastVitalFileInfoArr:
                vs_parameter = vinfo["par"]
                vs_startTm = vinfo["startTm"]
                vitalFilename = vinfo["filename"]
                vitalFileOut = VitalFile(vitalFilename, "r+")
                vitalFileOut.open()
                vitalFileInfo = {"par": vs_parameter, "startTm": startVitalTm, "vitalFileOut": vitalFileOut, "filename": vitalFilename}
                vitalFileInfoArr.append(vitalFileInfo)
                vitalParName2Info[vs_parameter] = vitalFileInfo

        xml_unit = ""
        xml_bed = ""

        infilePath = Path(xmlFile)
        if not infilePath.exists():
            raise XmlConverterError("Cannot open file: {0}".format(xmlFile))
        else:
            if len(x.lastBinFilename) > 0:
                binFileOut = BinFile(x.lastBinFilename, "r+")
                binFileOut.open()
                binFileOut.readHeader()
                filename = x.lastBinFilename
                firstBinFile = False
                # copy BIN file's channel info
                chanLabel = []
                tempChanLabel = []
                chanLabel2Index = {}
                tempChanLabel2Index = {}
                idx = 0
                for c in binFileOut.channels:
                    chanLabel.append(c.Title)
                    chanLabel2Index[c.Title] = idx
                    tempChanLabel.append(c.Title)
                    tempChanLabel2Index[c.Title] = idx
                    idx += 1
                self.header = binFileOut.header
                second = int(math.floor(self.header.Second))
                microsecond = int((self.header.Second - math.floor(self.header.Second)) * 100)
                self.headerStartDt = datetime.datetime(self.header.Year, self.header.Month, self.header.Day,
                                                       self.header.Hour, self.header.Minute, second, microsecond)
                numSamples = self.header.SamplesPerChannel
            # end-if len(x.lastBinFilename)
            if print_processing_fn:
                print("Processing XML file: {0}".format(infilePath.name))
            tree = ET.parse(xmlFile)
            root = tree.getroot()
            # print("root tag = {0}".format(root.tag))
            if root.tag == "BedMasterEx":
                for child1 in root:
                    if child1.tag == "FileInfo":
                        for child3 in child1:
                            if child3.tag == "Unit":
                                xml_unit = child3.text
                            elif child3.tag == "Bed":
                                xml_bed = child3.text
                    if child1.tag == "Segment":
                        for child3 in child1:
                            if child3.tag == "Waveforms":
                                collectionTime, collectionTimeUTC = self.processWaveforms(child3)
                                collectionTimeDt = parsetime(collectionTime)
                                tempChanInfo = []
                                tempChanLabel = []
                                tempChanLabel2Index = {}
                                if self.header is None:
                                    self.headerStartDt = collectionTimeDt
                                    self.header = CFWBINARY()
                                    self.header.setValue(1.0 / self.defaultSamplesPerSec, collectionTimeDt.year, collectionTimeDt.month,
                                                         collectionTimeDt.day, collectionTimeDt.hour, collectionTimeDt.minute, collectionTimeDt.second, 0, 0)
                                # print(collectionTime, collectionTimeUTC)
                                idx = 0
                                for child4 in child3:
                                    if (child4.tag == "WaveformData"):
                                        ID, channel, hz, points, uom, wave = self.processWaveformData(child4)
                                        if self.inChannelPatternList(channel):
                                            # print(channel, wave, points, pointsBytes, min_, max_, offset, gain, hz)
                                            wavedata = self.decodeWave(wave)
                                            # print(wavedata)
                                            if self.defaultSamplesPerSec != hz:
                                                wavedata = fixsamplingarr(wavedata, hz, self.defaultSamplesPerSec)
                                            tempChanInfo.append({"label": channel, "data": wavedata, "points": points, "hz": hz, "uom": uom})
                                            tempChanLabel.append(channel)
                                            tempChanLabel2Index[channel] = idx
                                            idx += 1
                                # end-for child4
                                if (firstBinFile is True) or (len(tempChanLabel) > 0 and self.channelChanged(chanLabel, tempChanLabel)):
                                    if firstBinFile is False:
                                        binFileOut.close()
                                        binFileOut = None
                                        # rename the file that we just closed (if filename pattern has {endtime})
                                        self.renameOutputFnWithEndtime(numSamples, tagsDict, x, filename)
                                        if not (x.lastBinFilename in self.outputFileSet):
                                            self.outputFileSet.add(x.lastBinFilename)
                                            self.outputFileList.append(x.lastBinFilename)
                                        # reset new headerStartDt
                                        self.headerStartDt = collectionTimeDt
                                        self.header.setValue(1.0 / self.defaultSamplesPerSec, collectionTimeDt.year, collectionTimeDt.month,
                                                             collectionTimeDt.day, collectionTimeDt.hour, collectionTimeDt.minute, collectionTimeDt.second, 0, 0)
                                    firstBinFile = False
                                    self.header.NChannels = len(tempChanInfo)
                                    fmt = self.outputFnTimeFormatDict.get("starttime", None) if (self.outputFnTimeFormatDict is not None) else None
                                    tagsDict["starttime"] = dtTimestampFormat(self.headerStartDt, fmt)
                                    fmt = self.outputFnTimeFormatDict.get("exetime", None) if (self.outputFnTimeFormatDict is not None) else None
                                    tagsDict["exetime"] = dtTimestampFormat(x.timestampTm, fmt)
                                    # we do not know the end at this point
                                    tagsDict["endtime"] = "tempendtime" + str(random.randint(10000, 100000))
                                    filename = getOutputFilename(self.outputDir, self.outputFnPattern, tagsDict, self.outputFnExt)
                                    x.lastBinFilename = filename
                                    binFileOut = BinFile(filename, "w")
                                    binFileOut.open()
                                    binFileOut.setHeader(self.header)
                                    chanData = []
                                    chanLabel = []
                                    for cinfo in tempChanInfo:
                                        label = cinfo["label"]
                                        cSettingInfo = self.getChannelInfo(label)
                                        if cSettingInfo is not None:
                                            uom = cSettingInfo.get("uom", cinfo.get("uom", ""))
                                            rangeLow = cSettingInfo.get("rangeLow", 0)
                                            rangeHigh = cSettingInfo.get("rangeHigh", 100)
                                            offset = cSettingInfo.get("offset", 0)
                                            scale = cSettingInfo.get("scale", 1)
                                        else:
                                            uom = cinfo.get("uom", "")
                                            rangeLow = 0
                                            rangeHigh = 100
                                            offset = 0
                                            scale = 1
                                        channel = CFWBCHANNEL()
                                        channel.setValue(label, uom, scale, offset, rangeLow, rangeHigh)
                                        binFileOut.addChannel(channel)
                                        chanData.append(cinfo["data"])
                                        chanLabel.append(label)
                                    binFileOut.writeHeader()
                                    firstMeasurement = False
                                    numSamples = binFileOut.writeChannelData(chanData)
                                    totalNumSamplesWritten += numSamples
                                    binFileOut.updateSamplesPerChannel(numSamples, True)
                                elif len(tempChanInfo) > 0:
                                    chanData = []
                                    chanLabel = []
                                    for cinfo in tempChanInfo:
                                        label = cinfo["label"]
                                        chanData.append(cinfo["data"])
                                        chanLabel.append(label)
                                    # gap handling
                                    endDt = self.headerStartDt + datetime.timedelta(seconds=int(numSamples / self.defaultSamplesPerSec))
                                    gap = collectionTimeDt - endDt
                                    actualGapInSec = int(gap.total_seconds())
                                    if (not self.ignoreGapBetweenSegs) and firstMeasurement:
                                        gapInSec = actualGapInSec
                                    elif not self.ignoreGap:
                                        gapInSec = actualGapInSec
                                    else:
                                        gapInSec = 0
                                    if self.warningOnGaps and (gapInSec != 0):
                                        print("Waveforms CollectionTime: {0} shows gap (or overlap) = {1} secs".format(collectionTime, gapInSec))
                                    firstMeasurement = False
                                    numSamplesWritten = binFileOut.writeChannelData(chanData, self.defaultSamplesPerSec, gapInSec)
                                    totalNumSamplesWritten += numSamplesWritten
                                    numSamples = numSamplesWritten + binFileOut.header.SamplesPerChannel
                                    binFileOut.header.SamplesPerChannel = numSamples
                                    self.header.SamplesPerChannel = numSamples
                                    binFileOut.updateSamplesPerChannel(numSamples, True)
                                # end-if firstBinFile
                            # end-if "measurement"
                            if child3.tag == "VitalSigns":
                                collectionTime, collectionTimeUTC = self.processVitalSigns(child3)
                                collectionTimeDt = parsetime(collectionTime)
                                vs_parameter = ""
                                vs_time = ""
                                vs_value = ""
                                vs_uom = ""
                                vs_alarmLimitLow = ""
                                vs_alarmLimitHigh = ""
                                for child4 in child3:
                                    if (child4.tag == "VitalSign"):
                                        vs_parameter, vs_time, vs_value, vs_uom, vs_alarmLimitLow, vs_alarmLimitHigh = self.processVitalSign(child4)
                                        # print(vs_parameter)
                                    if (vs_parameter is not None) and len(vs_parameter) > 0:
                                        vitalFileInfo = None
                                        if vs_parameter in vitalParName2Info:
                                            vitalFileInfo = vitalParName2Info.get(vs_parameter)
                                        else:
                                            vs_time_dt = parsetime(vs_time)
                                            fmt = self.outputFnTimeFormatDict.get("starttime", None) if (self.outputFnTimeFormatDict is not None) else None
                                            tagsDict["starttime"] = dtTimestampFormat(vs_time_dt, fmt)
                                            fmt = self.outputFnTimeFormatDict.get("exetime", None) if (self.outputFnTimeFormatDict is not None) else None
                                            tagsDict["exetime"] = dtTimestampFormat(x.timestampTm, fmt)
                                            # we do not know the end at this point
                                            tagsDict["endtime"] = "0000"  # "tempendtime" + str(random.randint(10000, 100000))
                                            # array of parName, startTm, vitalFileOut, filename
                                            vitalFilename = getOutputFilename(self.outputDir, self.outputFnPattern + "_" + vs_parameter, tagsDict, "vital")
                                            vitalFileOut = VitalFile(vitalFilename, "w")
                                            vitalFileOut.open()
                                            if startVitalTm == datetime.datetime.min:
                                                startVitalTm = vs_time_dt
                                            vitalFileInfo = {"par": vs_parameter, "startTm": startVitalTm, "vitalFileOut": vitalFileOut, "filename": vitalFilename}
                                            vitalFileInfoArr.append(vitalFileInfo)
                                            vitalParName2Info[vs_parameter] = vitalFileInfo
                                            vs_header = VITALBINARY(vs_parameter, vs_uom, xml_unit, xml_bed, collectionTimeDt.year, collectionTimeDt.month, startVitalTm.day, startVitalTm.hour, startVitalTm.minute, startVitalTm.second)
                                            vitalFileOut.setHeader(vs_header)
                                            vitalFileOut.writeHeader()
                                        if vitalFileInfo is not None:
                                            vs_value_num = DEFAULT_VS_LIMIT_LOW
                                            try:
                                                vs_value_num = float(vs_value)
                                            except:
                                                pass
                                            vs_time_dt = parsetime(vs_time)
                                            vs_offset_num = (vs_time_dt - vitalFileInfo["startTm"]).total_seconds()
                                            vs_low_num = DEFAULT_VS_LIMIT_LOW
                                            try:
                                                vs_low_num = float(vs_alarmLimitLow)
                                            except:
                                                pass
                                            vs_high_num = DEFAULT_VS_LIMIT_HIGH
                                            try:
                                                vs_high_num = float(vs_alarmLimitHigh)
                                            except:
                                                pass
                                            vitalFileOut = vitalFileInfo["vitalFileOut"]
                                            vitalFileOut.writeVitalData(vs_value_num, vs_offset_num, vs_low_num, vs_high_num)
                # end-for child1
            # end-if root
        # end-if
        if binFileOut is not None:
            binFileOut.close()
            binFileOut = None
            # rename the file that we just closed (if filename pattern has {endtime})
            self.renameOutputFnWithEndtime(numSamples, tagsDict, x, filename)
            if not (x.lastBinFilename in self.outputFileSet):
                self.outputFileSet.add(x.lastBinFilename)
                self.outputFileList.append(x.lastBinFilename)
        for vf in vitalFileInfoArr:
            vitalFileOut = vf["vitalFileOut"]
            if vitalFileOut is not None:
                vitalFileOut.close()
                vf["vitalFileOut"] = None
            x.addOrUpdateLastVitalFileInfo(vf["par"], vf["startTm"], vf["filename"])

        return totalNumSamplesWritten

    def renameChannels(self, print_rename_details: bool = False):
        # progress
        # need to support renaming of channel label
        # also, support use of Regex for channel label, and renameTo expression
        hasRenameTo = False
        for cSettingInfo in self.channelInfoList:
            renameTo = cSettingInfo.get("renameTo", "")
            if len(renameTo) > 0:
                hasRenameTo = True
                break
            # end-if
        # end-for
        if hasRenameTo:
            print("Renaming channels in output files...")
            numFilesChanged = 0
            for fn in self.outputFileList:
                print("processing {0}...".format(Path(fn).name))
                updatedFile = False
                with BinFile(fn, "r+") as f:
                    f.readHeader()
                    for c in f.channels:
                        cSettingInfo = self.getChannelInfo(c.Title)
                        if cSettingInfo is not None:
                            newLabel = cSettingInfo.get("renameTo", "")
                            if len(newLabel) > 0:
                                filename = Path(fn).name
                                oldLabel = c.Title
                                c.Title = newLabel
                                updatedFile = True
                                if print_rename_details:
                                    print("{0} file: {1} -> {2}".format(filename, oldLabel, newLabel))
                    if updatedFile:
                        f.writeHeader()
                        numFilesChanged += 1
                # end-with
            # end-for
            if numFilesChanged == 0:
                if print_rename_details:
                    print("No output files's channel labels need to be changed.")
        # end-if

    def processWaveforms(self, e: object):
        collectionTime = e.attrib.get("CollectionTime", "")
        collectionTimeUTC = e.attrib.get("CollectionTimeUTC", "")
        return collectionTime, collectionTimeUTC

    def processWaveformData(self, e: object):
        ID = int(e.attrib.get("ID", ""))
        channel = e.attrib.get("Label", "")
        hz = float(e.attrib.get("SampleRate", ""))
        points = int(e.attrib.get("Samples", ""))
        uom = e.attrib.get("UOM", "")
        wave = e.text
        return ID, channel, hz, points, uom, wave

    def decodeWave(self, x: str):
        val_list = x.split(",")
        a = array("h", (0,) * len(val_list))
        i = 0
        for v in val_list:
            a[i] = int(v)
            i += 1
        return a

    def processVitalSigns(self, e: object):
        collectionTime = e.attrib.get("CollectionTime", "")
        collectionTimeUTC = e.attrib.get("CollectionTimeUTC", "")
        return collectionTime, collectionTimeUTC

    def processVitalSign(self, e: object):
        vs_parameter = ""
        vs_time = ""
        vs_value = ""
        vs_uom = ""
        vs_alarmLimitLow = ""
        vs_alarmLimitHigh = ""
        for child in e:
            if child.tag == "Parameter":
                if child.text is not None:
                    vs_parameter = child.text
            elif child.tag == "Time":
                if child.text is not None:
                    vs_time = child.text
            elif child.tag == "Value":
                if child.text is not None:
                    vs_value = child.text
                vs_uom = child.attrib.get("UOM", "")
            elif child.tag == "AlarmLimitLow":
                if child.text is not None:
                    vs_alarmLimitLow = child.text
            elif child.tag == "AlarmLimitHigh":
                if child.text is not None:
                    vs_alarmLimitHigh = child.text
        return vs_parameter, vs_time, vs_value, vs_uom, vs_alarmLimitLow, vs_alarmLimitHigh

    def moveTempChanLabel(self, chanLabelArr: List[str], tempChanLabelArr: List[str]):
        chanLabelArr.clear()
        for l in tempChanLabelArr:
            chanLabelArr.append(l)

    def channelChanged(self, chanLabelArr: List[str], tempChanLabelArr: List[str]):
        changed = False
        if len(chanLabelArr) != len(tempChanLabelArr):
            changed = True
        else:
            for l, t in zip(chanLabelArr, tempChanLabelArr):
                if l != t:
                    changed = True
                    break
        return changed
