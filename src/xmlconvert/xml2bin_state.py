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
import datetime
from myutil import parsetime


class Xml2BinState:
    lastBinFilename = ""
    # array of parName, startTm, filename
    lastVitalFileInfoArr = []
    vitalParName2LastVitalFileInfo = {}
    xmlStartTm = None
    xmlEndTm = None
    timestampTm = ""

    def __init__(self):
        self.lastBinFilename = ""
        self.lastVitalFileInfoArr = []
        self.vitalParName2LastVitalFileInfo = {}
        self.xmlStartTm = None
        self.xmlEndTm = None
        self.timestampTm = None

    def freeXmlBinState(self):
        self.lastBinFilename = ""
        self.lastVitalFileInfoArr = []
        self.vitalParName2LastVitalFileInfo = {}
        self.xmlStartTm = None
        self.xmlEndTm = None
        self.timestampTm = None

    def setLastBinFilename(self, fn: str):
        self.lastBinFilename = fn

    def addOrUpdateLastVitalFileInfo(self, par: str, startTm: datetime.datetime, filename: str):
        if not (par in self.vitalParName2LastVitalFileInfo):
            vitalFileInfo = {"par": par, "startTm": startTm, "filename": filename}
            self.lastVitalFileInfoArr.append(vitalFileInfo)
            self.vitalParName2LastVitalFileInfo[par] = vitalFileInfo
        else:
            vitalFileInfo = self.vitalParName2LastVitalFileInfo[par]
            vitalFileInfo["startTm"] = startTm
            vitalFileInfo["filename"] = filename

    def isParInLastVitalFileInfo(self, par: str):
        return par in self.vitalParName2LastVitalFileInfo

    def getLastVitalFileInfo(self, par: str):
        if par in self.vitalParName2LastVitalFileInfo:
            return self.vitalParName2LastVitalFileInfo[par]
        return None

    def setTimestampTm(self, timestamp: object):
        if isinstance(timestamp, str):
            self.timestampTm = parsetime(timestamp)
        elif isinstance(timestamp, datetime.datetime):
            self.timestampTm = timestamp

    def setXmlStartEndTm(self, startTimeStr: str, endTimeStr: str):
        self.xmlStartTm = parsetime(startTimeStr)
        self.xmlEndTm = parsetime(endTimeStr)
