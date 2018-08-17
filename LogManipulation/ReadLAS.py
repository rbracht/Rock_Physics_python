#!/usr/bin/python

import Tkinter as tk
import MyModule as mm
import tkFileDialog
import sys

class PopulateBunchOfLogsFromLasFile:
	
	def __init__( self, parent = tk.Tk(), aBunchOfLogs = mm.BunchOfLogs() ):
		"""
		
		Creates a 'BunchOfLogs' object that contains the log information 
		chosen by the caller through a GUI dialog.
		
		"""
		self.numberOneSon = tk.Toplevel( parent )				# tk window under the parent from calling routine
		self.aBunchOfLogs = aBunchOfLogs						# BunchOfLogs object from calling function
		self.lasPath = str(mm.findFileGui()).split("'")[1]		# GUI to find las file and return it's path
		
		if self.lasPath == "":									# If no file selected, exit the process
			mm.myStdErrorMessage( "Read LAS error", "No file was select so will exit process.")
			sys.exit()
			
		self.lasFileName = self.ladPath.split("/")				# Retreive the LAS filename
		self.lasFileName = self.lasFileName[len(self.lasFileName)-1]..replace(".LAS", "").replace(".las", "")
		
		self.logsAvailableInFile = mm.CurvesInfoLAS()			# Retreive available log names in LAS file
		self.logsAvailableInFile = mm.ReadLogTypesFromLasFile( self.lasPath )
		
		self.logChoiceNumber = IntVar( master = self.numberOneSon )
		self.logChoiceNumber.set( 0 )							# Initialize radio button choice to first log curve
		
		self.theLog = mm.Log()									# Container for a particular log curve
		
		self.aBunchOfLogs.setnumLogs( self.logsAvailableInFile.numCurves )
		self.aBunchOfLogs.setlogType( self.logsAvailableInFile.logMnemonic )
		
		def readAllLogsFromLasFile( self , filePath ):
			
			"""
		
			Attempts to read all detected logs in a LAS file pointed to
			by 'filePath' into the 'aBunchOfLogs' object
		
			"""
			self.lasFile = open( self.lasPath, "r")
		
			self.flag = 0
			
			for self.i in range( self.logsAvailableInFile.numCurves ):
				self.aBunchOfLogs.theLog[ self.i ].ValueUnit = self.logsAvailableInFile.logUnit[ self.i ]
		
				for self.lasRow in self.lasFile:
			
			if self.lasRow[0] == '~A':
				self.flag = 1
				
			if self.flag == 1:
				for i in range( self.logsAvailableInFile.numCurves ):
					self.aBuchOfLogs.TheLog[i].Value.append( float( self.lasRow.split()[ self.i ] ) )
				