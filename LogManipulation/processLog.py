import Tkinter as tk
import MyModule as mm

class processLog:
	"""
	
	A container for general log processing algorithms. A new log
	is generally created and appended to the input 'Bunch of Logs'
	
	"""
	def __init__( self, root, aBunchOfLogs = mm.BunchOfLogs ):
		"""
		Main setup for log processing.
		"""
		self.root = root
		self.localBunchOfLogs = aBunchOfLogs
		
		print "+------------------------------------------+"
		print "|                                          |"
		print "|     Start of log processing function     |"
		print "|                                          |"
		print "+------------------------------------------+"
		
		self.branch1 = tk.Toplevel( self.root )
		self.branch1.wm_title( "Basic Log Processing Function" )
		
		self.logChoice = mm.getALogFromBunchOfLogs()
		
		self.rowCounter = 0
		
		self.basicDescriptionLable = tk.Label( self.branch1,
								text = "Basic Log Processing Functions",
								font=("Comic Sans MS",16), 
								fg="Blue" )
		self.basicDescriptionLable.grid( row = self.rowCounter, 
								columnspan = 2,
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 30 )
								
		self.rowCounter += 1
		
		self.unitConversionLabel = tk.Label( self.branch1, 
								text = "Log unit conversion",
								font=("Comic Sans MS",16), 
								fg="Blue" )
		self.unitConversionLabel.grid( row = self.rowCounter, 
								column = 0,
								sticky = tk.W, 
								padx = 20, 
								pady = 5 )
		
		self.unitConversionButton = tk.Button( self.branch1,
								text = "Convert Unit",
								command = self.unitConvert,
								font=("Comic Sans MS",14), 
								fg="Black" )
		self.unitConversionButton.grid( row = self.rowCounter, 
								column = 1, 
								sticky = tk.E+tk.W,
								padx = 20, 
								pady = 5 )
								
		self.rowCounter += 1
		
		self.sonicVelocityTransformLabel = tk.Label( self.branch1, 
								text = "Transform sonic travel time to velocity",
								font=("Comic Sans MS",16), 
								fg="Blue" )
		self.sonicVelocityTransformLabel.grid( row = self.rowCounter, 
								column = 0, 
								sticky = tk.W,
								padx = 20, 
								pady = 5 )
		
		self.sonicVelocityTransformButton = tk.Button( self.branch1,
								text = "Time to Vel",
								command = self.sonicVelocityTransform,
								font=("Comic Sans MS",14), 
								fg="Black" )
		self.sonicVelocityTransformButton.grid( row = self.rowCounter, 
								column = 1, 
								sticky = tk.E+tk.W,
								padx = 20, 
								pady = 5 )
								
		self.rowCounter += 1
		
		self.castagnaVpVsTransformLabel = tk.Label( self.branch1, 
								text = "Castagna Vp to Vs transform",
								font=("Comic Sans MS",16), 
								fg="Blue" )
		self.castagnaVpVsTransformLabel.grid( row = self.rowCounter, 
								column = 0, 
								sticky = tk.W,
								padx = 20, 
								pady = 5 )
		
		self.castagnaVpVsTransformButton = tk.Button( self.branch1,
								text = "Vp to Vs",
								command = self.castagnaVpVsTransform,
								font=("Comic Sans MS",14), 
								fg="Black" )
		self.castagnaVpVsTransformButton.grid( row = self.rowCounter, 
								column = 1, 
								sticky = tk.E+tk.W,
								padx = 20, 
								pady = 5 )
								
		self.rowCounter += 1
		
		self.gardnerDensityVpTransformLabel = tk.Label( self.branch1, 
								text = "Gardner's Density to Vp transform",
								font=("Comic Sans MS",16), 
								fg="Blue" )
		self.gardnerDensityVpTransformLabel.grid( row = self.rowCounter, 
								column = 0, 
								sticky = tk.W,
								padx = 20, 
								pady = 5 )
		
		self.gardnerDensityVpTransformButton = tk.Button( self.branch1,
								text = "Density to Vp",
								command = self.gardnerDensityVpTransform,
								font=("Comic Sans MS",14), 
								fg="Black" )
		self.gardnerDensityVpTransformButton.grid( row = self.rowCounter, 
								column = 1, 
								sticky = tk.E+tk.W,
								padx = 20, 
								pady = 5 )
		
		self.rowCounter += 1
		
		self.exitButton = tk.Button( self.branch1,
								text = "Exit", 
								font = ("Comic Sans MS",16), 
								fg = "Red", 
								command = self.branch1.destroy )
		self.exitButton.grid( row = self.rowCounter, 
								columnspan = 2, 
								# sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 20 )
		
		self.rowCounter += 1
		
		
# **********************************************
# **********************************************
# ***                                        ***
# ***   Start of unit conversion functions   ***
# ***                                        ***
# **********************************************
# **********************************************
		
	def unitConvert( self ):
		"""
		Unit conversion function
		"""
		
		self.FirstTimeOnly = 0
		
		self.branch2 = tk.Toplevel( self.branch1 )
		self.branch2.wm_title( "Unit Conversion Funation")
		
		if self.FirstTimeOnly == 0:
			self.logChoice = mm.getALogFromBunchOfLogs( self.branch2 )
			self.FirstTimeOnly = 1
		
		self.fromUnit = tk.StringVar( self.branch2 )
		self.fromUnit.set( "None" )
		self.toUnit = tk.StringVar( self.branch2 )
		self.toUnit.set( "None" )
		self.displayLogId = tk.StringVar( self.branch2 )
		self.displayLogId.set( str( self.logChoice.logId ) )
		self.displayLogType = tk.StringVar( self.branch2 )
		self.displayLogType.set( str( self.logChoice.logType ) )
		self.displayLogUnit = tk.StringVar( self.branch2 )
		self.displayLogUnit.set( str( self.logChoice.logUnit ) )
		
		# Set up allowed units for conversion. Not all members of this list
		# can be converted to the others, as can be seen by its members.
		# Current convention is to keep all the characters in lower case.
		self.unitList = [ 	'm', 
							'ft',
							'km',
							'm/s',
							'ft/s',
							'km/s'  ]
							
		self.conversionScalarDictionary = {}
		for u in self.unitList:						# Initialize conversion matrix to null values
			for v in self.unitList:
				self.conversionScalarDictionary[ (u,v) ] = None
				
		# Setting actual possible conversion scalars within conversion matrix
		
		#  For Lengths
		self.conversionScalarDictionary[ ('m','ft')  ] = 3.280839895
		self.conversionScalarDictionary[ ('m','m')   ] = 1.0
		self.conversionScalarDictionary[ ('m','km')  ] = 0.001
		self.conversionScalarDictionary[ ('ft','ft') ] = 1.0
		self.conversionScalarDictionary[ ('ft','m')  ] = 1.0 / self.conversionScalarDictionary[('m','ft')]
		self.conversionScalarDictionary[ ('ft','km') ] = 0.001 / self.conversionScalarDictionary[('m','ft')]
		self.conversionScalarDictionary[ ('km','ft') ] = 1.0 / self.conversionScalarDictionary[('ft','km')]
		self.conversionScalarDictionary[ ('km','m')  ] = 1000.0
		self.conversionScalarDictionary[ ('km','km') ] = 1.0
		
		# For Velocities
		self.conversionScalarDictionary[ ('m/s','ft/s')  ] = 3.280839895
		self.conversionScalarDictionary[ ('m/s','km/s')  ] = 0.001
		self.conversionScalarDictionary[ ('m/s','m/s')   ] = 1.0
		self.conversionScalarDictionary[ ('ft/s','m/s')  ] = 1.0 / self.conversionScalarDictionary[('m/s','ft/s')]
		self.conversionScalarDictionary[ ('ft/s','ft/s') ] = 1.0
		self.conversionScalarDictionary[ ('ft/s','km/s') ] = 0.001 / self.conversionScalarDictionary[('m/s','ft/s')]
		self.conversionScalarDictionary[ ('km/s','ft/s') ] = 1.0 / self.conversionScalarDictionary[('ft/s','km/s')]
		self.conversionScalarDictionary[ ('km/s','m/s')  ] = 1000.0
		self.conversionScalarDictionary[ ('km/s','km/s') ] = 1.0
		
		for i in range( len( self.unitList ) ):
			if str( self.unitList[i] ) == str( self.logChoice.logUnit ).lower() :
				self.fromUnit.set( str( self.unitList[i] ) )
				
		self.rowCounter = 0
		
		self.unitConvertMainLable = tk.Label( self.branch2,
								text = "Unit conversion function",
								font = ("Comic Sans MS",16), 
								fg = "Blue" )
		self.unitConvertMainLable.grid( row = self.rowCounter, 
								columnspan = 3, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 10 )
								
		self.rowCounter += 1
		
		self.chooseALogButton = tk.Button( self.branch2,
								text = "Click to select log",
								command = self.selectingLogForConversion,
								font=("Comic Sans MS",16), 
								fg="Black")
		self.chooseALogButton.grid( row = self.rowCounter,
								columnspan = 3, 
								padx = 20, 
								pady = 10)
		
		self.rowCounter += 1
		
		self.logLabel = tk.Label( self.branch2,
								text = "Chosen Log for Conversion",
								font=("Comic Sans MS",16), 
								fg="Green")
		self.logLabel.grid(row = self.rowCounter, 
								columnspan = 3, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5 )
		
		self.rowCounter += 1
		
		self.chosenLogIdLabel = tk.Label( self.branch2,
								text = "Log ID",
								font=("Comic Sans MS",14), 
								fg="Green")
		self.chosenLogIdLabel.grid(row = self.rowCounter,
		 						column = 0, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5 )
								
		self.chosenLogTypeLabel = tk.Label( self.branch2,
								text = "Log Type",
								font=("Comic Sans MS",14), 
								fg="Green")
		self.chosenLogTypeLabel.grid(row = self.rowCounter,
		 						column = 1, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5 )
								
		self.chosenLogUnitLabel = tk.Label( self.branch2,
								text = "Log Unit",
								font=("Comic Sans MS",14), 
								fg="Green" )
		self.chosenLogUnitLabel.grid( row = self.rowCounter,
								column = 2, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5)
		
		self.rowCounter += 1
		
		self.chosenLogId = tk.Label( self.branch2, 
								textvariable = self.displayLogId,
								font=("Comic Sans MS",14), 
								fg="Black" )
		self.chosenLogId.grid( row = self.rowCounter,
								column = 0,
								sticky = tk.E+tk.W, 
								padx = 20 )
								
		self.chosenLogType = tk.Label( self.branch2,
								textvariable = self.displayLogType,
								font=("Comic Sans MS",14), 
								fg="Black" )
		self.chosenLogType.grid( row = self.rowCounter,
								column = 1,
								sticky = tk.E+tk.W, 
								padx = 20 )
								
		self.chosenLogUnit = tk.Label( self.branch2,
								textvariable = self.displayLogUnit,
								font=("Comic Sans MS",14), 
								fg="Black" )
		self.chosenLogUnit.grid( row = self.rowCounter,
								column = 2,
								sticky = tk.E+tk.W, 
								padx = 20)
		
		self.rowCounter += 1
		
		self.lineLabel = tk.Label( self.branch2,
								text = "___________________________________________",
								font = ("Comic Sans MS",12), 
								fg = "Blue")
		self.lineLabel.grid( row = self.rowCounter,
								columnspan = 3,
								sticky = tk.E+tk.W, 
								padx = 20 )
		
		self.rowCounter += 1
		
		self.convetChoiceLabel = tk.Label( self.branch2,
								text = "Choose units for conversion",
								font=("Comic Sans MS",16), 
								fg="Green")
		self.convetChoiceLabel.grid(row = self.rowCounter, 
								columnspan = 3, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5 )
								
		self.rowCounter += 1
		
		self.fromLabel = tk.Label( self.branch2,
								text = "From unit below",
								font=("Comic Sans MS",16), 
								fg="Green")
		self.fromLabel.grid(row = self.rowCounter, 
								column = 0, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5 )
								
		self.toLabel = tk.Label( self.branch2,
								text = "To unit below",
								font = ("Comic Sans MS",16), 
								fg = "Green")
		self.toLabel.grid(row = self.rowCounter, 
								column = 2, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5 )
								
		self.rowCounter += 1
		
		for self.unit in self.unitList:
			
			self.fromButton = tk.Radiobutton( self.branch2,
								text = self.unit,
								variable = self.fromUnit,
								value = self.unit,
								font=("Comic Sans MS",14), 
								fg="Black" )
			self.fromButton.grid(row = self.rowCounter, 
								column = 0, 
								sticky = tk.E+tk.W, 
								padx = 20 )
			
			self.pointerLabel = tk.Label( self.branch2,
								text = "--------->",
								font=("Comic Sans MS",14), 
								fg="Blue" )
			self.pointerLabel.grid( row = self.rowCounter,
								column = 1, 
								sticky = tk.E+tk.W, 
								padx = 20)
			
			self.toButton = tk.Radiobutton( self.branch2,
								text = self.unit,
								variable = self.toUnit,
								value = self.unit,
								font=("Comic Sans MS",14), 
								fg="Black" )
			self.toButton.grid(row = self.rowCounter, 
								column = 2, 
								sticky = tk.E+tk.W, 
								padx = 20 )
								
			self.rowCounter += 1
		
		self.fromLabel = tk.Label( self.branch2,
								text = "Chosen unit to convert",
								font=("Comic Sans MS",16), 
								fg="Green")
		self.fromLabel.grid(row = self.rowCounter, 
								column = 0, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5 )

		self.toLabel = tk.Label( self.branch2,
								text = "Chosen resultant unit",
								font = ("Comic Sans MS",16), 
								fg = "Green")
		self.toLabel.grid(row = self.rowCounter, 
								column = 2, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 5 )

		self.rowCounter += 1
		
		self.finalUnitChoosen = tk.Label( self.branch2,
								textvariable = self.fromUnit,
								font=("Comic Sans MS",14), 
								fg="Black")
		self.finalUnitChoosen.grid( row = self.rowCounter, 
								column = 0, 
								sticky = tk.E+tk.W, 
								padx = 20 )
								
		self.finalUnitChoosen = tk.Label( self.branch2,
								textvariable = self.toUnit,
								font=("Comic Sans MS",14), 
								fg="Black")
		self.finalUnitChoosen.grid( row = self.rowCounter, 
								column = 2, 
								sticky = tk.E+tk.W, 
								padx = 20 )
		self.rowCounter += 1
		
		self.lineLabel = tk.Label( self.branch2,
								text = "___________________________________________",
								font = ("Comic Sans MS",12), 
								fg = "Blue")
		self.lineLabel.grid( row = self.rowCounter,
								columnspan = 3,
								sticky = tk.E+tk.W, 
								padx = 20 )
		
		self.rowCounter += 1
		
		self.submitForUnitChange = tk.Button( self.branch2,
								text = "Submit for Conversion",
								command = self.applyUnitChange,
								font=("Comic Sans MS",14), 
								fg="Black")
		self.submitForUnitChange.grid( row = self.rowCounter,
								column = 1, 
								sticky = tk.E+tk.W, 
								padx = 20)
								
		self.rowCounter += 1
		
		self.exitButton = tk.Button( self.branch2,
								text = "Exit", 
								font = ("Comic Sans MS",14), 
								fg = "Red", 
								command = self.branch1.destroy )
		self.exitButton.grid( row = self.rowCounter, 
								column = 1, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 10 )
								
		self.branch2.mainloop()
		
	def applyUnitChange( self ):
		"""
		Front end to apply the changes and append the new log
		"""
		self.newLog = self.localBunchOfLogs.theLog[ self.logChoice.logNumber ]
		self.newLog.setValueUnit( self.toUnit )
		
		self.branch3 = tk.Toplevel( self.branch2 )
		self.branch3.wm_title( "Set new log ID if desired" )
		
		self.newLogID = tk.StringVar( self.branch3 )
		self.newLogID.set( str( self.localBunchOfLogs.MyLogLabel[ self.logChoice.logNumber ] ) )
		
		self.rowCounter = 0
		
		self.idChangeLabel = tk.Label( self.branch3,
								text = "Change the log ID below if desired then click Submit",
								font = ("Comic Sans MS",16), 
								fg = "Green")
		self.idChangeLabel.grid( row = self.rowCounter,
								column = 0,
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 10 )
		self.rowCounter += 1
		
		self.idChangeEntry = tk.Entry( self.branch3,
								textvariable = self.newLogID,
								font=("Comic Sans MS",14), 
								fg="Black" )
		self.idChangeEntry.grid( row = self.rowCounter,
								column = 0,
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 10 )
								
		self.rowCounter += 1
			
		self.actuallySubmitUnitChange = tk.Button( self.branch3,
								text = "Submit and Exit",
								command = self.applyUnitChange,
								font=("Comic Sans MS",14), 
								fg="Black")
		self.actuallySubmitUnitChange.grid( row = self.rowCounter,
								column = 0, 
								sticky = tk.E+tk.W, 
								padx = 20)
			
		self.rowCounter += 1
		
		self.exitButton = tk.Button( self.branch3,
								text = "Exit", 
								font = ("Comic Sans MS",14), 
								fg = "Red", 
								command = self.branch3.destroy )
		self.exitButton.grid( row = self.rowCounter, 
								column = 1, 
								sticky = tk.E+tk.W, 
								padx = 20, 
								pady = 10 )
								
		self.branch3.mainloop()
		
	def actuallySubmitUnitChange( self ):
		"""
		Perform the mechanical work of updating the BunchofLogs
		"""
		
		print "\n    Start conversion from (" +str( self.fromUnit.get() ) + ") to (" + str( self.toUnit.get() ) + ")"
		print "    Conversion scalar value : " + str( self.conversionScalarDictionary[ ( self.fromUnit.get(), self.toUnit.get() ) ] )
		
		self.chosenScalar = self.conversionScalarDictionary[ ( self.fromUnit.get(), self.toUnit.get() ) ]
		self.newLog.Value[:] = [ x*self.chosenScalar for x in self.newLog.Value ]
		
		aBunchOfLog.appendALog( self.newLog, self.newLogID.get() )
		
		print " Completed unit conversion of log and added to project logs\n"
		
		self.branch3.destroy()
		
	def selectingLogForConversion( self ):
		"""
		Selection of a log from a bunch of logs for 
		"""
		
		self.logChoice = mm.getALogFromBunchOfLogs( self.branch2, self.localBunchOfLogs )
		
		
		print "\n Update information for chosen log:\n"
		print "     Log ID   : " + str(self.logChoice.logId)
		print "     Log Type : " + str(self.logChoice.logType)
		print "     Log Unit : " + str(self.logChoice.logUnit) + "\n"
		
		
		self.displayLogId.set( self.logChoice.logId )
		self.displayLogType.set( self.logChoice.logType )
		self.displayLogUnit.set( self.logChoice.logUnit )
	
#
# ********************************************
# ********************************************
# ***                                      ***
# ***   End of unit conversion functions   ***
# ***                                      ***
# ********************************************
# ********************************************
#

#
# *****************************************************************
# *****************************************************************
# ***                                                           ***
# ***   Start of transit time to velocity transform functions   ***
# ***                                                           ***
# *****************************************************************
# *****************************************************************
#
	def sonicVelocityTransform( self ):
		"""
		Nothing yet
		"""	
		print "Transit time to velocity transform not yet implimented"
#
# *****************************************************************
# *****************************************************************
# ***                                                           ***
# ***     End of transit time to velocity transform functions   ***
# ***                                                           ***
# *****************************************************************
# *****************************************************************
#


#
# *************************************************
# *************************************************
# ***                                           ***
# ***   Start of Vp to Vs transform functions   ***
# ***                                           ***
# *************************************************
# *************************************************
#
	def castagnaVpVsTransform( self ):
		"""
		Nothin yet
		"""
		print "Catagna's transform not yet implimented"
#
#
# *************************************************
# *************************************************
# ***                                           ***
# ***     End of Vp to Vs transform functions   ***
# ***                                           ***
# *************************************************
# *************************************************
#

#
# ***********************************************************************
# ***********************************************************************
# ***                                                                 ***
# ***   Start density to compressional velocity transform functions   ***
# ***                                                                 ***
# ***********************************************************************
# ***********************************************************************
#
	def gardnerDensityVpTransform( self ):
		"""
		Nothing yet
		"""
		print "Gardner's transform not yet implimented"
#
# ***********************************************************************
# ***********************************************************************
# ***                                                                 ***
# ***     End density to compressional velocity transform functions   ***
# ***                                                                 ***
# ***********************************************************************
# ***********************************************************************
#