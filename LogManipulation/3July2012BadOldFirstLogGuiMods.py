#!/usr/bin/python

from Tkinter import *
from MyModule import *
import tkFileDialog
from SimPy.SimPlot import *

class BunchOfLogs:
	def __init__( self, MyLogLabel=[], LogType=[] , TheLog=[Log()], NumLogs=0 ):
		"""
		This class is ment to hold a bunch of logs for further processing.
		BuchOfLogs takes three arguments (lists), none of which have to be specified
		"""
		self.MyLogLabel = MyLogLabel					# An arbitrary ID for a log like its well name
		self.LogType = LogType							# The log type associated with corresponding label
		self.TheLog = TheLog							# The values and other information of corresponding log
		self.NumLogs = NumLogs							# Number of logs in this bunch of Logs
		
		def setMyLogLabel( self, MyLogLabel ):
			"""
			Takes one list, the list of user defined ID's for the associated log.
			"""
			self.MyLogLabel = MyLogLabel
			
		def setLogType( self, LogType ):
			"""
			Takes one list, the list of log types corresponding to defined ID's list.
			"""
			self.LogType = LogType
			
		def setTheLog( self, TheLog ):
			"""
			Takes one list, the list of values of the logs and additonal information 
			that can be found in class "Log" in "MyLogsMod" module.
			"""
			self.TheLog = TheLog
			
		def setNumLogs( self, NumLogs ):
			"""
			Takes one argument, the number of logs in this bunch of logs.
			"""
			self.NumLogs = NumLogs
			
def findFileGui():
	"""
	Opens up a dialog box for finding files and returns the result
	"""
	branch1 = Tk()
	branch1.withdraw()
	
	pathFile = tkFileDialog.askopenfile(parent=branch1,title='Choose a file')
	
	if pathFile != None:
		return pathFile
	else:
		myStdErrorMessage("File read error", "No file was detected")
		
	branch1.destroy
	
def append2BunchOfLogs( logsAreUs=BunchOfLogs(), AddLog=Log(), AddLogID="" ):
	"""
	Append a new log to a "BunchOfLogs" object.
	"""
	print "+-----------------------------------------------------+"
	print "|                                                     |"
	print "|           Entering Append a 'Log' Object            |"
	print "|           to a 'BunchOfLogs' Object Module          |"
	print "|                                                     |"
	print "+-----------------------------------------------------+"

	if AddLogID == "":
		AddLogID = AddLog.ValueID
		
	logsAreUs.NumLogs = logsAreUs.NumLogs + 1
	logsAreUs.MyLogLabel.append(AddLogID)
	logsAreUs.LogType.append(AddLog.ValueID)
	logsAreUs.TheLog.insert(int(logsAreUs.NumLogs - 1), AddLog)
	
	myStdErrorMessage( "Append Log Report", 
						"Number of logs: " + str(logsAreUs.NumLogs) + 
						"\n Appended log ID: " + str(logsAreUs.MyLogLabel) + 
						"\n Appended log Type: " + str(logsAreUs.LogType))
	
	print ""
	print "New number of Logs    : " + str( logsAreUs.NumLogs )
	print "New Log Identifier    : " + str( logsAreUs.MyLogLabel[ logsAreUs.NumLogs - 1 ] )
	print "Type of log appended  : " + str( logsAreUs.LogType[ logsAreUs.NumLogs - 1 ] )
	print "Units for appended log: " + str( AddLog.ValueUnit )
	print "Units after append    : " + str( logsAreUs.TheLog[ logsAreUs.NumLogs -1].ValueUnit )

def findReadLogLasFile(logsAreUs=BunchOfLogs()):
	import sys
	global CurveChoiceNumber, chosenCurveNumber, chosenLogId, curvesMnemonic, chosenCurveActualValues, chosenCurveValueUnit
	"""
	From user chosen LAS file, read all the available curves in the file,
	then allow user to choose one to be output into a user input Log object 
	(The class Log can be found in the "MyLogsMod" module)
	"""
	print "+-------------------------------------------------------------+"
	print "|                                                             |"
	print "|           Entering Find and Read LAS File Module            |"
	print "|                                                             |"
	print "+-------------------------------------------------------------+"
	
	LasPath = str(findFileGui()).split("'")[1]				# Get the path to the LAS file

	if LasPath == "":										# If no file selected, exit the process
		myStdErrorMessage( "Read LAS error", "No file was select so will exit process.")
		sys.exit()
		
															# Find the name of the LAS file from LasPath
	splitLasPath = LasPath.split("/")						# the .las or .LAS portion will be removed
	LasName = splitLasPath[len(splitLasPath) - 1].replace(".LAS", "").replace(".las", "")		
	
	LogsAvailable = CurvesInfoLAS()
	LogsAvailable = ReadLogTypesFromLasFile(LasPath)		# Import the curve names available in choosen file
	
	branch1 = Tk()
	branch1.withdraw
	branch1.wm_title("Curve selection menu")

	CurveChoiceNumber = IntVar(master = branch1)
	CurveChoiceNumber.set(0)								# Initialize radio button choice to first curve
	
	theLog = Log											# Initialeze a "Log" object to hold chosen curve values
	

	
	def setupVariablesPrint():
		global chosenCurveNumber, chosenLogId, curvesMnemonic, chosenCurveActualValues, chosenCurveValueUnit
		chosenCurveNumber = CurveChoiceNumber.get()
		chosenLogId = LogIdEntry.get()
		curvesMnemonic = LogsAvailable.LogMnemonic[chosenCurveNumber]
		chosenCurveActualValues = Read1LogFromLasFileNoWrap( LasPath, curvesMnemonic )
		chosenCurveValueUnit = LogsAvailable.LogUnit[chosenCurveNumber]
		chosenCurveActualValues.ValueUnit = chosenCurveValueUnit
		print "\n\nChosen curve from LAS file below.\n"
		print "Curve number: " + str(chosenCurveNumber + 1)
		print "Curve name  : " + str(curvesMnemonic)
		print "Cure unit   : " + str(chosenCurveActualValues.ValueUnit)
		print "Log ID	   : " + str(chosenLogId)

	
	# Defining the widgets in the window and its grid geometry
	
	RowCounter = 0											# Initialize row for first row widget

	instructLabel = Label( branch1, text="Enter a log Id and select an available curve", 
					font=("Comic Sans MS",18), fg="Blue")
	instructLabel.grid( row=RowCounter, columnspan=3, sticky=E+W, padx=50, pady=20 )
	RowCounter += 1											# Incrementing couter for the next row of widgets
	
	LogIdLablel = Label( branch1, text="Enter a Log ID: ", font=("Comic Sans MS",14), fg="Blue")
	LogIdLablel.grid( row=RowCounter, column=0, sticky=E+S, pady=10 )
	
	LogIdEntry = Entry( branch1)							# Entry box for a user determined ID name 
	LogIdEntry.delete(0, END)
	LogIdEntry.insert(0, LasName)							# Set default ID to LAS filename
	LogIdEntry.grid ( row=RowCounter, column=1, columnspan=2, sticky=E+W+N+S, pady=10, padx=20 )
	RowCounter += 1

	
	LogIdDefaultLabel1 = Label( branch1, text="		   ")
	LogIdDefaultLabel1.grid( row=RowCounter, column=0)
	
	LogIdDefaultLabel2 = Label( branch1, text="Default will be LAS file name")
	LogIdDefaultLabel2.grid( row=RowCounter, column=1, columnspan=2, sticky=W, padx=100)
	RowCounter += 1
	
	
	LogNameLabel = Label(branch1, text="Pick a curve below", font=("Comic Sans MS",16), fg="Blue")
	LogNameLabel.grid( row=RowCounter, column=0, sticky=W, padx=20, pady=10)
	RowCounter += 1
	
	NumberOfCurvesInLasFile = int(LogsAvailable.NumCurves)
	
	for i in range(NumberOfCurvesInLasFile):
		Radiobutton( branch1, 
					text = str(i+1)+": " + str(LogsAvailable.LogMnemonic[i]),
					padx = 20,
					variable = CurveChoiceNumber, 
					command = setupVariablesPrint,
					value = int(i) ).grid( row=RowCounter, column=0, sticky=W, padx=20)
		LabelRadB1 = Label( branch1, text="( " + str(LogsAvailable.LogDescription[i]) + " )" )
		LabelRadB1.grid( row=RowCounter, column=1, sticky=W, padx=20)
		LabelRadB2 = Label( branch1, text="( " + str(LogsAvailable.LogUnit[i]) + " )" )
		LabelRadB2.grid( row=RowCounter, column=2, sticky=W, padx=20)
		RowCounter += 1

	chosenCurveNumber = CurveChoiceNumber.get()
	chosenLogId = LogIdEntry.get()
	curvesMnemonic = LogsAvailable.LogMnemonic[chosenCurveNumber]
	chosenCurveActualValues = Read1LogFromLasFileNoWrap( LasPath, curvesMnemonic )
	chosenCurveActualValues.ValueUnit = LogsAvailable.LogUnit[chosenCurveNumber]
	
	LineSpaceLabel1 = Label(branch1, text="			  ")
	LineSpaceLabel1.grid( row=RowCounter, columnspan=2, sticky=E+W, padx=20, pady=20 )
	
	SubmitButton = Button( branch1, text="	  Submit		", 
					command=lambda: append2BunchOfLogs( logsAreUs, chosenCurveActualValues, chosenLogId ) )
	SubmitButton.grid( row=RowCounter, column=1, sticky=E+W, padx=20 )
	RowCounter += 1
	
	ExitButton = Button( branch1, text="  Exit Process	  ", command=branch1.destroy )
	ExitButton.grid(row=RowCounter, column=1, sticky=E+W, padx=20 )
	RowCounter += 1
	
	LineSpaceLabel2 = Label(branch1, text="			  ")
	LineSpaceLabel2.grid( row=RowCounter, columnspan=2, sticky=E+W, padx=20, pady=20 )

	
	branch1.mainloop()
	
def DisplayContentsOfBunchOfLogs(theBunchOfLogs = BunchOfLogs()):

	print "+-------------------------------------------------------------+"
	print "|                                                             |"
	print "|           Entering program to list all available            |"
	print "|           logs usable within the current session.           |"
	print "|           The logs are assumed to be stored in a            |"
	print "|           'BunchOfLogs' object which stores a bunch         |"
	print "|           of 'Log' objects.                                 |"
	print "|                                                             |"
	print "+-------------------------------------------------------------+"

	
	branch = Tk()
	branch.withdraw()
	
	branch1 = Toplevel(branch)
	rowCounter = 0
	
	if theBunchOfLogs.NumLogs == 0 :
		
		branch1.wm_title("Display Available Logs")
		
		noLogsLabel = Label(branch1, text="It appears there are no logs in current session.", font=("Comic Sans MS",16), fg="Blue")
		noLogsLabel.grid(row=rowCounter, columnspan=3, sticky=E+W, padx=20, pady=20)
		rowCounter += 1
		
		ExitButton = Button( branch1, text="  Exit Process	  ", command=branch1.destroy, font=("Comic Sans MS",14), bg="Pink" )
		ExitButton.grid(row=rowCounter, column=1, padx=20  )
		rowCounter += 1 
			
	else:
		
		branch1.wm_title("Display Available Logs")
		
		descriptionLabel = Label(branch1, text=str(theBunchOfLogs.NumLogs)+" Log curves available for current session.", font=("Comic Sans MS",16), fg="Blue")
		
		
		descriptionLabel.grid(row=rowCounter, columnspan=3, sticky=E+W, padx=20, pady=5)
		rowCounter += 1
		
		logID_label = Label( branch1, text="Log ID", font=("Comic Sans MS",14), fg="Green" )
		logType_label = Label( branch1, text="Log Type", font=("Comic Sans MS",14), fg="Green" )
		logUnit_label = Label( branch1, text="Log Unit", font=("Comic Sans MS",14), fg="Green" )
		
		logID_label.grid( row=rowCounter, column=0, sticky=E+W, padx=20 )
		logType_label.grid( row=rowCounter, column=1, sticky=E+W, padx=20 )
		logUnit_label.grid( row=rowCounter, column=2, sticky=E+W, padx=20 )
		
		rowCounter += 1
		
		for i in range(int(theBunchOfLogs.NumLogs)):
			
			logID_value = Label( branch1, text=str(theBunchOfLogs.MyLogLabel[i]), font=("Comic Sans MS",12) )
			logType_value = Label( branch1, text=str(theBunchOfLogs.LogType[i]), font=("Comic Sans MS",12) )
			logUnit_value = Label( branch1, text=str(theBunchOfLogs.TheLog[i].ValueUnit), font=("Comic Sans MS",12) )

			logID_value.grid( row=rowCounter, column=0, sticky=E+W, padx=20 )
			logType_value.grid( row=rowCounter, column=1, sticky=E+W, padx=20 )
			logUnit_value.grid( row=rowCounter, column=2, sticky=E+W, padx=20 ) 
			
			rowCounter += 1
			
		ExitButton = Button( branch1, text="  Exit Process	  ", command=branch1.destroy, font=("Comic Sans MS",14), fg="Red" )
		ExitButton.grid(row=rowCounter, column=1, padx=20, pady=20  )
		rowCounter += 1		
	
	branch.mainloop()
	
chosenNumber = -999
def ChooseALog( aBunchOfLogs = BunchOfLogs() ):
	""" 
	
	This function sets up a gui to allow the choice of a particular
	log within a 'BunchOfLogs' object. The idea is to return the
	chosen log number within the 'BunchOfLogs' collection of logs
	
	"""
	print "+--------------------------------------------------------+"
	print "|                                                        |"
	print "|   Entering log selection from a bunch of logs object   |"
	print "|                                                        |"
	print "+--------------------------------------------------------+"
	
	
	
	
	branch1 = Tk()
	#branch1.withdraw()
	
	global logChoiceNumber, chosenNumber	
	logChoiceNumber = IntVar()
	logChoiceNumber.set(-999)
	chosenNumber = -999
		
	rowCounter = 0
	
	
	def setChosenNumber():
		
		if logChoiceNumber.get() == -999:
			print "No choice has been detected yet."
			print "Press a radio button to make a choice."
			
		else:
			print "You have chosen log number: " + str( logChoiceNumber.get() )
			
			
		chosenNumber = logChoiceNumber.get()
		# return chosenNumber
		print "I have at least exited the if block with chosenNumber: "+str(chosenNumber)
	
	
	if aBunchOfLogs.NumLogs == 0 :
		
		branch1.wm_title("Choose An Available Log")
		
		noLogsLabel = Label(branch1, text="It appears, there are no logs in current session.", font=("Comic Sans MS",16), fg="Blue")
		noLogsLabel.grid(row=rowCounter, columnspan=3, sticky=E+W, padx=20, pady=20)
		rowCounter += 1
		
		ExitButton = Button( branch1, text="  Exit Process	  ", command=chooseLogExit, font=("Comic Sans MS",14), bg="Pink" )
		ExitButton.grid(row=rowCounter, column=1, padx=20  )
		rowCounter += 1 
			
	else:
		
		branch1.wm_title("Choose An Available Log")
		
		descriptionLabel = Label(branch1, text="Choose one of "+str(aBunchOfLogs.NumLogs)+" Log curves available for current session.", font=("Comic Sans MS",16), fg="Blue")
		descriptionLabel.grid(row=rowCounter, columnspan = 3, sticky=E+W, padx=20, pady=5)
		rowCounter += 1
		
		logID_label = Label( branch1, text="Log ID", font=("Comic Sans MS",14), fg="Green" )
		logType_label = Label( branch1, text="Log Type", font=("Comic Sans MS",14), fg="Green" )
		logUnit_label = Label( branch1, text="Log Unit", font=("Comic Sans MS",14), fg="Green" )
		
		logID_label.grid( row=rowCounter, column=0, sticky=E+W, padx=20 )
		logType_label.grid( row=rowCounter, column=1, sticky=E+W, padx=20 )
		logUnit_label.grid( row=rowCounter, column=2, sticky=E+W, padx=20 )
		
		rowCounter += 1
		
		for i in range(int(aBunchOfLogs.NumLogs)):
			
			global logChoiceNumber
			
			logID_value = Label( branch1, text=str( aBunchOfLogs.MyLogLabel[i] ), font=("Comic Sans MS",12) )
			logTypeValueButton = Radiobutton( branch1, text=str( aBunchOfLogs.LogType[i] ), font=("Comic Sans MS",12),
			 						variable = logChoiceNumber, value=int(i), command=setChosenNumber )
			logUnit_value = Label( branch1, text=str(aBunchOfLogs.TheLog[i].ValueUnit), font=("Comic Sans MS",12) )
	
			logID_value.grid( row=rowCounter, column=0, sticky=E+W, padx=20 )
			logTypeValueButton.grid( row=rowCounter, column=1, sticky=E+W, padx=20 )
			logUnit_value.grid( row=rowCounter, column=2, sticky=E+W, padx=20 ) 
			
			rowCounter += 1
		
		ExitButton = Button( branch1, text="  Exit Process	  ", command=branch1.destroy, font=("Comic Sans MS",14), fg="Red" )
		ExitButton.grid(row=rowCounter, column=1, padx=20, pady=20  )
		rowCounter += 1
	
	print "Does it get here ???"+str(logChoiceNumber.get())
	chosenNumber = logChoiceNumber.get()
	return chosenNumber
	
	branch1.mainloop()

def PlotALogWithinBunchOfLogs( availableBunchOfLogs = BunchOfLogs() ):
	"""
	
	This function will allow the user to select two logs stored in a 'BunchOfLogs' object
	and produce a cross-plot of the two logs. An example would be, letting x values come from
	sonic travel time values and y values from the measured depth log, this would then produce a
	familiar log section with depth. Alternately we could plot any log againts any other log.
	
	"""
	print "+--------------------------------------------------+"
	print "|                                                  |"
	print "|    Entering function to created an (x,y) plot    |"
	print "|    from two user selected log curves available   |"
	print "|    in the current session.                       |"
	print "|                                                  |"
	print "+--------------------------------------------------+"
	print ""
	
	branch1 = Tk()
	branch1.withdraw
	branch1.wm_title("Log X-Y Plot")

	global logX, logY, xNumber, yNumber, xLabelTypeText, xLabelIdText, yLabelTypeText, yLabelIdText
	
	xLabelTypeText = StringVar( master = branch1 )
	xLabelTypeText.set("None")
	
	xLabelIdText = StringVar( master = branch1 )
	xLabelIdText.set("None")
	
	yLabelTypeText = StringVar( master = branch1 )
	yLabelTypeText.set("None")
	
	yLabelIdText = StringVar( master = branch1 )
	yLabelIdText.set("None")
	
	xNumber = -999									# Setting value to indicate no x-value log has been chosen yet
	yNumber = -999									# Setting value to indicate no y-value log has been chosen yet
	

	
	def GetX_Log(branch):
		
		global logX, xNumber, xLabelIdText, xLabelTypeText
		
		print "Getting log for X-values"
		#xNumber = ChooseALog(availableBunchOfLogs)
		ChooseALog(availableBunchOfLogs)
		xNumber = chosenNumber
		print "Hello Dolly X: "+str(xNumber)
		
		if xNumber != -999:
			print "X log number: "+str(xNumber)+"\n"
			logX = availableBunchOfLogs.TheLog[ int(xNumber) ].Value
			print "logX has been set to values of log number: "+str(xNumber)+"\n"
			xLabelIdText.set( str( availableBunchOfLogs.MyLogLabel[ int(xNumber) ] ) )
			xLabelTypeText.set( str( availableBunchOfLogs.LogType[ int(xNumber) ] ) )
		
	def GetY_Log(branch):
		
		global logY, yNumber, yLabelIdText, yLabelTypeText
		
		print "Getting log for Y-values"
		#yNumber = ChooseALog(availableBunchOfLogs)
		ChooseALog(availableBunchOfLogs)
		yNumber = chosenNumber
		
		if yNumber != -999:
			print "Y log number: "+str(yNumber)+"\n"
			logY = availableBunchOfLogs.TheLog[ int(yNumber) ].Value
			print "Hello Dolly Y: "+str(xNumber)
			print "logY has been set to values of log number: "+str(xNumber)+"\n"
			yLabelIdText.set( str( availableBunchOfLogs.MyLogLabel[ int(yNumber) ] ) )
			yLabelTypeText.set( str( availableBunchOfLogs.LogType[ int(yNumber) ] ) )
		
	def MyXYPlot(X, Y):
		
		Z = zip(X,Y)
		plt=SimPlot()
		plt.plotLine(Z)
		plt.mainloop()
		
	rowCounter = 0
	
	mainLabel = Label( branch1, text="Choose logs for (x,y) values, then press submit to plot", font=("Comic Sans MS",16), fg="Blue")
	mainLabel.grid( row = rowCounter, columnspan=3, sticky = E+W, padx = 20, pady = 20 )
	rowCounter += 1
	
	idLabel = Label( branch1, text="Log ID", font=("Comic Sans MS",14), fg="Green" )
	idLabel.grid( row = rowCounter, column=1, sticky = E+W, padx = 20, pady = 5 )
	typeLabel = Label( branch1, text="Log type", font=("Comic Sans MS",14), fg="Green" )
	typeLabel.grid( row = rowCounter, column=2, sticky = E+W, padx = 20, pady = 5 )
	rowCounter += 1
	
	xButton = Button( branch1, text = "Choose Log for X value", command =lambda: GetX_Log(branch1), font=("Comic Sans MS",14) )
	xButton.grid( row = rowCounter, column=0, sticky = E+W, padx = 20, pady = 5 )
	xLabel_Id = Label( branch1, text="default", textvariable = xLabelIdText, font=("Comic Sans MS",14))
	xLabel_Id.grid( row = rowCounter, column=1, sticky = E+W, padx = 20, pady = 5 )
	xLabel_Type = Label( branch1, text="default", textvariable = xLabelTypeText, font=("Comic Sans MS",14))
	xLabel_Type.grid( row = rowCounter, column=2, sticky = E+W, padx = 20, pady = 5 )
	rowCounter += 1
	
	yButton = Button( branch1, text = "Choose Log for Y value", command =lambda: GetY_Log(branch1), font=("Comic Sans MS",14) )
	yButton.grid( row = rowCounter, column=0, sticky = E+W, padx = 20, pady = 5 )
	yLabel_Id = Label( branch1, text = "default", textvariable = yLabelIdText, font=("Comic Sans MS",14))
	yLabel_Id.grid( row = rowCounter, column=1, sticky = E+W, padx = 20, pady = 5 )
	yLabel_Type = Label( branch1, text = "default", textvariable = yLabelTypeText, font=("Comic Sans MS",14))
	yLabel_Type.grid( row = rowCounter, column=2, sticky = E+W, padx = 20, pady = 5 )
	rowCounter += 1	
	
	if xNumber == -999 and yNumber == -999:
		print "No logs chosen to be plotted yet."
	elif xNumber == -999 or yNumber == -999:
		print "Only one log chosen, both need to be chosen for plotting."
	else:
		print "Chosen x-log number : " + str( xNumber )
		print "Chosen x-log ID     : " + str( availableBunchOfLogs.MyLogLabel[ int(yNumber) ] )
		print "Chosen x-log type   : " + str( availableBunchOfLogs.LogType[ int(xNumber) ] )
		print ""
		print "Chosen y-log number : " + str( yNumber )
		print "Chosen y-log ID     : " + str( availableBunchOfLogs.MyLogLabel[ int(yNumber) ] )
		print "Chosen y-log type   : " + str( availableBunchOfLogs.LogType[ int(yNumber) ] )
	
	submitButton = Button( branch1, text="Submit (x,y) for Plotting", font=("Comic Sans MS",14), command = lambda: MyXYPlot(logX,logY) )
	submitButton.grid( row = rowCounter, column=1, sticky = E+W, padx = 20, pady = 20 )
	rowCounter += 1
	
	ExitButton = Button( branch1, text="  Exit Process	  ", command=branch1.destroy, font=("Comic Sans MS",14), fg="Red" )
	ExitButton.grid(row=rowCounter, column=1, padx=20  )
	rowCounter += 1
	
	branch1.mainloop()

	
