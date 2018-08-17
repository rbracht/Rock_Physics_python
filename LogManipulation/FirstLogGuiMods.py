#!/usr/bin/python

from Tkinter import *
from MyModule import *
import tkFileDialog
# from SimPy.SimPlot import *

# def findFileGui():
# 	"""
# 	Opens up a dialog box for finding files and returns the result
# 	"""
# 	branch1 = Tk()
# 	branch1.withdraw()
# 	
# 	pathFile = tkFileDialog.askopenfile(parent=branch1,title='Choose a file')
# 	
# 	if pathFile != None:
# 		return pathFile
# 	else:
# 		myStdErrorMessage("File read error", "No file was detected")
# 		
# 	branch1.destroy
	
def append2BunchOfLogs( logsAreUs=BunchOfLogs(), AddLog=Log(), AddLogID="" ):
	"""
	Append a new log to a "BunchOfLogs" object.
	"""
	print "+-----------------------------------------------------+"
	print "|                                                     |"
	print "|           Entering Append a 'Log' Object            |"
	print "|           to a 'BunchOfLogs' Object function        |"
	print "|                                                     |"
	print "+-----------------------------------------------------+"

	if AddLogID == "":
		AddLogID = AddLog.ValueID
		
	logsAreUs.numLogs = logsAreUs.numLogs + 1
	logsAreUs.myLogLabel.append(AddLogID)
	logsAreUs.logType.append(AddLog.ValueID)
	logsAreUs.theLog.insert(int(logsAreUs.numLogs - 1), AddLog)
	
	myStdErrorMessage( "Append Log Report", 
						"Number of logs: " + str(logsAreUs.numLogs) + 
						"\n Appended log ID: " + str(logsAreUs.myLogLabel) + 
						"\n Appended log Type: " + str(logsAreUs.logType))
	
	print ""
	print "New number of Logs    : " + str( logsAreUs.numLogs )
	print "New Log Identifier    : " + str( logsAreUs.myLogLabel[ logsAreUs.numLogs - 1 ] )
	print "Type of log appended  : " + str( logsAreUs.logType[ logsAreUs.numLogs - 1 ] )
	print "Units for appended log: " + str( AddLog.ValueUnit )
	print "Units after append    : " + str( logsAreUs.theLog[ logsAreUs.numLogs -1].ValueUnit )

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
		curvesMnemonic = LogsAvailable.logMnemonic[chosenCurveNumber]
		chosenCurveActualValues = Read1LogFromLasFileNoWrap( LasPath, curvesMnemonic )
		chosenCurveValueUnit = LogsAvailable.logUnit[chosenCurveNumber]
		chosenCurveActualValues.ValueUnit = chosenCurveValueUnit
		print "\n\nChosen curve from LAS file below.\n"
		print "Curve number: " + str(chosenCurveNumber + 1)
		print "Curve name  : " + str(curvesMnemonic)
		print "Curve unit   : " + str(chosenCurveActualValues.ValueUnit)
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
	
	NumberOfCurvesInLasFile = int(LogsAvailable.numCurves)
	
	for i in range(NumberOfCurvesInLasFile):
		Radiobutton( branch1, 
					text = str(i+1)+": " + str(LogsAvailable.logMnemonic[i]),
					padx = 20,
					variable = CurveChoiceNumber, 
					command = setupVariablesPrint,
					value = int(i) ).grid( row=RowCounter, column=0, sticky=W, padx=20)
		LabelRadB1 = Label( branch1, text="( " + str(LogsAvailable.logDescription[i]) + " )" )
		LabelRadB1.grid( row=RowCounter, column=1, sticky=W, padx=20)
		LabelRadB2 = Label( branch1, text="( " + str(LogsAvailable.logUnit[i]) + " )" )
		LabelRadB2.grid( row=RowCounter, column=2, sticky=W, padx=20)
		RowCounter += 1
	
	Radiobutton( branch1, 
				text = "Select all logs above",
				padx = 20,
				variable = CurveChoiceNumber, 
				command = setupVariablesPrint,
				value = -999 ).grid( row=RowCounter, column=0, sticky=W, padx=20)	
				
	chosenCurveNumber = CurveChoiceNumber.get()
	chosenLogId = LogIdEntry.get()
	curvesMnemonic = LogsAvailable.logMnemonic[chosenCurveNumber]
	chosenCurveActualValues = Read1LogFromLasFileNoWrap( LasPath, curvesMnemonic )
	chosenCurveActualValues.ValueUnit = LogsAvailable.logUnit[chosenCurveNumber]
	
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
	
	if theBunchOfLogs.numLogs == 0 :
		
		branch1.wm_title("Display Available Logs")
		
		noLogsLabel = Label(branch1, text="It appears there are no logs in current session.", font=("Comic Sans MS",16), fg="Blue")
		noLogsLabel.grid(row=rowCounter, columnspan=3, sticky=E+W, padx=20, pady=20)
		rowCounter += 1
		
		ExitButton = Button( branch1, text="  Exit Process	  ", command=branch1.destroy, font=("Comic Sans MS",14), bg="Pink" )
		ExitButton.grid(row=rowCounter, column=1, padx=20  )
		rowCounter += 1 
			
	else:
		
		branch1.wm_title("Display Available Logs")
		
		descriptionLabel = Label(branch1, text=str(theBunchOfLogs.numLogs)+" Log curves available for current session.", font=("Comic Sans MS",16), fg="Blue")
		
		
		descriptionLabel.grid(row=rowCounter, columnspan=3, sticky=E+W, padx=20, pady=5)
		rowCounter += 1
		
		logID_label = Label( branch1, text="Log ID", font=("Comic Sans MS",14), fg="Green" )
		logType_label = Label( branch1, text="Log Type", font=("Comic Sans MS",14), fg="Green" )
		logUnit_label = Label( branch1, text="Log Unit", font=("Comic Sans MS",14), fg="Green" )
		
		logID_label.grid( row=rowCounter, column=0, sticky=E+W, padx=20 )
		logType_label.grid( row=rowCounter, column=1, sticky=E+W, padx=20 )
		logUnit_label.grid( row=rowCounter, column=2, sticky=E+W, padx=20 )
		
		rowCounter += 1
		
		for i in range(int(theBunchOfLogs.numLogs)):
			
			logID_value = Label( branch1, text=str(theBunchOfLogs.myLogLabel[i]), font=("Comic Sans MS",12) )
			logType_value = Label( branch1, text=str(theBunchOfLogs.logType[i]), font=("Comic Sans MS",12) )
			logUnit_value = Label( branch1, text=str(theBunchOfLogs.theLog[i].ValueUnit), font=("Comic Sans MS",12) )

			logID_value.grid( row=rowCounter, column=0, sticky=E+W, padx=20 )
			logType_value.grid( row=rowCounter, column=1, sticky=E+W, padx=20 )
			logUnit_value.grid( row=rowCounter, column=2, sticky=E+W, padx=20 ) 
			
			rowCounter += 1
			
		ExitButton = Button( branch1, text="  Exit Process	  ", command=branch1.destroy, font=("Comic Sans MS",14), fg="Red" )
		ExitButton.grid(row=rowCounter, column=1, padx=20, pady=20  )
		rowCounter += 1		
	
	branch.mainloop()
	


