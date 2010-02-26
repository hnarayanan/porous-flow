-- Part of the build-system for managing my doctoral dissertation.
-- (c) Harish Narayanan, 2007 -- 2010

on run argv

	tell application "Finder" to get folder of (path to me) as Unicode text
	set workingPath to POSIX path of result
	set inFile to workingPath & item 1 of argv
	set myPDFSettingsPath to item 2 of argv
	
	tell application "Acrobat Distiller"
		Distill sourcePath inFile destinationPath workingPath adobePDFSettingsPath myPDFSettingsPath
		run
		quit
	end tell
	
end run
