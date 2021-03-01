========================================================================
    CONSOLE APPLICATION : 204101059_yes_no Project Overview
========================================================================


Vikash Raghuwanshi
204101059

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Idea:
	First calculate STE and ZCR of speech waveform for each frame of 320 samples.
	Then observe the patterns between STE and ZCR of silence, vowels and fricative_sh.
	We observed that, STE of silence is very much low as compared to vowels and fricative_sh.
	And the ZCR of fricative_sh is very much higher as compared to vowels.

	Now, "NO" ends with a vowel "O" and "YES" ends with fricative_sh, which is used for classifying the word.


	Hence, Firstly to make the boundary between silence and a word (either Yes or NO) we use STE. By the observation, 
	some threshold is set for making the boundary.

	Now after getting the boundary, we will store the zcr values of words in the vector.

	Since, we know that the zcr has a very good role at the end of both the words "YES" and "NO" i.e. "sh" and "o" respectively,
	therefore we will take the last 40% data(zcr values) of the word and calculate a zcrCount which counts how many zcr have higher
	values and also takes sum of those values, if we have some higher values at the end and their sum is larger than some observed value 
	then it will be a "YES"

	if it has not too much higher values and sum less than some observed value, then the word is "NO"
	else Ambiguity There

	Some other constraints also applied like if the number of data frames are too less then discard 
	that word, since those may be some other sounds.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
AppWizard has created this 204101059_yes_no application for you.

This file contains a summary of what you will find in each of the files that
make up your 204101059_yes_no application.


204101059_yes_no.vcxproj
    This is the main project file for VC++ projects generated using an Application Wizard.
    It contains information about the version of Visual C++ that generated the file, and
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

204101059_yes_no.vcxproj.filters
    This is the filters file for VC++ projects generated using an Application Wizard. 
    It contains information about the association between the files in your project 
    and the filters. This association is used in the IDE to show grouping of files with
    similar extensions under a specific node (for e.g. ".cpp" files are associated with the
    "Source Files" filter).

204101059_yes_no.cpp
    This is the main application source file.

/////////////////////////////////////////////////////////////////////////////
Other standard files:

StdAfx.h, StdAfx.cpp
    These files are used to build a precompiled header (PCH) file
    named 204101059_yes_no.pch and a precompiled types file named StdAfx.obj.

Input Files
	This folder contains text files that are used to give input to the program.

STEandZCRfiles
	This folder contains text files of STEs and ZCRs given by program as output.

/////////////////////////////////////////////////////////////////////////////
Other notes:

AppWizard uses "TODO:" comments to indicate parts of the source code you
should add to or customize.

/////////////////////////////////////////////////////////////////////////////
