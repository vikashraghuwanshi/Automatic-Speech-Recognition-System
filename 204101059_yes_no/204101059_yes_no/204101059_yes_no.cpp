

//****************************************** Yes/No Classification Using STE(Short Term Energy) and ZCR (Zero Crossing Rate) *****************************************



//ReadMe file contains idea and other details about files.


#include "stdafx.h"
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<math.h>
using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{


	while(true) {

		//Initialising variables

	string readText, fileName;				//fileName is used to store name of the file and readText is used to store value taken from file
	long double shortTermEnergy = 0;			//shortTermEnergy is used to calculate Short Term Energy for each frame
	vector < long long int > zeroCrossingRates;		//zeroCrossingRates is a vector used to store zero crossing rates for each word
	int sign = 0, flag = 0, maxFlag = 0, kFlag = 0;		//these variables are used as flags
	long long int zcr = 0, zcrCount = 0, sumZCR = 0, trailer = 0, i = 0, j = 0;

			/* zcr is used to calculate zcr values for a frame, it increments by 1 when there is change in signs (i.e. +, -)
			zcrCount is used to count number of zcr values that are greater than some threshold
			sumZCR is used to store addition of zcr values for a single word
			i and j variables used for loops everywhere
			trailer is used on vector of zcr values								
			*/
	
	system("Recording_Module.exe 3 InputFiles/input_file.wav InputFiles/yesnos.txt");

	//Taking filenames and opening files

	cout << "Enter FileName to Recognize Either Yes/No OR 0 to Exit :=" << endl;
	//cin >> fileName;			//filename is stored in the variable fileName
	fileName = "yesnos.txt";

	fstream fileRead, writeSTE, writeZCR;		//fstream variables are used to perform read and write operations on file

	fileRead.open("InputFiles\\" + fileName, ios::in);		//opening the file for reading the amplitudes of speech waveform
	fileName.erase(fileName.length()-4);						//removing .txt from 

	//for storing STE and ZCR files into separate folder "STEandZCRfiles\\" is added while opening
	fileName += "ste.txt";									//adding ste.txt to filname to make new file for storing ZCR values
	writeSTE.open("STEandZCRfiles\\" + fileName, ios::out);						//opening the file for writing STEs
	fileName.erase(fileName.length()-7);						//removing ste.txt from filename
	fileName += "zcr.txt";									//adding zcr.txt to filname to make new file for storing ZCR values
	writeZCR.open("STEandZCRfiles\\" + fileName, ios::out);						//opening the file for writing ZCRs




	//start reading file, calculating ste and zcr and then predicting YES/NO

	fileRead >> readText;					//reading first value of first sample to get the initial sign for zcr calculation
	if(stod(readText) < 0) sign = 1;		//getting the sign as + or - and stod() takes string and convert it to double


	while (fileRead >> readText)			//running a loop until file reached at end
	{ 

		if(stod(readText) < 0 && sign == 0) {		//increment zcr when sign changes from - to +
			zcr += 1;
			sign = 1;
		}
		else if(stod(readText) >= 0 && sign == 1) {		//increment zcr when sign changes from + to -
			zcr += 1;
			sign = 0;
		}

		shortTermEnergy += (pow(stod(readText), 2));	//adding squares of amplitudes to calculate ShortTermEnergy for each frame of 320 samples each
		i += 1;



		if(i == 320) {				//when total samples becomes 320

			writeSTE << fixed << shortTermEnergy/320 << endl;   //writing STE value into ste file

			if((shortTermEnergy / 320) > 4000) flag = 1;		//if(normalized STE of any frame > 30000) make flag as 1, which means its a word boundary from silence
			if((shortTermEnergy / 320) <= 4000) flag = 0;		//if(normalized STE of any frame < 10000) make flag as 0, which means its a silence STE



			if(flag == 1) {			//if flag equals to 1 then we will calculate zcr for word and store them in a vector zeroCrossingRates
				writeZCR << zcr << endl;										//writing zcr values to zcr files
				if(zcr > 100) maxFlag = 1;										//if(zcr > 100) then its a zcr for fricative_sh therefore set maxflag = 1
				if(maxFlag == 1 && zcr < 80) flag = 0;							//if at somepoint zcr becomes less than 80 that means fricative_sh is ended
																						//and silence started because zcr of fricative_sh quite high i.e > 100
				zeroCrossingRates.push_back(zcr);			//pushing zcr of words into vector
				kFlag = 1;				//kFlag is used to tell that the word has already started 
			}



			if(flag == 0 && kFlag == 1) {		//when flag = 0 then word boundary ended and kFlag = 1 means there is a word to predict therefore predict Yes/No


				if(zeroCrossingRates.size() >= 10) {		//if word has atleast 10 frames

					trailer = ceil(zeroCrossingRates.size() * 0.40);		//point to last 40% data

							//in this loop, calculating the sum of last 40% data and count how many zcr values greater than 50
					for(j = zeroCrossingRates.size() - 1; j > zeroCrossingRates.size() - trailer + 1; j --) {
						sumZCR += zeroCrossingRates[j];
						if(zeroCrossingRates[j] > 50) zcrCount ++;
					}
							
					if(zcrCount >= 3 || (sumZCR/trailer) > 80) {		//if avg of last 40% data greater than 80 and atleast 3 Zcr values greater than 50

						cout << "YES" << endl;			//the word is "YES"
					}
					else if(zcrCount <= 2 && (sumZCR/trailer) < 60) {	//if avg less than 60 and atmost 2 values greater than 50 than 

						cout << "NO" << endl;				//the word is "NO"
					}


							
				}
				else {
							//if no. of zcr values less than 10 then check how many values greater than 70
					for(j = 0; j < zeroCrossingRates.size(); j ++) {
						if(zeroCrossingRates[j] > 70) zcrCount ++;
					}


					if(zcrCount >= 4) {				//if atleast 4 values greater than 70 in this much small data then more possibility of having "YES"
						cout << "YES" << endl;
					}


					else if(zeroCrossingRates.size() > 6) {		//if there is data and atmost 3 values greater than 70 then "NO"
						cout << "NO" << endl;
					}

							//but if data has less 6 values then discard that data, there is Ambiguity
				}
					

				//clearing vector after word predicted and putting 0 in all other variables
				zcrCount = 0;
				kFlag = 0;
				sumZCR = 0;
				maxFlag = 0;
				zeroCrossingRates.clear();
			}


			zcr = 0;
			i = 0;
			shortTermEnergy = 0;


					//taking first sample value of next frame to initialize sign for ZCR
			fileRead >> readText;
			if(stod(readText) < 0) sign = 1;
			else sign = 0;

		}
	}

			//if last frame is less than 320 samples then to store STE and ZCR

	if(shortTermEnergy != 0) writeSTE << shortTermEnergy/320 << endl;
	if(zcr != 0) writeZCR << zcr << endl;



			//if last frame has sufficient data to predict a word
	if(zeroCrossingRates.size() > 6) {				//if avg of last 40% data greater than 80 and atleast 3 Zcr values greater than 50
		if(zcrCount >= 3 || (sumZCR/trailer) > 80) {
			cout << "YES" << endl;
		}
		else if(zcrCount <= 2 && (sumZCR/trailer) < 60) {	//if avg less than 60 and atmost 2 values greater than 50 than 
			cout << "NO" << endl;
		}				
	}

			
			//closing all the files opened for reading and writing
	fileRead.close();
	writeSTE.close();
	writeZCR.close();

	char c;
	cout << "Press 0 to exit OR any character other than 0 to continue..." << endl;
	cin >> c;

	if(c == '0') exit(0);

	}

	return 0;
}

