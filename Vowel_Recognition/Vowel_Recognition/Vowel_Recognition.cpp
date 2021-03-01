
//***************************************************************VOWEL RECOGNITION ASSIGNMENT***********************************************************************

#include "stdafx.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include <sstream>
using namespace std;

#define ld long double
#define li long int
#define dcCheckSamples 500			//indicates number of samples used for DC SHIFT
#define maxNormalizedValue 10000	//indicates maximum value of normalization
#define pastSamples_p 12			//indicates past speech samples used
#define frameSamples 320			//indicates samples in a frame
#define silenceThreshold 3000		//indicates threshold for silence noise boundary


vector<ld> speechSamples; //used to store amplitudes of speech waveform
ld hammingWeights[frameSamples], raisedSinWeights[pastSamples_p + 1], R[pastSamples_p + 1];
							/*hammingWeights used to store weights of hamming window
							 raisedSinWeights used to store weights of raised Sine window
							 R is used to store autocorrelation vectors*/

ld lpCoefficients_ai[pastSamples_p + 1], cepstrals_Ci[pastSamples_p + 1], testCepstrals[60];
							/*lpCoefficients_ai used to store our linear predictive coefficients
							  cepstrals_Ci used to store cepstral coefficients
							  testCepstrals used to store cepstral coefficients of 5 steady frames i.e total 60 values of test files
							  */

ld tokhuraWeights[] = {0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
							//tokhuraWeights are used to store the tokhura weights for reference files

ld aCepstrals[5][pastSamples_p], eCepstrals[5][pastSamples_p], iCepstrals[5][pastSamples_p], oCepstrals[5][pastSamples_p], uCepstrals[5][pastSamples_p];
							// these arrays are used to store cepstral references for steady frames

li leaveSamples = 0, totalSamples = 0, fileNumber = 0, fileCount = 0, test = 0, accuracy = 0;
							/*leaveSamples used to store the count of samples that are leaved initially.
							  totalSamples used to store the count of total samples in the file.
							  fileNumber is used to store the number of file i.e. a-1, e-2, i-3, o-4, u-5.
							  fileCount is used to store count of file of particular vowel.
							  test is used to store the count of cepstral number in the test files.
							  accuracy is used to count number of correct predictions*/

string readText, fileName, fileType = "reference", rollNumber="204101059";
							/*fileName is used to store name of the file and readText is used to store value taken from file
							  fileType indicates reference file or test file
							  rollNumber is used to give initial name to filename
							  readText is used to store read text from the file*/


void dcShift() {	//applying dc shift by taking some initial samples

	ld samplesAverage = 0;

	for(li m = 0; m < dcCheckSamples; m ++) {
		samplesAverage += speechSamples[m];
	}

	samplesAverage = (samplesAverage/dcCheckSamples);

	if(samplesAverage == 0) return;

	for(li i = 0; i < totalSamples; i ++) {
		speechSamples[i] -= samplesAverage;
	}

	return;
}

void normalization() {	//normalising the speech amplitudes by converting max value to maxNormalized value

	ld maxSampleValue = 0;

	for(li i= 0; i < totalSamples; i ++) {
		if(abs(speechSamples[i]) > maxSampleValue) {
			maxSampleValue = abs(speechSamples[i]);
		}
	}

	maxSampleValue = maxNormalizedValue/maxSampleValue;

	for(li m = 0; m < totalSamples; m ++) {
		speechSamples[m] *= maxSampleValue;
	}

	return;
}

void hammingWindow() {		//applying hamming window to normalized input by taking sizes of each frame

	for(li m = 0; m < frameSamples; m ++) {
		hammingWeights[m] = 0.54 - (0.46 * cos(2 * 3.14 * m / (frameSamples - 1)));
	}

	for(li n = 0; n < totalSamples/frameSamples; n ++) {
		for(int m = 0; m < frameSamples; m ++) {
			speechSamples[m + (n * frameSamples)] *= hammingWeights[m];
		}
	}

	return;
}

void raisedSinWindow() {		//Applying Raised Sin Window on Cepstral Coefficients and storing those coefficients into text files

	fstream writeCepstrals;

			//if fileType is reference then only we need to store in text files for later use and 
			//if fileType is test then we only need to calculate at real time to check accuracy

	if(fileType == "reference") {

										//opening the file for writing Cepstral Coefficients

		if(fileNumber == 1)
			writeCepstrals.open("Cepstral_Files\\aCepstrals.txt", ios::app);						
		else if(fileNumber == 2)
			writeCepstrals.open("Cepstral_Files\\eCepstrals.txt", ios::app);
		else if(fileNumber == 3)
			writeCepstrals.open("Cepstral_Files\\iCepstrals.txt", ios::app);
		else if(fileNumber == 4)
			writeCepstrals.open("Cepstral_Files\\oCepstrals.txt", ios::app);
		else
			writeCepstrals.open("Cepstral_Files\\uCepstrals.txt", ios::app);

	}

	for(int m = 1; m <= pastSamples_p; m ++) {
		cepstrals_Ci[m] *= raisedSinWeights[m];
		if(fileType == "reference") writeCepstrals << fixed << cepstrals_Ci[m] << " ";	//writing into file if refercing is done
		else testCepstrals[test++] = cepstrals_Ci[m];									//calcuting raised ci's not writing in file
	}
	if(fileType == "reference") {
		writeCepstrals << endl;
		writeCepstrals.close();
	}

	return;
}

void cepstralCoefficients() {		//computing cepstral coefficients by using LP coefficients returned by durbins algorithm

	//for(int m = 1; m <= pastSamples_p; m ++) {	//inverting the ai's before calculating ci's
	//	lpCoefficients_ai[m] *= -1;
	//}

	cepstrals_Ci[0] = 2 * log(R[0]);

	for(int m = 1; m <= pastSamples_p; m ++) {
		cepstrals_Ci[m] = lpCoefficients_ai[m];
		for(int k = 1; k < m; k ++) {
			cepstrals_Ci[m] += ((k * cepstrals_Ci[k] * lpCoefficients_ai[m-k])/m);
		}
	}

	raisedSinWindow();		//after calculating cepstral coefficients calling raisedSinWindow() to apply raised Sine Window

	return;
}

void durbinsAlgorithm() {		//computing Linear Predictive coefficients by using autoCorrelation Vector returned by AutoCorrelation method algorithm

	ld E[pastSamples_p + 1], K[pastSamples_p + 1], sum;
	vector<ld> alpha[pastSamples_p + 1];	//making an array of vectors for storing alphas for diffent i values

	E[0] = R[0];

	for(int i = 1; i <= pastSamples_p; i ++) {

		alpha[i].push_back(0);
		
		sum = 0.0;
		
		for(int j = 1; j <= i - 1; j ++) {
			sum += (alpha[i-1][j] * R[i-j]);
		}
		
		K[i] = (R[i] - sum)/E[i-1];

		for(int j = 1; j <= i-1; j ++) {
			alpha[i].push_back(alpha[i-1][j] - (K[i] * alpha[i-1][i-j]));
		}

		alpha[i].push_back(K[i]);

		E[i] = (1 - (K[i] * K[i])) * E[i-1];

	}

	for(int i = 1; i <= pastSamples_p; i ++) {
		lpCoefficients_ai[i] = alpha[pastSamples_p][i];
	}

	cepstralCoefficients();		//calling cepstralCoefficients() to calculate Cepstral Coefficients

	return;
}

void autoCorrelation() {

	ld maxRi_Energy = 0, maxxRi_index = 0;

									//this loop is used to find the index of frame which has highest autocorrelation with lag zero i.e R[0]

	for(int n = 0; n < totalSamples/frameSamples; n ++) {
		for(int k = 0; k <= pastSamples_p; k ++) {
			R[k] = 0;
			for(int m = 0; m < frameSamples - k; m ++) {
				R[k] += (speechSamples[m + (frameSamples * n)] * speechSamples[m + (frameSamples * n) + k]);
			} 

			if(R[0] > maxRi_Energy) {
				maxRi_Energy = R[0];
				maxxRi_index = n;
			}

		}
		
	}
					//if test file then count which cepstral number is calculating
	if(fileType == "test") test = 0;

					//Now by index use 5 steady frames 2 before 2 after 1 itself apply durbins on them
	for(int n = maxxRi_index - 2; n < maxxRi_index + 3; n ++) {
		for(int k = 0; k <= pastSamples_p; k ++) {
			R[k] = 0;
			for(int m = 0; m < frameSamples - k; m ++) {
				R[k] += (speechSamples[m + (frameSamples * n)] * speechSamples[m + (frameSamples * n) + k]);
			} 
		}

		if(R[0] >= silenceThreshold) {
			durbinsAlgorithm();
		}
		
	}
	return;
}



int _tmain(int argc, _TCHAR* argv[])
{	

	fstream fileRead, fileEE, fileII, fileOO, fileUU, writeCepstrals, writeEE, writeII, writeOO, writeUU;		//fstream variables are used to perform read and write operations on file

		//opening and closing all reference files to make them empty

	writeCepstrals.open("Cepstral_Files\\aCepstrals.txt", ios::out);
	writeCepstrals.close();
	writeCepstrals.open("Cepstral_Files\\eCepstrals.txt", ios::out);
	writeCepstrals.close();
	writeCepstrals.open("Cepstral_Files\\iCepstrals.txt", ios::out);
	writeCepstrals.close();
	writeCepstrals.open("Cepstral_Files\\oCepstrals.txt", ios::out);
	writeCepstrals.close();
	writeCepstrals.open("Cepstral_Files\\uCepstrals.txt", ios::out);
	writeCepstrals.close();

	fileName = rollNumber + "_a_";

		//iterating through all the type of files to calculate reference files

	while(fileNumber ++ != 5) {

				//for each type of file using first 10 files to calculate reference files
		
		while(fileCount ++ != 10) {

			speechSamples.assign(0, totalSamples);

			leaveSamples = 0;
			totalSamples = 0;

			ostringstream temp;
			temp << fileCount;		//making fileCount a string

			fileRead.open("InputFiles\\" + fileName + temp.str() + ".txt", ios::in);		//opening the file for reading the amplitudes of speech waveform

			while (fileRead >> readText) {			//running a loop until file reached at end and store in vector after leaving some initial samples
				if(leaveSamples++ > 100 ) {
					speechSamples.push_back(stod(readText));
					totalSamples ++;
				}
			}

			fileRead.close();

			//cout << "Total Samples :- " << totalSamples << endl << endl;

								//calculating raised sine weights
			for(li m = 1; m <= pastSamples_p; m ++) {
				raisedSinWeights[m] = 1 + ((pastSamples_p / 2) * sin(3.14 * m / (pastSamples_p)));
			}

			dcShift();
			normalization();
			hammingWindow();
			autoCorrelation();

			speechSamples.clear();

		}
		fileCount = 0;
							//after one type file change to another type file
		if(fileNumber == 1) fileName = rollNumber + "_e_";
		else if(fileNumber == 2) fileName = rollNumber + "_i_";
		else if(fileNumber == 3) fileName = rollNumber + "_o_";
		else fileName = rollNumber + "_u_";

	}

				//opening files for reading cepstral coefficients
	fileRead.open("Cepstral_Files\\aCepstrals.txt", ios::in);
	fileEE.open("Cepstral_Files\\eCepstrals.txt", ios::in);
	fileII.open("Cepstral_Files\\iCepstrals.txt", ios::in);
	fileOO.open("Cepstral_Files\\oCepstrals.txt", ios::in);
	fileUU.open("Cepstral_Files\\uCepstrals.txt", ios::in);
	
				//calculating sums of all cepstral coefficients with respect to each frame and ci's
	for(int i = 0; i < 10; i ++) {
		for(int j = 0; j < 5; j ++) {
			for(int k = 0; k < 12; k ++) {
				fileRead >> readText;
				aCepstrals[j][k] += stoi(readText);
				fileEE >> readText;
				eCepstrals[j][k] += stoi(readText);
				fileII >> readText;
				iCepstrals[j][k] += stoi(readText);
				fileOO >> readText;
				oCepstrals[j][k] += stoi(readText);
				fileUU >> readText;
				uCepstrals[j][k] += stoi(readText);
			}
		}
	}

	fileRead.close();
	fileEE.close();
	fileII.close();
	fileOO.close();
	fileUU.close();

		//opening files to write final cepstral coefficients

	writeCepstrals.open("Cepstral_Files\\aCepstrals.txt", ios::out);
	writeEE.open("Cepstral_Files\\eCepstrals.txt", ios::out);
	writeII.open("Cepstral_Files\\iCepstrals.txt", ios::out);
	writeOO.open("Cepstral_Files\\oCepstrals.txt", ios::out);
	writeUU.open("Cepstral_Files\\uCepstrals.txt", ios::out);


		//write averages of final cepstral coefficients into files
	for(int j = 0; j < 5; j ++) {
		for(int k = 0; k < 12; k ++) {
			writeCepstrals << fixed << aCepstrals[j][k]/10 << " ";
			writeEE << fixed << eCepstrals[j][k]/10 << " ";
			writeII << fixed << iCepstrals[j][k]/10 << " ";
			writeOO << fixed << oCepstrals[j][k]/10 << " ";
			writeUU << fixed << uCepstrals[j][k]/10 << " ";
		}
		writeCepstrals << endl;
		writeEE << endl;
		writeII << endl;
		writeOO << endl;
		writeUU << endl;
	}

	writeCepstrals.close();
	writeEE.close();
	writeII.close();
	writeOO.close();
	writeUU.close();


	//test file calculations based on 10 test files of each type

	fileNumber = 0;
	fileCount = 10;
	fileType = "test";
	fileName = rollNumber + "_a_";

	while(fileNumber ++ != 5) {
		
		while(fileCount ++ != 20) {

			leaveSamples = 0;
			totalSamples = 0;

			ostringstream temp;
			temp << fileCount;

			fileRead.open("InputFiles\\" + fileName + temp.str() + ".txt", ios::in);		//opening the file for reading the amplitudes of speech waveform

			while (fileRead >> readText) {			//running a loop until file reached at end
				if(leaveSamples++ > 100) {
					speechSamples.push_back(stod(readText));
					totalSamples ++;
				}
			}

			fileRead.close();

			//cout << "Total Samples :- " << totalSamples << endl << endl;

			dcShift();
			normalization();
			hammingWindow();
			autoCorrelation();

						//opening files for reading cepstral coefficients
			fileRead.open("Cepstral_Files\\aCepstrals.txt", ios::in);
			fileEE.open("Cepstral_Files\\eCepstrals.txt", ios::in);
			fileII.open("Cepstral_Files\\iCepstrals.txt", ios::in);
			fileOO.open("Cepstral_Files\\oCepstrals.txt", ios::in);
			fileUU.open("Cepstral_Files\\uCepstrals.txt", ios::in);

					//calculating tokhura distaces for a, e, i, o, u vowels by given tokhura weights
			ld distAA = 0, distEE = 0, distII = 0, distOO = 0, distUU = 0;
			int weight = 1;

			for(int j = 0; j < pastSamples_p * 5; j ++) {

				fileRead >> readText;

				distAA += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

				fileEE >> readText;

				distEE += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

				fileII >> readText;

				distII += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

				fileOO >> readText;

				distOO += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

				fileUU >> readText;

				distUU += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

				weight ++;
				if(weight == 13) weight = 1;


			}

						//finding minimum distance from all the distances
			int index = 0;
			ld findMin[] = {distAA, distEE, distII, distOO, distUU}, minValue = 100000000;

			for(int i = 0; i < 5; i ++) {
				if(minValue > findMin[i]) {
					minValue = findMin[i];
					index = i;

				}
			}

				//printing minimum distance vowel

			for(int i = 0; i < 5; i ++) {
				if(minValue == findMin[i]) {
					if(i == 0) {
						cout << "a" << endl;
						if(fileNumber == 1) accuracy ++;
					}
					else if(i == 1) {
						cout << "e" << endl;
						if(fileNumber == 2) accuracy ++;
					}
					else if(i == 2) {
						cout << "i" << endl;
						if(fileNumber == 3) accuracy ++;
					}
					else if(i == 3) {
						cout << "o" << endl;
						if(fileNumber == 4) accuracy ++;
					}
					else {
						cout << "u" << endl;
						if(fileNumber == 5) accuracy ++;
					}
					break;
				}
		}

		fileRead.close();
		fileEE.close();
		fileII.close();
		fileOO.close();
		fileUU.close();

		speechSamples.clear();

		}
		cout << endl << "*************" << endl << endl;

				//similarly testing for another vowels

		fileCount = 10;
		if(fileNumber == 1) fileName = rollNumber + "_e_";
		else if(fileNumber == 2) fileName = rollNumber + "_i_";
		else if(fileNumber == 3) fileName = rollNumber + "_o_";
		else fileName = rollNumber + "_u_";

	}

	cout << endl << accuracy*2 << "%" << endl;	//printing accuracy

	//live recording testing
	//for vowel e little jerk with high pitch is required

	while(true) {

		cout << endl << "Press 0 to exit or Any Key to Record" << endl << endl;

		system("Recording_Module.exe 3 InputFiles/input_file.wav InputFiles/input_file.txt");

		fileName = "input_file.txt";
		leaveSamples = 0;
		totalSamples = 0;

		ostringstream temp;
		temp << fileCount;

		fileRead.open("InputFiles\\" + fileName, ios::in);		//opening the file for reading the amplitudes of speech waveform

		while (fileRead >> readText) {			//running a loop until file reached at end
			if(leaveSamples++ > 100) {
				speechSamples.push_back(stod(readText));
				totalSamples ++;
			}
		}

		fileRead.close();

		//cout << "Total Samples :- " << totalSamples << endl << endl;

		dcShift();
		normalization();
		hammingWindow();
		autoCorrelation();

		fileRead.open("Cepstral_Files\\aCepstrals.txt", ios::in);
		fileEE.open("Cepstral_Files\\eCepstrals.txt", ios::in);
		fileII.open("Cepstral_Files\\iCepstrals.txt", ios::in);
		fileOO.open("Cepstral_Files\\oCepstrals.txt", ios::in);
		fileUU.open("Cepstral_Files\\uCepstrals.txt", ios::in);

		//calculating tokhura distaces for a, e, i, o, u vowels by given tokhura weights
		ld distAA = 0, distEE = 0, distII = 0, distOO = 0, distUU = 0;
		int weight = 1;

		for(int j = 0; j < pastSamples_p * 5; j ++) {

			fileRead >> readText;

			distAA += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

			fileEE >> readText;

			distEE += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

			fileII >> readText;

			distII += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

			fileOO >> readText;

			distOO += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

			fileUU >> readText;

			distUU += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

			weight ++;
			if(weight == 13) weight = 1;


		}

					//finding minimum distance from all the distances
		int index = 0;
		ld findMin[] = {distAA, distEE, distII, distOO, distUU}, minValue = 100000000;

		for(int i = 0; i < 5; i ++) {
			if(minValue > findMin[i]) {
				minValue = findMin[i];
				index = i;

			}
		}

				//printing minimum distance vowel

		for(int i = 0; i < 5; i ++) {
			if(minValue == findMin[i]) {
				if(i == 0) {
					cout << "a" << endl;
				}
				else if(i == 1) {
					cout << "e" << endl;
				}
				else if(i == 2) {
					cout << "i" << endl;
				}
				else if(i == 3) {
					cout << "o" << endl;
				}
				else {
					cout << "u" << endl;
				}
				break;
			}
		}

		fileRead.close();
		fileEE.close();
		fileII.close();
		fileOO.close();
		fileUU.close();

		speechSamples.clear();	//clearing vector for next iteration

	}

	return 0;
}

