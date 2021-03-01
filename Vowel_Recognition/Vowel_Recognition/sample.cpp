
//***************************************************************VOWEL RECOGNITION ASSIGNMENT***********************************************************************

#include "stdafx.h"
#include<iostream>
#include<conio.h>
#include<fstream>
#include<vector>
#include<string>
#include <sstream>
using namespace std;

#define ld long double
#define li long int
#define dcCheckSamples 500
#define maxNormalizedValue 10000
#define pastSamples_p 12
#define frameSamples 320
#define silenceThreshold 3000


vector<ld> speechSamples;
ld hammingWeights[frameSamples], raisedSinWeights[pastSamples_p + 1], R[pastSamples_p + 1];
ld lpCoefficients_ai[pastSamples_p + 1], cepstrals_Ci[pastSamples_p + 1];
ld tokhuraWeights[] = {0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
int arr[5][pastSamples_p], brr[5][pastSamples_p], testCepstrals[60];
li leaveSamples = 0, totalSamples = 0, fileNumber = 0, fileCount = 0, test = 0;
string readText, fileName, fileType = "reference";				

void dcShift() {

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

void normalization() {

	

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

void hammingWindow() {

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

void raisedSinWindow() {

	

	fstream writeCepstrals;

	if(fileType == "reference") {

		if(fileNumber == 1)
			writeCepstrals.open("STEandZCRfiles\\aaCepstrals.txt", ios::app);						
		else
			writeCepstrals.open("STEandZCRfiles\\eeCepstrals.txt", ios::app);

	}

	
	for(int m = 1; m <= pastSamples_p; m ++) {
		cepstrals_Ci[m] *= raisedSinWeights[m];
		//cout << cepstrals_Ci[m] << "  ";
		if(fileType == "reference") writeCepstrals << fixed << cepstrals_Ci[m] << " ";
		else testCepstrals[test++] = cepstrals_Ci[m];
	}
	if(fileType == "reference") writeCepstrals << endl;
	

	if(fileType == "reference") writeCepstrals.close();

	

	return;
}

void cepstralCoefficients() {

	
	for(int m = 1; m <= pastSamples_p; m ++) {
		lpCoefficients_ai[m] *= -1;
	}

	cepstrals_Ci[0] = 2 * log(R[0]);

	//cout << "Cepstral Coefficients Ci Vector :- " << endl;
	for(int m = 1; m <= pastSamples_p; m ++) {
		cepstrals_Ci[m] = lpCoefficients_ai[m];
		for(int k = 1; k < m; k ++) {
			cepstrals_Ci[m] += ((k * cepstrals_Ci[k] * lpCoefficients_ai[m-k])/m);
		}
		
	}
	

	raisedSinWindow();

	return;
}

void durbinsAlgorithm() {

	

	ld E[pastSamples_p + 1], K[pastSamples_p + 1], sum;
	vector<ld> alpha[pastSamples_p + 1];

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
	 

	

	cepstralCoefficients();

	return;
}

void autoCorrelation() {

	ld maxRi_Energy = 0, maxxRi_index = 0;

	
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

	if(fileType == "test") test = 0;

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

	fstream fileRead, fileEE, writeCepstrals, writeEE;		
	writeCepstrals.open("aaCepstrals.txt", ios::out);
	writeCepstrals.close();
	writeCepstrals.open("eeCepstrals.txt", ios::out);
	writeCepstrals.close();

	fileName = "204101059_aa_";

	while(fileNumber ++ != 2) {
		
		while(fileCount ++ != 10) {

			leaveSamples = 0;
			totalSamples = 0;

			ostringstream temp;
			temp << fileCount;

			fileRead.open(fileName + temp.str() + ".txt", ios::in);		//opening the file for reading the amplitudes of speech waveform

			
			while (fileRead >> readText) {			
				if(leaveSamples++ > 1000) {
					speechSamples.push_back(stod(readText));
					totalSamples ++;
				}
			}

			fileRead.close();

			
			for(li m = 1; m <= pastSamples_p; m ++) {
				raisedSinWeights[m] = 1 + ((pastSamples_p / 2) * sin(3.14 * m / (pastSamples_p)));
				
			}
			
			dcShift();
			normalization();
			hammingWindow();
			autoCorrelation();
		}
		fileCount = 0;
		fileName = "204101059_ee_";

	}

	fileRead.open("aaCepstrals.txt", ios::in);
	fileEE.open("eeCepstrals.txt", ios::in);
	
	for(int i = 0; i < 10; i ++) {
		for(int j = 0; j < 5; j ++) {
			for(int k = 0; k < 12; k ++) {
				fileRead >> readText;
				arr[j][k] += stoi(readText);
				fileEE >> readText;
				brr[j][k] += stoi(readText);
			}
		}
	}

	fileRead.close();
	fileEE.close();

	writeCepstrals.open("STEandZCRfiles\\aaCepstrals.txt", ios::out);
	writeCepstrals.close();
	writeCepstrals.open("STEandZCRfiles\\eeCepstrals.txt", ios::out);
	writeCepstrals.close();

	writeCepstrals.open("aaCepstrals.txt", ios::app);
	writeEE.open("eeCepstrals.txt", ios::out);
	for(int j = 0; j < 5; j ++) {
		for(int k = 0; k < 12; k ++) {
			writeCepstrals << fixed << arr[j][k]/10 << " ";
			writeEE << fixed << brr[j][k]/10 << " ";
		}
		writeCepstrals << endl;
		writeEE << endl;
	}

	writeCepstrals.close();
	writeEE.close();


	fileNumber = 0;
	fileCount = 11;

	while(fileNumber ++ != 2) {
		
		while(fileCount ++ != 20) {

			leaveSamples = 0;
			totalSamples = 0;

			ostringstream temp;
			temp << fileCount;

			fileRead.open(fileName + temp.str() + ".txt", ios::in);	

			

			while (fileRead >> readText) {			
				if(leaveSamples++ > 1000) {
					speechSamples.push_back(stod(readText));
					totalSamples ++;
				}
			}

			fileRead.close();

			
			for(li m = 1; m <= pastSamples_p; m ++) {
				raisedSinWeights[m] = 1 + ((pastSamples_p / 2) * sin(3.14 * m / (pastSamples_p)));
				
			}
			

			dcShift();
			normalization();
			hammingWindow();
			autoCorrelation();

			fileRead.open("aaCepstrals.txt", ios::in);
			fileEE.open("eeCepstrals.txt", ios::in);
			

			ld distAA = 0, distEE = 0;
			int weight = 1;

			for(int j = 0; j < pastSamples_p * 5; j ++) {
				fileRead >> readText;

				distAA += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

				fileEE >> readText;

				distEE += (tokhuraWeights[weight] * ((stod(readText) - testCepstrals[j]) * (stod(readText) - testCepstrals[j])));

				weight ++;
				if(weight == 13) weight = 1;


			}

			cout << fixed << distAA << " " << distEE << endl;

			if(distEE < distAA) cout << "Vowel ee" << endl;
			else if(distEE > distAA) cout << "Vowel aa" << endl;
			else cout << "Ambigous" << endl;

		fileRead.close();
		fileEE.close();

		}
		fileCount = 0;
		fileName = "204101059_ee_";

	}

	return 0;
}

