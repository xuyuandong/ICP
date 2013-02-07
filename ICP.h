/*
 * ICP.h
 *
 *  Created on: Dec 24, 2008
 *      Author: a
 */

#ifndef ICP_H_
#define ICP_H_

#include <vector>
#include <string>

using namespace std;

struct Vertex{
	double coord[3];
};

class ICP
{
public:
	ICP(int controlnum = 1000,double thre = 0.01,int iter = 10);
	~ICP();
	void readfile(string firstname,string secondname);
	void run();
	void writefile(string name);

private:
	void initransmat(); //init translate and rotate matrix
	void sample(); //sample control points
	double closest(); // find corresponding points and return error
	void getcenter(); //get center from two control points
	void rmcontcenter(); //remove center from two control points
	void transform(); //get transform (rotate) matrix
	void uprotate(); //change quaternious number to matrix and update whole rotate matrix
	void uptranslate(); // update whole translate matrix
	void updata();  // update control points coordinate
	void applyall();

private:
	double distance(Vertex a,Vertex b);
	void printTT();
	void printTR();

private:
	int cono; // control points number
	int iterate; //iterate number
	double threshold; //stop threshold
	vector<Vertex> VarrP; //original points
	vector<Vertex> VarrQ;
	Vertex meanP; // control points center
	Vertex meanQ;
	Vertex *contP; //control points in P
	Vertex *contQ;
	Vertex *rmcoP; //control points after removing center
	Vertex *rmcoQ;
	int *index;	//use when sampling control points and in finding corresponding points index
	double TT[3]; //translate
	double TR[3][3]; //rotate
	double Rw[3][3];//step rotate
	double Tw[3]; //step translate
	double quad[4];



};

#endif /* ICP_H_ */
