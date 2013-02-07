/*
 * ICP.cpp
 *
 *  Created on: Dec 24, 2008
 *      Author: a
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <math.h>
#include <time.h>
#include <newmat/newmat.h>
#include <newmat/newmatap.h>
#include "ICP.h"

using namespace std;
using namespace NEWMAT;

ICP::ICP(int controlnum,double thre,int iter)
{
	cono = controlnum;
	threshold = thre;
	iterate = iter;
	//
	contP = new Vertex[cono];
	assert(contP!=NULL);
	contQ = new Vertex[cono];
	assert(contQ!=NULL);
	rmcoP = new Vertex[cono];
	assert(rmcoP!=NULL);
	rmcoQ = new Vertex[cono];
	assert(rmcoQ!=NULL);
	index = new int[cono];
	assert(index!=NULL);
}

ICP::~ICP()
{
	delete [] contP;
	delete [] contQ;
	delete [] rmcoP;
	delete [] rmcoQ;
	delete [] index;
}

void ICP::readfile(string firstname,string secondname)
{
	cout<<"read two clouds of points from obj files"<<endl;
	char buf[1024];
	ifstream inobj;

	inobj.open(firstname.c_str());
	while(inobj.getline(buf,sizeof(buf)))
	{
		string s = buf;
		if(s.find_first_of("v") == 0)
		{
			istringstream is(s);
			string title;
			Vertex v;
			is>>title>>v.coord[0]>>v.coord[1]>>v.coord[2];
			VarrP.push_back(v);
		}
	}
	inobj.close();
	//
	inobj.open(secondname.c_str());
	while(inobj.getline(buf,sizeof(buf)))
	{
		string s = buf;
		if(s.find_first_of("v") == 0)
		{
			istringstream is(s);
			string title;
			Vertex v;
			is>>title>>v.coord[0]>>v.coord[1]>>v.coord[2];
			VarrQ.push_back(v);
		}
	}
	inobj.close();
}

void ICP::run()
{
	initransmat();
	sample();
	//
	double err = closest();
	cout<<"initial error = "<<err<<endl;
	//
	for(int i=0;i<iterate;i++)
	{
		getcenter();
		rmcontcenter();
		transform();
		uprotate();
		uptranslate();
		updata();
		double newerr = closest();
		cout<<"iterate times = "<<i<<endl;
		cout<<"error = "<<newerr<<endl;
		double delta = fabs(err-newerr)/cono;
		cout<<"delta = "<<delta<<endl;
		if(delta<threshold)
				break;
		err = newerr;

	}
	printTR();
	applyall();
}

void ICP::initransmat()
{
	//cout<<"initial translate and rotate matrix"<<endl;
	//
	for(int i=0;i<3;i++)
		TT[i] = 0.0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			if(i!=j)
				TR[i][j]=0.0;
			else
				TR[i][j]=1.0;
		}
	}
}

void ICP::sample()
{
	//cout<<"sample control points from P"<<endl;
	//
	int N = VarrP.size();
	bool *flag = new bool[N];
	assert(flag!=NULL);
	for(int i=0;i<N;i++)
		flag[i] = false;
	//sample control points index
	srand((unsigned)time(NULL));
	for(int i=0;i<cono;i++)
	{
		while(true)
		{
			int sam =rand()%N;
			if(!flag[sam])
			{
				index[i] = sam;
				flag[sam] = true;
				break;
			}
		}
	}
	//
	//cout<<"store control points into contP"<<endl;
	//
	for(int i=0;i<cono;i++)
	{
		Vertex v = VarrP[index[i]];
		for(int j=0;j<3;j++)
		{
			contP[i].coord[j] = v.coord[j];
		}
	}

	delete [] flag;
}

double ICP::closest()
{
	//cout<<"find closest points and error"<<endl;
	//
	double error = 0.0;
	for(int i=0;i<cono;i++)
	{
		double maxdist = 100.0;
		index[i] = 0;
		for(unsigned int j=0;j<VarrQ.size();j++)
		{
			double dist = distance(contP[i],VarrQ[j]);
			if(dist < maxdist)
			{
				maxdist = dist;
				index[i] = j;
			}
		}
		Vertex v = VarrQ[index[i]];
		for(int j=0;j<3;j++)
		{
			contQ[i].coord[j] = v.coord[j];
		}
		error += maxdist;
	}
	return error;
}

void ICP::getcenter()
{
	//cout<<"get center from two clouds of control points"<<endl;
	//
	for(int i=0;i<3;i++)
		meanP.coord[i] = 0.0 ;
	for(int i=0;i<cono;i++)
	{
		for(int j=0;j<3;j++)
			meanP.coord[j] += contP[i].coord[j];
	}
	for(int i=0;i<3;i++)
		meanP.coord[i] = meanP.coord[i]/cono;
	//
	for(int i=0;i<3;i++)
		meanQ.coord[i] = 0.0;
	for(int i=0;i<cono;i++)
	{
		for(int j=0;j<3;j++)
			meanQ.coord[j] += contQ[i].coord[j];
	}
	for(int i=0;i<3;i++)
		meanQ.coord[i] = meanQ.coord[i]/cono;
}

void ICP::rmcontcenter()
{
	cout<<"move clouds of control points to their correspond points center"<<endl;
	//
	for(int i=0;i<cono;i++)
	{
		for(int j=0;j<3;j++)
		{
			rmcoP[i].coord[j] = contP[i].coord[j] - meanP.coord[j];
		}
	}
	//
	for(int i=0;i<cono;i++)
	{
		for(int j=0;j<3;j++)
		{
			rmcoQ[i].coord[j] = contQ[i].coord[j] - meanQ.coord[j];
		}
	}
}

void ICP::transform()
{
	cout<<"get transform matrix"<<endl;
	//
	Matrix B(4,4);
	B = 0;
	double u[3]; //di+di'
	double d[3]; //di-di'
	for(int i=0;i<cono;i++)
	{
		for(int j=0;j<3;j++)
		{
			u[j] = rmcoP[i].coord[j]+rmcoQ[i].coord[j];
			d[j] = rmcoP[i].coord[j]-rmcoQ[i].coord[j];
		}
		double uM[16] = {0,    -d[0],  -d[1],  -d[2],
						 d[0], 0,      -u[2],  -u[1],
						 d[1], -u[2],  0,      u[0],
						 d[2], u[1],   -u[0],  0   };
		/*
		double uM[16] = {0,    -u[2] , -u[1], d[0],
						-u[2], 0    ,  u[0], d[1],
						u[1] , -u[0],   0  , d[2],
						-d[0], -d[1], -d[2],  0 };*/
		Matrix Ai(4,4) ;
		Ai << uM ;
		B += Ai*Ai.t();
	}
	//
	Matrix U;
	Matrix V;
	DiagonalMatrix D;
	SVD(B,D,U,V);
	//

	for(int i=0;i<4;i++)
	{
		quad[i] = V.element(i,3);
	}
	//
	B.Release();
	U.Release();
	V.Release();
	D.Release();
}

void ICP::uprotate()
{
	//cout<<"change quadternious to 3*3 rotate matrix"<<endl;
	//
	Rw[0][0] = quad[0]*quad[0]+quad[1]*quad[1]-quad[2]*quad[2]-quad[3]*quad[3];
	Rw[0][1] = 2*(-quad[0]*quad[3]+quad[1]*quad[2]);
	Rw[0][2] = 2*(quad[0]*quad[2]+quad[1]*quad[3]);
	Rw[1][0] = 2*(quad[0]*quad[3]+quad[1]*quad[2]);
	Rw[1][1] = quad[0]*quad[0]-quad[1]*quad[1]+quad[2]*quad[2]-quad[3]*quad[3];
	Rw[1][2] = 2*(-quad[0]*quad[1]+quad[2]*quad[3]);
	Rw[2][0] = 2*(-quad[0]*quad[2]+quad[1]*quad[3]);
	Rw[2][1] = 2*(quad[0]*quad[1]+quad[2]*quad[3]);
	Rw[2][2] = quad[0]*quad[0]-quad[1]*quad[1]-quad[2]*quad[2]+quad[3]*quad[3];
	//
	//Rn+1 = R * Rn
	double tmp[3][3];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			tmp[i][j] = 0.0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				tmp[i][j] += Rw[i][k]*TR[k][j];
			}
		}
	}
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			TR[i][j] = tmp[i][j];
}

void ICP::uptranslate()
{
	//Tw = P' -Rw*P
	double tmp[3] = {0,0,0};
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			tmp[i] += Rw[i][j]*meanP.coord[j];
		}
	}
	for(int i=0;i<3;i++)
	{
		Tw[i] = meanQ.coord[i] - tmp[i];
	}
	//Tn+1 = R*Tn + Tw
	double temp[3] = {0,0,0};
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			temp[i] += Rw[i][j]*TT[j];
		}
	}
	for(int i=0;i<3;i++)
	{
		TT[i] = temp[i] + Tw[i];
	}

}

void ICP::updata()
{
	//cout<<"update control points in P"<<endl;
	//
	//rotate + translate
	for(int i=0;i<cono;i++)
	{
		double tmp[3] = {0,0,0};
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				tmp[j] += Rw[j][k]*contP[i].coord[k];
			}
		}
		for(int j=0;j<3;j++)
			contP[i].coord[j] = tmp[j] + Tw[j];

	}
}

void ICP::applyall()
{
	cout<<"transform to all data in P"<<endl;
	// make rotate
	for(vector<Vertex>::iterator p=VarrP.begin();p!=VarrP.end();p++)
	{
		Vertex v = *p;
		double tmp[3] = {0,0,0};
		for(int i=0;i<3;i++)
		{
			for(int k=0;k<3;k++)
			{
				tmp[i] += TR[i][k]*v.coord[k];
			}
		}
		for(int i=0;i<3;i++)
		{
			v.coord[i] = tmp[i] + TT[i];
		}
		*p = v;
	}
	//
}

void ICP::writefile(string name)
{
	cout<<"output clouds of points P after transform"<<endl;
	//
	ofstream outobj;
	outobj.open(name.c_str());
	outobj<<"# Geomagic Studio"<<endl;
	int num = 1;
	for(vector<Vertex>::const_iterator p=VarrP.begin();p!=VarrP.end();p++)
	{
		Vertex v;
		v = *p;
		outobj<<"v "<<v.coord[0]<<" "<<v.coord[1]<<" "<<v.coord[2]<<endl;
		outobj<<"p "<<num++<<endl;
	}
	//
	outobj.close();
}

double ICP::distance(Vertex a,Vertex b)
{
	double dist = 0.0;
	for(int i=0;i<3;i++)
	{
		dist += (a.coord[i]-b.coord[i])*(a.coord[i]-b.coord[i]);
	}
	return dist;
}

void ICP::printTT()
{
	cout<<"Translate Matrix = "<<endl;
	for(int i=0;i<3;i++)
	{
		cout<<TT[i]<<"  ";
	}
	cout<<endl;
}

void ICP::printTR()
{
	cout<<"Rotate Matrix = "<<endl;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			cout<<TR[i][j]<<" ";
		}
		cout<<endl;
	}
}
