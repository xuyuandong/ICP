/*
 * main.cpp
 *
 *  Created on: Dec 24, 2008
 *      Author: a
 */

#include "ICP.h"

int main()
{
	ICP myicp(1000,0.0001,50);
	myicp.readfile("1.obj","2.obj");
	myicp.run();
	myicp.writefile("out.obj");
	return 0;
}
