#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>

using namespace std;

void read_pgm(string fname, int &height, int &width, int &gmax, unsigned short **gval){
	ifstream file(fname.c_str());
	string str;
	char arr[128];
	char *pstr;

  int linecount = 0;
	//int idx;
	height = 128;
	pstr = arr;
	while (file.getline(pstr,height)) {
		if(linecount == 0)
    {
			if(arr[0] != 'P' || arr[1] != '5')
			{
				cout << "Error in read_pgm: not a pgm file" << endl;
				return;
			}
			char *p;
			for(p = arr;(*p!=' ') && (*p!=0);p++)
				;
			if(*p == ' ')
			{
				linecount = 1;
				pstr = p + 1;
			}
    }
    if(linecount == 1)
    {
			char *p;
			for(p = pstr;*p!=' ';p++)
				;
			*p = 0;
			height = atoi(pstr);
			width = atoi(p + 1);
//			cout << height << endl;
//			cout << width << endl;
			*gval = new unsigned short[height * width];
			for(p = p + 1;(*p!=' ') && (*p!=0);p++)
				;
			if(*p == ' ')
			{
				linecount = 2;
				pstr = p + 1;
			}
    }
    if(linecount == 2)
    {
			gmax = atoi(pstr);
			if(gmax <= 255)
				file.read((char*)*gval, height * width);
			else if(gmax <= 65535)
			{
				file.read((char*)*gval, height * width * 2);
				for(int i = 0;i < height * width;i++)
				{
					unsigned short g = (*gval)[i], g0;
					g0 = (g >> 8) & 0xff;
					g0 = g0 | (((g & 0xff) << 8) & 0xff00);
					(*gval)[i] = g0;
				}
			}
			else
			{
				cout << "gmax error" << endl;
				return;
			}
			break;
		}

	}
}
