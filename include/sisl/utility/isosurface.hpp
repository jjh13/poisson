#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_

#include <sisl/sisl.hpp>
#include <sisl/utility/ply_writer.hpp>

#include <tuple>
#include <unordered_map>
#include <map>

namespace sisl{
namespace utility{

template<class T>
class marchingCubes{
public:
	template<class L, class I, class O>
	void marchLattice(
			L *f,
			L *dx, 
			L *dy, 
			L *dz,
			const O &isoValue, 
			const I &stepSize,
			sisl::vector3<I> origin,
			sisl::vector3<I> boundary ){

		// Namespaces
		using namespace std;
		using namespace sisl;

		int xres = (int)floor((boundary.i - origin.i)/stepSize);
		int yres = (int)floor((boundary.j - origin.j)/stepSize);
		int zres = (int)floor((boundary.k - origin.k)/stepSize);
		int res = xres*2;

		map<int, sisl::vertex3<T> > vMap;

		#pragma omp parallel for
		for(int i = 0; i < xres; i++) {
			float x = origin.i + i * stepSize;
			int j = 0;

			O 	v000 = 0, v001 = 0, v010 = 0, v100 = 0, 
				v011 = 0, v110 = 0, v101 = 0, v111 = 0;
			// Break up the inner loops over
			for(float y = origin.j; y < boundary.j; y+=stepSize,j++) {
				int idx = 0;
				int k = 0;
				// Eval these only once
				v000 = f->f(x,				y,			origin.k);
				v100 = f->f(x+stepSize, 	y, 			origin.k);
				v110 = f->f(x+stepSize, 	y+stepSize,	origin.k);
				v010 = f->f(x,				y+stepSize,	origin.k);

				idx |=	(v000 < isoValue) ? 16 : 0; 	// {0,0,0}, 0001000
				idx |=	(v100 < isoValue) ? 32 : 0; 	// {1,0,0}, 0010000
				idx |=	(v110 < isoValue) ? 64 : 0; 	// {1,1,0}, 0100000
				idx |=	(v010 < isoValue) ? 128 : 0; 	// {0,1,0}, 1000000

				for(float z = origin.k; z < boundary.k; z+=stepSize, k++) {
					v001 = f->f(x,				y, 			z+stepSize);
					v101 = f->f(x+stepSize,		y, 			z+stepSize);
					v111 = f->f(x+stepSize,		y+stepSize,	z+stepSize);
					v011 = f->f(x,				y+stepSize,	z+stepSize);

					idx >>= 4;
					idx |= (v001 < isoValue) ? 16 : 0 ;  // {0,0,1}
					idx |= (v101 < isoValue) ? 32 : 0 ;  // {1,0,1}
					idx |= (v111 < isoValue) ? 64 : 0 ;  // {1,1,1}
					idx |= (v011 < isoValue) ? 128 : 0 ; // {0,1,1}

					// Cube is entirely in/out of the surface
					if (!edges[idx]){
						v000 = v001;
						v100 = v101;
						v110 = v111;
						v010 = v011;
						continue;
					}

					// Find the vertices where the surface intersects the cube
					int edgeBit = 1;
					int vertexIndexList[12];
					#define MAP_INDEX(x,y,z) (( (((x+i)*res) + (y+j))*res)+(z+k)-1)

					#define ss stepSize
					for(int c = 0; c < 12; c++) {
						if(edges[idx] & edgeBit){
							switch(c){
								case 0:	
									vertexIndexList[c] = MAP_INDEX(0,0,0) * 3;//vIndex*3;
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x,y,z,1,1,1),
										vertex3<T>(x+ss,y,z,1,1,1),
										v000,v100,isoValue,ss,f,dx,dy,dz);
									break;
								
								case 2:	
									vertexIndexList[c] = MAP_INDEX(0,1,0) * 3; //(vIndex+zres)*3; 
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x,y+ss,z,1,1,1),
										vertex3<T>(x+ss,y+ss,z,1,1,1),
										v010,v110,isoValue,ss,f,dx,dy,dz);
									break;
								case 4: 
									vertexIndexList[c] = MAP_INDEX(0,0,1) * 3; //(vIndex+1)*3;
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x+ss,y,z+ss,1,1,1),
										vertex3<T>(x,y,z+ss,1,1,1),
										v101,v001,isoValue,ss,f,dx,dy,dz);
									break;
								case 6: 
									vertexIndexList[c] = MAP_INDEX(0,1,1) * 3; //(vIndex+yres+1)*3;  
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x+ss,y+ss,z+ss,1,1,1),
										vertex3<T>(x,y+ss,z+ss,1,1,1),
										v111,v011,isoValue,ss,f,dx,dy,dz);
									break;
								/////////////

								case 1:	
									vertexIndexList[c] = MAP_INDEX(1,0,0) * 3 + 1;//(vIndex+yres*zres)*3+1;
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x+ss,y,z,1,1,1),
										vertex3<T>(x+ss,y+ss,z,1,1,1),
										v100,v110,isoValue,ss,f,dx,dy,dz);
									break; 

								case 3: 
									vertexIndexList[c] = MAP_INDEX(0,0,0) * 3 + 1;//vIndex*3+1;
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x,y,z,1,1,1),
										vertex3<T>(x,y+ss,z,1,1,1),
										v000,v010,isoValue,ss,f,dx,dy,dz);
									break; // 0, 0, 0

								case 5: 
									vertexIndexList[c] = MAP_INDEX(1,0,1) * 3 + 1;//(vIndex+yres*zres+1)*3+1;  
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x+ss,y,z+ss,1,1,1),
										vertex3<T>(x+ss,y+ss,z+ss,1,1,1),
										v101,v111,isoValue,ss,f,dx,dy,dz);
									break; // 1, 0, 1
								case 7: 
									vertexIndexList[c] = MAP_INDEX(0,0,1) * 3 + 1;//(vIndex+1)*3+1; 
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x,y,z+ss,1,1,1),
										vertex3<T>(x,y+ss,z+ss,1,1,1),
										v001,v011,isoValue,ss,f,dx,dy,dz);
									break; //0,0,1

								//////////////////////////////

								case 8: 
									vertexIndexList[c] = MAP_INDEX(0,0,0) * 3 + 2;//vIndex*3+2;  
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x,y,z,1,1,1),
										vertex3<T>(x,y,z+ss,1,1,1),
										v000,v001,isoValue,ss,f,dx,dy,dz);
									break; //0,0,0
								case 9: 
									vertexIndexList[c] = MAP_INDEX(1,0,0) * 3 + 2;//(vIndex+yres*zres)*3+2; 
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x+ss,y,z,1,1,1),
										vertex3<T>(x+ss,y,z+ss,1,1,1),
										v100,v101,isoValue,ss,f,dx,dy,dz);
									break; //1,0,0
								case 10:
									vertexIndexList[c] = MAP_INDEX(1,1,0) * 3 + 2;//(vIndex+yres*zres+zres)*3+2;
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x+ss,y+ss,z,1,1,1),
										vertex3<T>(x+ss,y+ss,z+ss,1,1,1),
										v110,v111,isoValue,ss,f,dx,dy,dz);
									break; //1,1,0
								case 11:
									vertexIndexList[c] = MAP_INDEX(0,1,0) * 3 + 2;//(vIndex+yres)*3+2;
									cacheVertex<L,I,O>(&vMap,vertexIndexList[c], 
										vertex3<T>(x,y+ss,z,1,1,1),
										vertex3<T>(x,y+ss,z+ss,1,1,1),
										v010,v011,isoValue,ss,f,dx,dy,dz);
									break; //0,1,0
							}
						}
						edgeBit<<=1;
					}

					// Create the triangles
					for(int c = 0; triangles[idx][c] != -1;){

						#pragma omp critical (add_face2)
						{
							faces.push_back(tuple<int,int,int>(
								vertexIndexList[triangles[idx][c++]],
								vertexIndexList[triangles[idx][c++]],
								vertexIndexList[triangles[idx][c++]])); 
						}
					}

					// Setup for next cube
					v000 = v001;
					v100 = v101;
					v110 = v111;
					v010 = v011;
				}
			}
		}
		processFaces(&vMap, faces);
	}

	bool writeSurface(const std::string &out) const {
		return plyfile.writePly(out);
	}

	void applyTransform(const matrix4x4<T> &Tr) {
		plyfile.transformMesh(Tr);
	}

private:
	ply_writer<T> plyfile;
	std::vector<std::tuple<int,int,int>> faces;
	
	template<class L, class I, class O>
	inline void cacheVertex(
			std::map<int, sisl::vertex3<T> > *vMap, 
			const int &vIndex, 
			const sisl::vertex3<T> &v1,
			const sisl::vertex3<T> &v2,
			const T &iVal1,
			const T &iVal2,
			const T &isoValue,
			const T &ss,
			L *l, 
			L *dx, 
			L *dy, 
			L *dz) {
		#pragma omp critical(cache_vertex)
		{
			if(vMap->find(vIndex) == vMap->end()) {		
				T intr = (isoValue -iVal1)/(iVal2-iVal1);
				if(fabs(iVal2-iVal1) < 0.0001)
					intr = 0.5;
				sisl::vector3<T> vnew = v1.p*(1.-intr) + v2.p*(intr);
				sisl::vector3<T> na = l->grad_f(vnew);
				
				vMap->insert(
					{
						vIndex,
						sisl::vertex3<T>(vnew, na)
					});
			}
		}
	}
	void processFaces(
			std::map<int, sisl::vertex3<T> > *vMap,
			const std::vector<std::tuple<int,int,int>> &faces){
		 std::map<int, int> faceMap;

		for(auto itr = vMap->begin(); itr != vMap->end(); ++itr) {
			sisl::vertex3<T> v = itr->second;
			int vIndex = itr->first;

			faceMap[vIndex] = plyfile.addVertex(v);
		}

		for(auto itr = faces.begin(); itr != faces.end(); ++itr) {
			std::tuple<int, int, int> t = *itr;

			int i1 = faceMap[std::get<0>(t)]; 
			int i2 = faceMap[std::get<1>(t)];
			int i3 = faceMap[std::get<2>(t)];

			plyfile.addTriangle(i1,i2,i3);
		}
	}

	const int edges[256]={
		   0,  265,  515,  778, 1030, 1295, 1541, 1804, 
		2060, 2309, 2575, 2822, 3082, 3331, 3593, 3840, 
		 400,  153,  915,  666, 1430, 1183, 1941, 1692, 
		2460, 2197, 2975, 2710, 3482, 3219, 3993, 3728, 
		 560,  825,   51,  314, 1590, 1855, 1077, 1340, 
		2620, 2869, 2111, 2358, 3642, 3891, 3129, 3376, 
		 928,  681,  419,  170, 1958, 1711, 1445, 1196, 
		2988, 2725, 2479, 2214, 4010, 3747, 3497, 3232, 
		1120, 1385, 1635, 1898,  102,  367,  613,  876, 
		3180, 3429, 3695, 3942, 2154, 2403, 2665, 2912, 
		1520, 1273, 2035, 1786,  502,  255, 1013,  764, 
		3580, 3317, 4095, 3830, 2554, 2291, 3065, 2800, 
		1616, 1881, 1107, 1370,  598,  863,   85,  348, 
		3676, 3925, 3167, 3414, 2650, 2899, 2137, 2384, 
		1984, 1737, 1475, 1226,  966,  719,  453,  204, 
		4044, 3781, 3535, 3270, 3018, 2755, 2505, 2240, 
		2240, 2505, 2755, 3018, 3270, 3535, 3781, 4044, 
		 204,  453,  719,  966, 1226, 1475, 1737, 1984, 
		2384, 2137, 2899, 2650, 3414, 3167, 3925, 3676, 
		 348,   85,  863,  598, 1370, 1107, 1881, 1616, 
		2800, 3065, 2291, 2554, 3830, 4095, 3317, 3580, 
		 764, 1013,  255,  502, 1786, 2035, 1273, 1520, 
		2912, 2665, 2403, 2154, 3942, 3695, 3429, 3180, 
		 876,  613,  367,  102, 1898, 1635, 1385, 1120, 
		3232, 3497, 3747, 4010, 2214, 2479, 2725, 2988, 
		1196, 1445, 1711, 1958,  170,  419,  681,  928, 
		3376, 3129, 3891, 3642, 2358, 2111, 2869, 2620, 
		1340, 1077, 1855, 1590,  314,   51,  825,  560, 
		3728, 3993, 3219, 3482, 2710, 2975, 2197, 2460, 
		1692, 1941, 1183, 1430,  666,  915,  153,  400, 
		3840, 3593, 3331, 3082, 2822, 2575, 2309, 2060, 
		1804, 1541, 1295, 1030,  778,  515,  265,    0, 
	};
	const int triangles[256][16] = {
		{ -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   3,   8,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   0,   9,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   9,   1,   8,   1,   3,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,   1,  10,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   3,   8,   2,   1,  10,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,  10,   2,   9,   2,   0,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   9,  10,   8,  10,   2,   8,   2,   3,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   2,  11,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   8,   0,  11,   0,   2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   0,   9,   3,   2,  11,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   8,   9,  11,   9,   1,  11,   1,   2,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,  11,   3,  10,   3,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,  11,   8,  10,   8,   0,  10,   0,   1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,  10,  11,   9,  11,   3,   9,   3,   0,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   9,  10,   8,  10,  11,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   7,   4,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  7,   4,   0,   3,   7,   0,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  7,   4,   8,   1,   0,   9,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   3,   7,   9,   1,   7,   4,   9,   7,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,   1,  10,   8,   7,   4,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,   1,  10,   7,   4,   0,   3,   7,   0,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  7,   4,   8,   9,  10,   2,   9,   2,   0,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,  10,   4,   4,  10,   2,   4,   2,   7,   7,   2,   3,  -1,  -1,  -1, -1},
		{  2,  11,   3,   4,   8,   7,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,   0,   2,   7,   4,   2,  11,   7,   2,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   0,   9,   3,   2,  11,   8,   7,   4,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   4,   9,   1,   7,   4,   1,   2,   7,   2,  11,   7,  -1,  -1,  -1, -1},
		{  4,   8,   7,   3,   1,  10,   3,  10,  11,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   4,   0,   7,   4,   1,  10,   7,   1,  11,   7,  10,  -1,  -1,  -1, -1},
		{  9,  10,  11,   9,  11,   3,   9,   3,   0,   8,   7,   4,  -1,  -1,  -1, -1},
		{  9,  10,   4,   4,  10,   7,  11,   7,  10,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   4,   5,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   9,   4,   3,   8,   0,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   4,   5,   0,   5,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   1,   3,   4,   5,   3,   8,   4,   3,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   9,   4,   1,  10,   2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,   1,  10,   0,   3,   8,   9,   4,   5,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,   0,   4,   2,   4,   5,   2,   5,  10,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,   5,  10,   2,   4,   5,   2,   3,   4,   3,   8,   4,  -1,  -1,  -1, -1},
		{  3,   2,  11,   9,   4,   5,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   9,   4,   0,   2,  11,   0,  11,   8,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   2,  11,   4,   5,   1,   0,   4,   1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,  11,   8,   2,  11,   4,   5,   2,   4,   1,   2,   5,  -1,  -1,  -1, -1},
		{  5,   9,   4,  11,   3,   1,  10,  11,   1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,  11,   8,  10,   8,   0,  10,   0,   1,   9,   4,   5,  -1,  -1,  -1, -1},
		{ 10,  11,   5,   5,  11,   3,   5,   3,   4,   4,   3,   0,  -1,  -1,  -1, -1},
		{ 10,  11,   5,   5,  11,   4,   8,   4,  11,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   9,   8,   7,   5,   8,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   7,   5,   0,   3,   5,   9,   0,   5,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  7,   5,   1,   8,   7,   1,   0,   8,   1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   7,   5,   1,   3,   5,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,  10,   2,   8,   7,   5,   9,   8,   5,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   7,   5,   0,   3,   5,   9,   0,   5,  10,   2,   1,  -1,  -1,  -1, -1},
		{  8,   2,   0,  10,   2,   8,   7,  10,   8,   5,  10,   7,  -1,  -1,  -1, -1},
		{ 10,   7,   5,   2,   7,  10,   7,   2,   3,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,  11,   3,   7,   5,   9,   7,   9,   8,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   2,   9,   9,   2,  11,   9,  11,   5,   5,  11,   7,  -1,  -1,  -1, -1},
		{  7,   5,   1,   8,   7,   1,   0,   8,   1,   2,  11,   3,  -1,  -1,  -1, -1},
		{  2,   5,   1,  11,   5,   2,   5,  11,   7,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,  11,   3,  10,   3,   1,   5,   9,   8,   7,   5,   8,  -1,  -1,  -1, -1},
		{  0,   7,   9,   9,   7,   5,   7,   0,   1,   1,  10,  11,   1,  11,   7, -1},
		{  8,   5,   0,   7,   5,   8,   3,   0,   5,  10,  11,   3,   5,  10,   3, -1},
		{ 10,  11,   7,  10,   7,   5,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   5,   6,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   3,   8,  10,   5,   6,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   1,   0,  10,   5,   6,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   5,   6,   3,   8,   9,   1,   3,   9,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   2,   1,   6,   1,   5,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   3,   8,   5,   6,   2,   1,   5,   2,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   2,   0,   6,   0,   9,   6,   9,   5,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   8,   9,   3,   8,   5,   6,   3,   5,   2,   3,   6,  -1,  -1,  -1, -1},
		{  2,  11,   3,   6,  10,   5,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,  10,   5,   8,   0,   2,  11,   8,   2,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   2,  11,   1,   0,   9,  10,   5,   6,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   8,   9,  11,   9,   1,  11,   1,   2,  10,   5,   6,  -1,  -1,  -1, -1},
		{  3,   1,   5,   3,   5,   6,   3,   6,  11,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   8,   6,   6,   8,   0,   6,   0,   5,   5,   0,   1,  -1,  -1,  -1, -1},
		{  3,   6,  11,   3,   5,   6,   3,   0,   5,   0,   9,   5,  -1,  -1,  -1, -1},
		{ 11,   8,   6,   6,   8,   5,   9,   5,   8,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   5,   6,   8,   7,   4,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   5,   6,   4,   0,   3,   4,   3,   7,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   5,   6,   8,   7,   4,   1,   0,   9,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   3,   7,   9,   1,   7,   4,   9,   7,   6,  10,   5,  -1,  -1,  -1, -1},
		{  8,   7,   4,   6,   2,   1,   6,   1,   5,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   2,   1,   6,   1,   5,   7,   4,   0,   3,   7,   0,  -1,  -1,  -1, -1},
		{  6,   2,   0,   6,   0,   9,   6,   9,   5,   4,   8,   7,  -1,  -1,  -1, -1},
		{  9,   3,   4,   4,   3,   7,   3,   9,   5,   5,   6,   2,   5,   2,   3, -1},
		{  8,   7,   4,  10,   5,   6,   3,   2,  11,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,   0,   2,   7,   4,   2,  11,   7,   2,  10,   5,   6,  -1,  -1,  -1, -1},
		{  1,   0,   9,   3,   2,  11,   8,   7,   4,  10,   5,   6,  -1,  -1,  -1, -1},
		{ 10,   5,   6,   1,   4,   9,   1,   7,   4,   1,   2,   7,   2,  11,   7, -1},
		{  3,   1,   5,   3,   5,   6,   3,   6,  11,   7,   4,   8,  -1,  -1,  -1, -1},
		{ 11,   1,   6,   6,   1,   5,   1,  11,   7,   7,   4,   0,   7,   0,   1, -1},
		{  8,   7,   4,   3,   6,  11,   3,   5,   6,   3,   0,   5,   0,   9,   5, -1},
		{  9,   5,  11,  11,   5,   6,  11,   7,   9,   7,   4,   9,  -1,  -1,  -1, -1},
		{  6,  10,   9,   4,   6,   9,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   8,   0,   4,   6,  10,   4,  10,   9,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   4,   6,   0,   6,  10,   0,  10,   1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,   6,   8,   8,   6,  10,   8,  10,   3,   3,  10,   1,  -1,  -1,  -1, -1},
		{  4,   6,   2,   4,   2,   1,   4,   1,   9,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,   6,   2,   4,   2,   1,   4,   1,   9,   0,   3,   8,  -1,  -1,  -1, -1},
		{  0,   4,   6,   0,   6,   2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,   6,   8,   8,   6,   3,   2,   3,   6,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,  11,   3,   9,   4,   6,  10,   9,   6,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   8,   0,  11,   0,   2,   6,  10,   9,   4,   6,   9,  -1,  -1,  -1, -1},
		{  0,   4,   6,   0,   6,  10,   0,  10,   1,   2,  11,   3,  -1,  -1,  -1, -1},
		{  1,   4,  10,  10,   4,   6,   4,   1,   2,   2,  11,   8,   2,   8,   4, -1},
		{ 11,   4,   6,   9,   4,  11,   3,   9,  11,   1,   9,   3,  -1,  -1,  -1, -1},
		{  1,  11,   0,   0,  11,   8,  11,   1,   9,   9,   4,   6,   9,   6,  11, -1},
		{  0,   4,   3,   3,   4,  11,   6,  11,   4,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   8,   4,  11,   4,   6,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   9,   8,   6,  10,   8,   7,   6,   8,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,  10,   9,   6,  10,   0,   3,   6,   0,   7,   6,   3,  -1,  -1,  -1, -1},
		{ 10,   1,   6,   1,   0,   6,   0,   7,   6,   0,   8,   7,  -1,  -1,  -1, -1},
		{  6,   3,   7,  10,   3,   6,   3,  10,   1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   2,   7,   7,   2,   1,   7,   1,   8,   8,   1,   9,  -1,  -1,  -1, -1},
		{  0,   7,   9,   3,   7,   0,   1,   9,   7,   6,   2,   1,   7,   6,   1, -1},
		{  6,   2,   7,   7,   2,   8,   0,   8,   2,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   2,   3,   6,   3,   7,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   9,   8,   6,  10,   8,   7,   6,   8,   3,   2,  11,  -1,  -1,  -1, -1},
		{  6,   9,   7,  10,   9,   6,  11,   7,   9,   0,   2,  11,   9,   0,  11, -1},
		{  3,   2,  11,  10,   1,   6,   1,   0,   6,   0,   7,   6,   0,   8,   7, -1},
		{  1,  11,   7,   2,  11,   1,   7,  10,   1,   7,   6,  10,  -1,  -1,  -1, -1},
		{ 11,   1,   6,   3,   1,  11,   7,   6,   1,   9,   8,   7,   1,   9,   7, -1},
		{  7,   6,  11,   9,   0,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,  11,   0,   0,  11,   3,   0,   8,   6,   8,   7,   6,  -1,  -1,  -1, -1},
		{  7,   6,  11,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   6,   7,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   7,  11,   0,   3,   8,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   0,   9,  11,   6,   7,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   7,  11,   8,   9,   1,   8,   1,   3,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   2,   1,  11,   6,   7,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   3,   8,   2,   1,  10,  11,   6,   7,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   6,   7,   0,   9,  10,   2,   0,  10,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   9,  10,   8,  10,   2,   8,   2,   3,  11,   6,   7,  -1,  -1,  -1, -1},
		{  3,   2,   6,   7,   3,   6,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   2,   6,   8,   0,   6,   7,   8,   6,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   0,   9,   6,   7,   3,   2,   6,   3,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   9,   7,   7,   9,   1,   7,   1,   6,   6,   1,   2,  -1,  -1,  -1, -1},
		{  7,   3,   1,   7,   1,  10,   7,  10,   6,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   7,   8,   0,   6,   7,   0,   1,   6,   1,  10,   6,  -1,  -1,  -1, -1},
		{  6,   9,  10,   0,   9,   6,   7,   0,   6,   3,   0,   7,  -1,  -1,  -1, -1},
		{  8,   9,   7,   7,   9,   6,  10,   6,   9,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,   8,  11,   6,   4,  11,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   4,   0,  11,   6,   0,   3,  11,   0,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   0,   9,   8,  11,   6,   8,   6,   4,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   6,   4,  11,   6,   9,   1,  11,   9,   3,  11,   1,  -1,  -1,  -1, -1},
		{ 10,   2,   1,   4,   8,  11,   6,   4,  11,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   4,   0,  11,   6,   0,   3,  11,   0,   1,  10,   2,  -1,  -1,  -1, -1},
		{  9,  10,   2,   9,   2,   0,   4,   8,  11,   6,   4,  11,  -1,  -1,  -1, -1},
		{ 11,   4,   3,   6,   4,  11,   2,   3,   4,   9,  10,   2,   4,   9,   2, -1},
		{  2,   6,   4,   3,   2,   4,   8,   3,   4,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   2,   6,   4,   0,   6,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,   6,   4,   3,   2,   4,   8,   3,   4,   9,   1,   0,  -1,  -1,  -1, -1},
		{  9,   6,   4,   1,   6,   9,   6,   1,   2,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   4,  10,  10,   4,   8,  10,   8,   1,   1,   8,   3,  -1,  -1,  -1, -1},
		{  1,   4,   0,  10,   4,   1,   4,  10,   6,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   6,   8,   8,   6,   4,   6,   3,   0,   0,   9,  10,   0,  10,   6, -1},
		{  9,  10,   6,   9,   6,   4,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   6,   7,   9,   4,   5,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   4,   5,  11,   6,   7,   0,   3,   8,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   6,   7,   5,   1,   0,   5,   0,   4,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   1,   3,   4,   5,   3,   8,   4,   3,  11,   6,   7,  -1,  -1,  -1, -1},
		{ 11,   6,   7,   9,   4,   5,   2,   1,  10,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   6,   7,   9,   4,   5,   2,   1,  10,   0,   3,   8,  -1,  -1,  -1, -1},
		{  2,   0,   4,   2,   4,   5,   2,   5,  10,   6,   7,  11,  -1,  -1,  -1, -1},
		{ 11,   6,   7,   2,   5,  10,   2,   4,   5,   2,   3,   4,   3,   8,   4, -1},
		{  9,   4,   5,   7,   3,   2,   7,   2,   6,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   2,   6,   8,   0,   6,   7,   8,   6,   5,   9,   4,  -1,  -1,  -1, -1},
		{  0,   4,   5,   0,   5,   1,   3,   2,   6,   7,   3,   6,  -1,  -1,  -1, -1},
		{  4,   1,   8,   5,   1,   4,   7,   8,   1,   2,   6,   7,   1,   2,   7, -1},
		{  7,   3,   1,   7,   1,  10,   7,  10,   6,   5,   9,   4,  -1,  -1,  -1, -1},
		{  9,   4,   5,   0,   7,   8,   0,   6,   7,   0,   1,   6,   1,  10,   6, -1},
		{  6,   3,  10,   7,   3,   6,   5,  10,   3,   0,   4,   5,   3,   0,   5, -1},
		{ 10,   6,   8,   8,   6,   7,   8,   4,  10,   4,   5,  10,  -1,  -1,  -1, -1},
		{  9,   8,  11,   5,   9,  11,   6,   5,  11,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   0,   5,   0,   3,   5,   3,   6,   5,   3,  11,   6,  -1,  -1,  -1, -1},
		{  8,  11,   0,   0,  11,   6,   0,   6,   1,   1,   6,   5,  -1,  -1,  -1, -1},
		{ 11,   1,   3,   6,   1,  11,   1,   6,   5,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   8,  11,   5,   9,  11,   6,   5,  11,   2,   1,  10,  -1,  -1,  -1, -1},
		{  2,   1,  10,   9,   0,   5,   0,   3,   5,   3,   6,   5,   3,  11,   6, -1},
		{ 10,   0,   5,   2,   0,  10,   6,   5,   0,   8,  11,   6,   0,   8,   6, -1},
		{  5,   2,   3,  10,   2,   5,   3,   6,   5,   3,  11,   6,  -1,  -1,  -1, -1},
		{  3,   9,   8,   5,   9,   3,   2,   5,   3,   6,   5,   2,  -1,  -1,  -1, -1},
		{  5,   2,   6,   9,   2,   5,   2,   9,   0,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   5,   0,   0,   5,   1,   5,   8,   3,   3,   2,   6,   3,   6,   5, -1},
		{  1,   2,   6,   5,   1,   6,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   3,  10,  10,   3,   1,   3,   6,   5,   5,   9,   8,   5,   8,   3, -1},
		{  6,   9,   0,   5,   9,   6,   0,  10,   6,   0,   1,  10,  -1,  -1,  -1, -1},
		{  6,   5,  10,   8,   3,   0,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  6,   5,  10,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  7,  11,  10,   5,   7,  10,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,   3,   8,  11,  10,   5,  11,   5,   7,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   1,   0,   7,  11,  10,   5,   7,  10,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   9,   1,   8,   1,   3,   7,  11,  10,   5,   7,  10,  -1,  -1,  -1, -1},
		{  1,   5,   7,   1,   7,  11,   1,  11,   2,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   5,   7,   1,   7,  11,   1,  11,   2,   3,   8,   0,  -1,  -1,  -1, -1},
		{  2,   0,  11,  11,   0,   9,  11,   9,   7,   7,   9,   5,  -1,  -1,  -1, -1},
		{  2,   5,  11,  11,   5,   7,   5,   2,   3,   3,   8,   9,   3,   9,   5, -1},
		{  5,   7,   3,   5,   3,   2,   5,   2,  10,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   0,   2,   8,   0,  10,   5,   8,  10,   7,   8,   5,  -1,  -1,  -1, -1},
		{  5,   7,   3,   5,   3,   2,   5,   2,  10,   1,   0,   9,  -1,  -1,  -1, -1},
		{ 10,   7,   2,   5,   7,  10,   1,   2,   7,   8,   9,   1,   7,   8,   1, -1},
		{  3,   1,   5,   3,   5,   7,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   5,   0,   0,   5,   8,   7,   8,   5,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   7,   9,   9,   7,   0,   3,   0,   7,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   9,   5,   8,   5,   7,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,  11,  10,   4,   8,  10,   5,   4,  10,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,  10,   3,   3,  10,   5,   3,   5,   0,   0,   5,   4,  -1,  -1,  -1, -1},
		{  8,  11,  10,   4,   8,  10,   5,   4,  10,   1,   0,   9,  -1,  -1,  -1, -1},
		{  9,   3,   4,   1,   3,   9,   5,   4,   3,  11,  10,   5,   3,  11,   5, -1},
		{  2,   8,  11,   4,   8,   2,   1,   4,   2,   5,   4,   1,  -1,  -1,  -1, -1},
		{  2,   5,  11,   1,   5,   2,   3,  11,   5,   4,   0,   3,   5,   4,   3, -1},
		{  5,   2,   9,   9,   2,   0,   2,   5,   4,   4,   8,  11,   4,  11,   2, -1},
		{  5,   4,   9,  11,   2,   3,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   3,   4,   3,   2,   4,   2,   5,   4,   2,  10,   5,  -1,  -1,  -1, -1},
		{ 10,   0,   2,   5,   0,  10,   0,   5,   4,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   0,   9,   8,   3,   4,   3,   2,   4,   2,   5,   4,   2,  10,   5, -1},
		{  4,   1,   2,   9,   1,   4,   2,   5,   4,   2,  10,   5,  -1,  -1,  -1, -1},
		{  3,   1,   8,   8,   1,   4,   5,   4,   1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   4,   0,   1,   5,   0,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  5,   4,   3,   3,   4,   8,   3,   0,   5,   0,   9,   5,  -1,  -1,  -1, -1},
		{  5,   4,   9,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,  10,   9,   7,  11,   9,   4,   7,   9,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,  10,   9,   7,  11,   9,   4,   7,   9,   0,   3,   8,  -1,  -1,  -1, -1},
		{  1,  11,  10,   7,  11,   1,   0,   7,   1,   4,   7,   0,  -1,  -1,  -1, -1},
		{  4,   1,   8,   8,   1,   3,   1,   4,   7,   7,  11,  10,   7,  10,   1, -1},
		{ 11,   2,   7,   2,   1,   7,   1,   4,   7,   1,   9,   4,  -1,  -1,  -1, -1},
		{  0,   3,   8,  11,   2,   7,   2,   1,   7,   1,   4,   7,   1,   9,   4, -1},
		{  2,   0,  11,  11,   0,   7,   4,   7,   0,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,   7,   2,   2,   7,  11,   2,   3,   4,   3,   8,   4,  -1,  -1,  -1, -1},
		{ 10,   9,   2,   2,   9,   4,   2,   4,   3,   3,   4,   7,  -1,  -1,  -1, -1},
		{  8,   2,   7,   0,   2,   8,   4,   7,   2,  10,   9,   4,   2,  10,   4, -1},
		{ 10,   7,   2,   2,   7,   3,   7,  10,   1,   1,   0,   4,   1,   4,   7, -1},
		{  4,   7,   8,  10,   1,   2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  7,   3,   4,   4,   3,   9,   1,   9,   3,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   9,   7,   7,   9,   4,   7,   8,   1,   8,   0,   1,  -1,  -1,  -1, -1},
		{  0,   4,   7,   0,   7,   3,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  4,   7,   8,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,  11,  10,   9,   8,  10,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  0,  10,   9,   3,  10,   0,  10,   3,  11,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,  11,  10,   0,  11,   1,  11,   0,   8,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,  11,  10,   1,   3,  10,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,   8,  11,   1,   8,   2,   8,   1,   9,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   3,  11,   0,   3,   9,  11,   1,   9,  11,   2,   1,  -1,  -1,  -1, -1},
		{  0,   8,  11,   2,   0,  11,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 11,   2,   3,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  3,   9,   8,   2,   9,   3,   9,   2,  10,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  2,  10,   9,   0,   2,   9,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ 10,   0,   8,   1,   0,  10,   8,   2,  10,   8,   3,   2,  -1,  -1,  -1, -1},
		{ 10,   1,   2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  1,   9,   8,   3,   1,   8,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  9,   0,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{  8,   3,   0,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
		{ -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1},
	};
};
};
};

#endif // _ISOSURFACE_H_