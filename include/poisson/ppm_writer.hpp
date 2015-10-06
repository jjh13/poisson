#ifndef _PPM_WRITER_H_
#define _PPM_WRITER_H_

#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <cfloat>

#include <sisl/sisl.hpp>
#include <sisl/array.hpp>
#include <sisl/utility/scattered.hpp>

template <class I = float, class O = float>
class ppm_writer
{
private:
	sisl::array2<sisl::vector3<unsigned char>> *m_parrImageData;
	sisl::array2<sisl::vector3<O>> *m_parrImageData2;

	int m_iHeight; 
	int m_iWidth;
	
	sisl::vector2<O> calc_minmax(){
		using namespace sisl;
		using namespace std; 

		O min = FLT_MAX;
		O max = -FLT_MAX;

		for(int i = 0; i < m_iWidth; i++)
			for(int j = 0; j < m_iHeight; j++) {
				vector3<O> data = this->at(i,j);
				if(data.i < min) min = data.i;
				if(data.j < min) min = data.j;
				if(data.k < min) min = data.k;
				if(data.i > max) max = data.i;
				if(data.j > max) max = data.j;
				if(data.k > max) max = data.k;
			}
		printf("%f %f\n", min, max );
		return vector2<O>(min, max);
	}

public:
	ppm_writer(const int &width, const int &height) : m_iHeight(height), m_iWidth(width) {
		using namespace sisl;
		m_parrImageData = new array2<vector3<unsigned char> >(width, height);
		m_parrImageData2 = new array2<vector3<O> >(width, height);
	}

	~ppm_writer() { delete m_parrImageData; };

	sisl::vector3<O> &at(const int &x, const int &y) {
		return (*m_parrImageData2)(x,y);
	}


	// sample
	bool write(const std::string &out) {
		using namespace sisl;
		using namespace std; 

		vector2<O> min_max = this->calc_minmax();

	    ofstream fp(out.c_str(), ios::out | ios::binary);
	    if(!fp.good()) false;
	    fp << "P3" << endl;
	    fp << m_iWidth << " " << m_iHeight << endl;
    	fp << "255" << endl;

		for(int i = 0; i < m_iWidth; i++)
			for(int j = m_iHeight - 1; j >= 0; j--) {
				vector3<O> pix = (*m_parrImageData2)(i,j);

				pix.i = 255. * (pix.i - min_max.i)/(min_max.j - min_max.i);
				pix.j = 255. * (pix.j - min_max.i)/(min_max.j - min_max.i);
				pix.k = 255. * (pix.k - min_max.i)/(min_max.j - min_max.i);


	    		fp << (int)pix.i << " " << (int)pix.j << " " << (int)pix.k << endl;
			}
		fp.close();
		return true;
	}
	/* data */
};

#endif