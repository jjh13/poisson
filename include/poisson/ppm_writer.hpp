#ifndef _PPM_WRITER_H_
#define _PPM_WRITER_H_

#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>

#include <sisl/sisl.hpp>
#include <sisl/array.hpp>

template <class I = float, class O = float>
class ppm_writer
{
private:
	sisl::array2<sisl::vector3<unsigned char>> *m_parrImageData;

	int m_iHeight; 
	int m_iWidth;

public:
	ppm_writer(const int &width, const int &height) : m_iHeight(height), m_iWidth(width) {
		using namespace sisl;
		m_parrImageData = new array2<vector3<unsigned char> >(width, height);
	}

	~ppm_writer() { delete m_parrImageData; };

	sisl::vector3<unsigned char> &at(const int &x, const int &y) {
		return (*m_parrImageData)(x,y);
	}

	void set(const int &x, const int &y, sisl::vector3<unsigned char> value) {
		(*m_parrImageData)(x,y) = value;
	}

	// sample
	bool write(const std::string &out) {
		using namespace sisl;
		using namespace std; 

	    ofstream fp(out.c_str(), ios::out | ios::binary);
	    if(!fp.good()) false;
	    fp << "P3" << endl;
	    fp << m_iWidth << " " << m_iHeight << endl;
    	fp << "255" << endl;

		for(int i = 0; i < m_iWidth; i++)
			for(int j = m_iHeight - 1; j >= 0; j--) {
				vector3<unsigned char> pix = (*m_parrImageData)(i,j);
	    		fp << (int)pix.i << " " << (int)pix.j << " " << (int)pix.k << endl;
			}
		fp.close();
		return true;
	}

	/* data */
};

#endif