/*
 * legendre.hpp
 *
 * Class to implement associated legendre function and the associated theta/phi derivatives. It calculates only the theta elements as the phi part is identical for all of Y, dtY and dpY.
 *
 * The phi derivative is strictly imaginary but expressed here as a real. Actual calculations using spherical harmonics will have to account for this.
 *
 * Uncalculated values appear as not-a-number.
 *
 * Created on: 03/02/2012
 *     Author: Alex Stilgoe
 *
 * */

#include <limits>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#ifndef LEGENDRE_HPP_
#define LEGENDRE_HPP_

namespace ott {
	
#ifndef pi
#define pi 3.141592653589793238462643383279502884197169399375105820974944L
#endif
	
	template <class T>
	class legendre {
	private:
		
		std::vector<T> P_nm;
		T x; //x is theta!
		
		//This is to decide if the row or column recursion is a good idea. NYI: the col recursion.
		std::vector<int> n;
		std::vector<int> m;
		
		T dtRecursion(int n, int m, int nnp1);
		
	public:
		
		legendre(){this->n.assign(1,0),this->m.assign(1,0),this->P_nm.assign(1,std::numeric_limits<T>::quiet_NaN());}
		legendre(int nmax, T x);
		legendre(std::vector<int> n, std::vector<int> m, T x);
		
		std::vector<T> Y();
		std::vector<T> dtY();
		std::vector<T> dpY();
		
		T Y(int n, int m);
		T dtY(int n, int m);
		T dpY(int n, int m);
		
		std::vector<T> Y(std::vector<int> n, std::vector<int> m);
		std::vector<T> dtY(std::vector<int> n, std::vector<int> m);
		std::vector<T> dpY(std::vector<int> n, std::vector<int> m);
		
		std::vector<T> row(int n, T x);
		// NYI std::vector<T> col(T x);
		
		~legendre(){}
		
		void change_x(T x);
		
	};
	
	/** Constructor for a fully packed Legendre function with parameter x.
	 *
	 */
	template <class T>
	legendre<T>::legendre(int nmax, T x_in) {
		
		this->n.assign(nmax+1,0);
		for ( int i = 0; i< (int) nmax+1; i++ )
			n[i]=i;
		
		this->m.assign(nmax+1,0);
		for ( int i = 0; i< (int) nmax+1; i++ )
			m[i]=i;
		
		nmax++;
		
		this->P_nm.assign(nmax*(nmax+2),std::numeric_limits<T>::quiet_NaN());
		
		change_x(x_in);
		
	}
	
	/** Constructor for a Legendre function of vector n and m with parameter x. This currently won't save much effort because it uses the row recursion.
	 *
	 */
	template <class T>
	legendre<T>::legendre(std::vector<int> n_in, std::vector<int> m_in, T x_in) {
		
		std::vector<int>::iterator it;
		std::vector<int> ntemp (1,0);
		
		ntemp.insert(ntemp.end(),n_in.begin(),n_in.end());

		std::sort(ntemp.begin(),ntemp.end());
		it = std::unique(ntemp.begin(),ntemp.end());
		
		this->n.assign(ntemp.begin(),it);
		
		int nmax = (this->n.back()+1);
		
		std::vector<int> mtemp(m_in.begin(),m_in.end());
		
		std::sort(mtemp.begin(),mtemp.end());
		it = std::unique(mtemp.begin(),mtemp.end());
		
		this->m.assign(mtemp.begin(),it);
		
		
		this->P_nm.assign(nmax*(nmax+2),std::numeric_limits<T>::quiet_NaN());
		
		change_x(x_in);
		
	}
	
	/**	Re-compute the Y function at x.
	 *
	 */
	template <class T>
	void legendre<T>::change_x(T x_in) {
		
		this->x = x_in;
		
		//This is the part where we will compare the lengths of m and n and judge which set of recursions is quicker!
		//NYI
		
		
		// This is atrocious... need to get rid of all the copy spam!
		for ( int i = 1 ; i< (int) n.size(); i++ ) {
			
			std::vector<T> temp = row(this->n[i],x_in);
			
			for ( int j = 0; j< (int) temp.size(); j++ ) {
				
				this->P_nm[(n[i]-1)*(n[i]+1)+j] = temp[j];
				
			}
			
			if ( n[i+1] != n[i]+1) {
				std::vector<T> temp = row(this->n[i]+1,x_in);
				
				for ( int j = 0; j< (int) temp.size(); j++ ) {
					
					this->P_nm[(n[i])*(n[i]+2)+j] = temp[j];
					
				}
			}
			
		}
		
		
	}
	
	/**	Output the Y function.
	 *
	 */
	template <class T>
	std::vector<T> legendre<T>::Y() {
		
		std::vector<T> temp;
		
		if (P_nm.size()!=1) {
			
			temp.assign(this->P_nm.begin(),this->P_nm.end()-(2*this->n.back()+3));
			
		}
		else {
			
			temp.assign(1,std::numeric_limits<T>::quiet_NaN());
			
		}
		
		return temp;
	}
	
	/**	Output the Y function.
	 *
	 */
	template <class T>
	T legendre<T>::Y(int n_in, int m_in) {
		
		int ci = n_in*(n_in+1)+m_in;
		
		if ( n_in>n.back() ) {
			this->n.push_back(n_in);
			this->P_nm.resize((n_in+1)*(n_in+3),std::numeric_limits<T>::quiet_NaN());
			
			std::vector<T> temp = row(n_in,this->x);
			
			for ( int j = 0; j< (int) temp.size(); j++ ) {
				
				this->P_nm[(n_in-1)*(n_in+1)+j] = temp[j];
				
			}
			
		}
		else {
			if ( P_nm[ci-1]!=P_nm[ci-1] ) {
				std::vector<T> temp = row(n_in,this->x);
				
				for ( int j = 0; j< (int) temp.size(); j++ ) {
					
					this->P_nm[(n_in-1)*(n_in+1)+j] = temp[j];
					
				}
			}
		}
		
		
		return P_nm[ci-1];
		
	}
	
	/**	Output the Y function.
	 *
	 */
	template <class T>
	std::vector<T> legendre<T>::Y(std::vector<int> n_in, std::vector<int> m_in) {
		
		std::vector<T> temp_out (n_in.size(),(T)0);
		
		int ci;
		
		for ( int i=0; i< (int) n_in.size(); i++) {
			
			if ( n_in[i]>n.back() ) {
				
				this->n.push_back(n_in[i]);
				this->P_nm.resize((n_in[i]+1)*(n_in[i]+3),std::numeric_limits<T>::quiet_NaN());
			}
			
			ci=n_in[i]*(n_in[i]+1)+m_in[i];
			
			if (P_nm[ci-1]!=P_nm[ci-1]) {
				
				std::vector<T> temp = row(n_in[i],this->x);
				
				for ( int j = 0; j< (int) temp.size(); j++ ) {
					
					this->P_nm[(n_in[i]-1)*(n_in[i]+1)+j] = temp[j];
					
				}
				
			}
			
			temp_out[i]=P_nm[ci-1];
			
		}
		
		return temp_out;
		
	}
	
	/**	Calculate the theta derivative.
	 *
	 */
	template <class T>
	std::vector<T> legendre<T>::dtY() {
		std::vector<T> temp;
		
		if (P_nm.size()!=1) {
			temp.assign(P_nm.size()-(2*this->n.back()+3),std::numeric_limits<T>::quiet_NaN());
			
			int nt2np1;
			int nt;
			
			for ( int i = 1; i < (int) n.size(); i++) {
				
				nt = n[i];
				nt2np1= nt*(nt+1);
				
				for ( int mt = -nt; mt<=nt; mt++) {
					temp[nt2np1+mt-1] = dtRecursion(nt, mt, nt2np1);
				}
			}
			
		}
		else {
			temp.assign(1,std::numeric_limits<T>::quiet_NaN());
		}
		
		
		return temp;
		
	}
	
	/**	Calculate the theta derivative.
	 *
	 */
	template <class T>
	T legendre<T>::dtY(int n_in, int m_in) {
		
		int nt2np1=n_in*(n_in+1);
		int ci = nt2np1+m_in;
		
		if ( n_in>n.back() ) {
			this->n.push_back(n_in);
			this->P_nm.resize((n_in+1)*(n_in+3),std::numeric_limits<T>::quiet_NaN());
			
			std::vector<T> temp = row(n_in,this->x);
			
			for ( int j = 0; j< (int) temp.size(); j++ ) {
				
				this->P_nm[(n_in-1)*(n_in+1)+j] = temp[j];
				
			}
			
		}
		else {
			if ( P_nm[ci-1]!=P_nm[ci-1] ) {
				std::vector<T> temp = row(n_in,this->x);
				
				for ( int j = 0; j< (int) temp.size(); j++ ) {
					
					this->P_nm[(n_in-1)*(n_in+1)+j] = temp[j];
					
				}
			}
		}
		
		return dtRecursion(n_in, m_in, nt2np1);
		
	}
	
	/**	Calculate the theta derivative.
	 *
	 */
	template <class T>
	std::vector<T> legendre<T>::dtY(std::vector<int> n_in, std::vector<int> m_in) {
		
		std::vector<T> temp_out (n_in.size(),(T)0);
		
		int ci;
		
		for ( int i=0; i< (int) n_in.size(); i++) {
			
			if ( n_in[i]>n.back() ) {
				this->n.push_back(n_in[i]);
				this->P_nm.resize((n_in[i]+1)*(n_in[i]+3),std::numeric_limits<T>::quiet_NaN());
			}
			
			ci=n_in[i]*(n_in[i]+1)+m_in[i];
			
			if (P_nm[ci-1]!=P_nm[ci-1]) {
				
				std::vector<T> temp = row(n_in[i],this->x);
				
				for ( int j = 0; j< (int) temp.size(); j++ ) {
					
					this->P_nm[(n_in[i]-1)*(n_in[i]+1)+j] = temp[j];
					
				}
				
			}
			
			temp_out[i]=dtRecursion(n_in[i], m_in[i], ci-m_in[i]);
			
		}
		
		return temp_out;
		
	}
	
	
	/**	Calculate the phi derivative.
	 *
	 */
	template <class T>
	std::vector<T> legendre<T>::dpY() {
		std::vector<T> temp;
		if (P_nm.size()!=1) {
			temp.assign(P_nm.size()-(2*this->n.back()+3),std::numeric_limits<T>::quiet_NaN());
			
			int nt;
			int nt2np1;
			int ntnp1;
			
			T a_temp;
			
			for ( int i = 1; i < (int) n.size(); i++) {
				
				nt = n[i];
				nt2np1=(nt+1)*(nt+2);
				ntnp1=nt*(nt+1);
				
				a_temp = std::sqrt((T)(2*nt+1)/(T)(2*nt+3))*(T).5;
				
				for ( int mt = -nt; mt<=nt; mt++) {
					// normal
					temp[ntnp1+mt-1] = a_temp * (std::sqrt((T)(nt+mt+1)*(nt+mt+2)) * P_nm[nt2np1+mt] + std::sqrt((T)(nt-mt+1)*(nt-mt+2)) * P_nm[nt2np1+mt-2]);
				}
			}
			
		}
		else {
			temp.assign(n.size(),std::numeric_limits<T>::quiet_NaN());
		}
		
		return temp;
		
	}
	
	/**	Calculate the phi derivative.
	 *
	 */
	template <class T>
	T legendre<T>::dpY(int n_in, int m_in) {
		
		T a_temp = std::sqrt((T)(2*n_in+1)/(T)(2*n_in+3))*(T).5;
		int nt2np1=(n_in+1)*(n_in+2);
		int ci = (n_in+1)*(n_in+2)+m_in;
		
		if ( n_in>n.back() ) {
			this->n.push_back(n_in);
			this->P_nm.resize((n_in+1)*(n_in+3),std::numeric_limits<T>::quiet_NaN());
			
			std::vector<T> temp = row(n_in+1,this->x);
			
			for ( int j = 0; j< (int) temp.size(); j++ ) {
				
				this->P_nm[(n_in)*(n_in+2)+j] = temp[j];
				
			}
			
		}
		else {
			if ( P_nm[ci-1]!=P_nm[ci-1] ) {
				std::vector<T> temp = row(n_in,this->x);
				
				for ( int j = 0; j< (int) temp.size(); j++ ) {
					
					this->P_nm[(n_in)*(n_in+2)+j] = temp[j];
					
				}
			}
		}
		
		return a_temp * (std::sqrt((T)(n_in+m_in+1)*(n_in+m_in+2)) * P_nm[nt2np1+m_in] + std::sqrt((T)(n_in-m_in+1)*(n_in-m_in+2)) * P_nm[nt2np1+m_in-2]);
		
		
	}
	
	/**	Calculate the phi derivative.
	 *
	 */
	template <class T>
	std::vector<T> legendre<T>::dpY(std::vector<int> n_in, std::vector<int> m_in) {
		
		std::vector<T> temp_out (n_in.size(),(T)0);
		
		int ci;
		int nt2np1;
		T a_temp;
		
		for ( int i=0; i< (int) n_in.size(); i++) {
			
			if ( n_in[i]>n.back() ) {
				this->n.push_back(n_in[i]);
				this->P_nm.resize((n_in[i]+1)*(n_in[i]+3),std::numeric_limits<T>::quiet_NaN());
			}
			
			ci=(n_in[i]+1)*(n_in[i]+2)+m_in[i];
			nt2np1=(n_in[i]+1)*(n_in[i]+2);
			a_temp = std::sqrt((T)(2*n_in[i]+1)/(T)(2*n_in[i]+3))*(T).5;
			
			if (P_nm[ci-1]!=P_nm[ci-1]) {
				
				std::vector<T> temp = row(n_in[i]+1,this->x);
				
				for ( int j = 0; j< (int) temp.size(); j++ ) {
					
					this->P_nm[(n_in[i])*(n_in[i]+2)+j] = temp[j];
					
				}
				
			}
			
			temp_out[i]=a_temp * (std::sqrt((T)(n_in[i]+m_in[i]+1)*(n_in[i]+m_in[i]+2)) * P_nm[nt2np1+m_in[i]] + std::sqrt((T)(n_in[i]-m_in[i]+1)*(n_in[i]-m_in[i]+2)) * P_nm[nt2np1+m_in[i]-2]);
			
		}
		
		return temp_out;
		
	}
	
	
	/**	Row recursion implementation for Legendre functions.
	 *
	 */
	template <class T>
	std::vector<T> legendre<T>::row(int n_in, T x_in) {
		
		std::vector<T> temp (2*n_in+1,(T)0);
		
		if ( n_in==0 ) {
			
			temp[0] = (T)1;
			
		}
		else {
			
			T cx = std::cos(x_in);
			T sx = std::sin(x_in);
			
			T prod = (T)1;
			
			for ( int i=1; i<=n_in; i++) {
				
				prod *= (T)1 -  (T)0.5 / (T)i;
				
			}
			
			T Wnn =std::sqrt( (T)0.25*(T)(  2*n_in + 1  ) * prod / (T)pi );
			
			temp.back() = Wnn;
			temp[2*n_in-1] = std::sqrt( (T)(2 * n_in) ) * cx * Wnn;
			
			
			T sx2 = sx * sx;
			
			if ( n_in>1 ) {
				
				T a;
				T b;
				T I_n;
				
				for ( int i=n_in-2; i>=0; i-- ) {
					
					I_n = (T)1/(T)((n_in-i)*(n_in+i+1));
					
					a = std::sqrt((T)(4*(i+1)*(i+1))*I_n);
					b = std::sqrt((T)((n_in-i-1)*(n_in+i+2))*I_n);
					
					temp[n_in+i]=a*cx*temp[n_in + i + 1]-b*sx2*temp[n_in + i + 2];
					temp[n_in + i + 2] *= std::pow( sx, i+2 );
					
				}
				
			}
			
			temp[n_in + 1] *= sx;
			
		}
		
		T minus=(T)1;
		
		if (n_in % 2 != 0) {
			minus=(T)(-1);
		}
		
		for ( int i = 0; i<n_in; i++) {
			temp[i]=minus*temp[2*n_in-i];
			minus*=(T)(-1);
		}
		
		return temp;
		
	}
	
	
	/**	Recursion for theta derivatives of Legendre functions. This is implemented because it's needed for three cases.
	 *
	 */
	template <class T>
	T legendre<T>::dtRecursion(int nt, int mt, int nt2np1) {
		T out;
		// theta
		T b = std::sqrt((T)(nt-mt+1)*(nt+mt))*(T).5;
		T a = std::sqrt((T)(nt-mt)*(nt+mt+1))*(T).5;
		
		if ( ( abs(mt+1) > nt ) | ( abs(mt-1) > nt ) ) {
			// abnormal
			if ( mt < 0 ) {
				out = -a * P_nm[nt2np1+mt];
			}
			else {
				out = b * P_nm[nt2np1+mt-2];
			}
			
		}
		else {
			// normal
			out = -a * P_nm[nt2np1+mt] + b * P_nm[nt2np1+mt-2];
		}
		
		return out;
	}
	
}


#endif
