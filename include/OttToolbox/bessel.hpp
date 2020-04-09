/*
 * bessel.hpp:
 *
 * Generates the result of fractional order Bessel functions for use in the toolbox. It is capable
 * of calculating both Spherical Bessel and Ricatti--Bessel type functions and their associated
 * derivatives.
 *
 *
 *
 *
 *  Created on: 01/02/2012
 *      Author: Alex Stilgoe
 */

#include <limits>
#include <complex>
#include <vector>
#include <cmath>

#ifndef BESSEL_HPP_
#define BESSEL_HPP_

namespace ott {
	
	//We want the template to unwrap the type inside the complex so we can use hankel and riccati hankel functions with complex args.
	template <typename U>
	struct remove_complex
	{
		using type = U;
	};
	
	template <typename T>
	struct remove_complex<std::complex<T>>
	{
		using type = T;
	};
	
	template <typename W>
	using RemoveComplex = typename remove_complex<W>::type;
	
	template <class T>
	class bessel{
		using U = RemoveComplex<T>;
	private:
		//This is NOT what we should do as it is not a good C++11 pattern.
		std::vector<T> jn;
		std::vector<T> yn;
		
		T x;
		T invx;
		U absx;
		
	public:
		
		bessel(int nmax);
		bessel(int nmax, T x);
		
		~bessel(){}
		
		void change_x( T x );
		
		// prototype spherical bessel functions of the first and second kind
		const std::vector<T> j() const { return jn; }
		const std::vector<T> y() const { return yn; }
		
		//The functions below assume that RVO is performed at compile time.
		
		// prototype spherical bessel functions of the third kind
		std::vector<std::complex<U>> h1();
		std::vector<std::complex<U>> h2();
		
		// prototype of the first and second kind derivatives.
		std::vector<T> dj();
		std::vector<T> dy();
		
		// prototype of the third kind derivatives.
		std::vector<std::complex<U>> dh1();
		std::vector<std::complex<U>> dh2();
		
		// prototype ricatti bessel functions of the first and second kind
		std::vector<T> S();
		std::vector<T> C();
		
		// prototype ricatti bessel functions of the third kind
		std::vector<std::complex<U>> xi();
		std::vector<std::complex<U>> zeta();
		
		// prototype derivatives of ricatti bessel functions of the first and second kind
		std::vector<T> dS();
		std::vector<T> dC();
		
		// prototype derivatives of ricatti bessel functions of the third kind
		std::vector<std::complex<U>> dxi();
		std::vector<std::complex<U>> dzeta();
		
	};
	
	template <class T>
	bessel<T>::bessel(int n, T x_in) {
		
		this->jn.assign(n+1,(T)0);
		this->yn.assign(n+1,(T)0);
		
		change_x(x_in);
		
	}
	
	template <class T>
	bessel<T>::bessel(int n) {
		
		this->jn.assign(n+1,(T)0);
		this->yn.assign(n+1,(T)0);
		
	}
	
	
	
	
	/// Spherical Bessel Functions!
	
	
	
	/** Calculate the spherical bessel function of the third kind.
	 *
	 */
	template <class T>
	std::vector<std::complex<typename remove_complex<T>::type>> bessel<T>::h1() {
		using U = typename remove_complex<T>::type;
		std::vector<std::complex<U>> temp(yn.size(),std::complex<U>(0.,0.));
		
		if (absx == 0) {
			temp[0] = std::complex<U>(1.,-std::numeric_limits<U>::infinity());
			for ( int i = 1; i<yn.size(); i++) {
				temp[i] = std::complex<U>(0,-std::numeric_limits<U>::infinity());
			}
			
		} else {
			for ( int i = 0; i<yn.size(); i++) {
				temp[i] = jn[i]+ std::complex<U>(0.,1.)*yn[i];
			}
		}
		
		return temp;
	}
	
	/** Calculate the spherical bessel function of the third kind.
	 *
	 */
	template <class T>
	std::vector<std::complex<typename remove_complex<T>::type>> bessel<T>::h2() {
		using U = typename remove_complex<T>::type;
		std::vector<std::complex<U>> temp(yn.size(),std::complex<U>(0.,0.));
		
		if (absx == 0) {
			temp[0] = std::complex<U>(1.,std::numeric_limits<U>::infinity());
			for ( int i = 1; i<yn.size(); i++) {
				temp[i] = std::complex<U>(0,std::numeric_limits<U>::infinity());
			}
			
		} else {
			for ( int i = 0; i<yn.size(); i++) {
				temp[i] = jn[i]- std::complex<U>(0.,1.)*yn[i];
			}
		}
		return temp;
	}
	
	/** Calculate the derivative of the spherical bessel function of the first kind.
	 *
	 */
	template <class T>
	std::vector<T> bessel<T>::dj() {
		std::vector<T> temp(yn.size(),0.);
		
		if (absx!=0) {
			temp[0] = -yn[0]-jn[0]*this->invx;
		}
		
		if (yn.size()>1) {
			temp[1]=1./3.;
			if (absx!=0) {
				temp[1]=jn[0]-2.*jn[1]*this->invx;
			}
		}
		
		if (yn.size()>2 & absx!=0) {
			for ( int i = 2; i< yn.size(); i++ ) {
				temp[i]=jn[i-1]-(T)(i+1)*jn[i]*this->invx;
			}
		}
		return temp;
		
	}
	
	/** Calculate the derivative of the spherical bessel function of the second kind.
	 *
	 */
	template <class T>
	std::vector<T> bessel<T>::dy() {
		std::vector<T> temp(yn.size(),0.);
		
		temp[0] = jn[0]-yn[0]*this->invx;
		
		if (yn.size()>1) {
			if (absx == 0) {
				for ( int i = 1; i< yn.size(); i++ ) {
					temp[i]=std::numeric_limits<U>::infinity();
				}
			} else {
				for ( int i = 1; i< yn.size(); i++ ) {
					temp[i]=yn[i-1]-(T)(i+1)*yn[i]*this->invx;
				}
			}
		}
		return temp;
		
	}
	
	/** Calculate the derivative of the spherical bessel function of the third kind.
	 *
	 */
	template <class T>
	std::vector<std::complex<typename remove_complex<T>::type>> bessel<T>::dh1() {
		using U = typename remove_complex<T>::type;
		std::vector<std::complex<U>> temp(yn.size(),std::complex<U>(0.,0.));
		
		std::vector<T> djn = this->dj();
		std::vector<T> dyn = this->dy();
		
		if (absx == 0) {
			for ( int i = 0; i<yn.size(); i++) {
				temp[i] = std::complex<U>(0,std::numeric_limits<U>::infinity());
			}
			if (yn.size()>1)
				temp[1]=std::complex<U>(1./3.,std::numeric_limits<U>::infinity());
			
		} else {
			for ( int i = 0; i<yn.size(); i++) {
				temp[i] = djn[i]+ std::complex<U>(0.,1.)*dyn[i];
			}
		}
		
		return temp;
		
	}
	
	/** Calculate the derivative of the spherical bessel function of the third kind.
	 *
	 */
	template <class T>
	std::vector<std::complex<typename remove_complex<T>::type>> bessel<T>::dh2() {
		using U = typename remove_complex<T>::type;
		std::vector<std::complex<U>> temp(yn.size(),std::complex<U>(0.,0.));
		
		std::vector<T> djn = this->dj();
		std::vector<T> dyn = this->dy();
		
		if (absx == 0) {
			for ( int i = 0; i<yn.size(); i++) {
				temp[i] = std::complex<U>(0,-std::numeric_limits<U>::infinity());
			}
			if (yn.size()>1)
				temp[1]=std::complex<U>(1./3.,std::numeric_limits<U>::infinity());
			
		} else {
			for ( int i = 0; i<yn.size(); i++) {
				temp[i] = djn[i]- std::complex<U>(0.,1.)*dyn[i];
			}
		}
		
		return temp;
		
	}
	
	
	
	
	/// Ricatti--Bessel Functions!
	
	
	
	
	/** Calculate the ricatti bessel function of the first kind.
	 *
	 */
	template <class T>
	std::vector<T> bessel<T>::S() {
		std::vector<T> temp(jn);
		
		if ( absx != 0 ) {
			for (auto &i : temp) {
				i*=this->x;
			}
		}
		
		return temp;
		
	}
	
	/** Calculate the ricatti bessel function of the second kind.
	 *
	 */
	template <class T>
	std::vector<T> bessel<T>::C() {
		std::vector<T> temp(yn);
		
		if ( absx != 0 ) {
			for (auto &i : temp) {
				i*=this->x;
			}
		} else {
			temp[0]=-1.;
		}
		
		return temp;
		
	}
	
	/** Calculate the ricatti bessel function of the third kind.
	 *
	 */
	template <class T>
	std::vector<std::complex<typename remove_complex<T>::type>> bessel<T>::xi() {
		using U = typename remove_complex<T>::type;
		std::vector<std::complex<U>> temp(yn.size(),std::complex<U>(0.,0.));
		
		if (absx == 0) {
			temp[0] = std::complex<U>(0.,-1.);
			for ( int i = 1; i<yn.size(); i++) {
				temp[i] = std::complex<U>(0.,-std::numeric_limits<U>::infinity());
			}
			
		} else {
			for ( int i = 0; i<yn.size(); i++) {
				temp[i] = (jn[i]+ std::complex<U>(0.,1.)*yn[i])*this->x;
			}
		}
		
		return temp;
	}
	
	/** Calculate the ricatti bessel function of the third kind.
	 *
	 */
	template <class T>
	std::vector<std::complex<typename remove_complex<T>::type>> bessel<T>::zeta() {
		using U = typename remove_complex<T>::type;
		std::vector<std::complex<U>> temp(yn.size(),std::complex<U>(0.,0.));
		
		if (absx == 0) {
			temp[0] = std::complex<U>(0.,-1.);
			for ( int i = 1; i<yn.size(); i++) {
				temp[i] = std::complex<U>(0.,std::numeric_limits<U>::infinity());
			}
			
		} else {
			for ( int i = 0; i<yn.size(); i++) {
				temp[i] = (jn[i]- std::complex<U>(0.,1.)*yn[i])*this->x;
			}
		}
		
		return temp;
	}
	
	/** Calculate the derivative of the ricatti bessel function of the first kind.
	 *
	 */
	template <class T>
	std::vector<T> bessel<T>::dS() {
		std::vector<T> temp(yn.size(),0.);
		
		if (absx!=0) {
			temp[0] = -yn[0]*this->x;
			if (yn.size()>1 & absx!=0) {
				for ( int i = 1; i< yn.size(); i++ ) {
					temp[i]=jn[i-1]*this->x-(T)(i)*jn[i];
				}
			}
		} else {
			temp[0] = 1.;
		}
		
		
		return temp;
		
	}
	
	/** Calculate the derivative of the ricatti bessel function of the second kind.
	 *
	 */
	template <class T>
	std::vector<T> bessel<T>::dC() {
		std::vector<T> temp(yn.size(),0.);
		
		temp[0] = jn[0]*this->x;
		
		if (yn.size()>1) {
			if (absx == 0) {
				for ( size_t i = 1; i< yn.size(); i++ ) {
					temp[i]=std::numeric_limits<U>::infinity();
				}
			} else {
				for ( size_t i = 1; i< yn.size(); i++ ) {
					temp[i]=yn[i-1]*this->x-(T)(i)*yn[i];
				}
			}
		}
		return temp;
		
	}
	
	/** Calculate the derivative of the ricatti bessel function of the third kind.
	 *
	 */
	template <class T>
	std::vector<std::complex<typename remove_complex<T>::type>> bessel<T>::dxi() {
		using U = typename remove_complex<T>::type;
		std::vector<std::complex<U>> temp(yn.size(),std::complex<U>(0.,0.));
		
		std::vector<T> djn = this->dS();
		std::vector<T> dyn = this->dC();
		
		if (absx == 0) {
			
			for ( unsigned i = 0; i<yn.size(); i++) {
				temp[i] = std::complex<U>(0.,std::numeric_limits<U>::infinity());
			}
			temp[0] = std::complex<U>(1.,0.);
		} else {
			for ( size_t i = 0; i<yn.size(); i++) {
				temp[i] = djn[i]+ std::complex<U>(0.,1.)*dyn[i];
			}
		}
		
		return temp;
		
	}
	
	/** Calculate the derivative of the ricatti bessel function of the third kind.
	 *
	 */
	template <class T>
	std::vector<std::complex<typename remove_complex<T>::type>> bessel<T>::dzeta() {
		using U = typename remove_complex<T>::type;
		std::vector<std::complex<U>> temp(yn.size(),std::complex<U>(0.,0.));
		
		std::vector<T> djn = this->dS();
		std::vector<T> dyn = this->dC();
		
		if (absx == 0) {
			for ( size_t i = 0; i<yn.size(); i++) {
				temp[i] = std::complex<U>(0.,-std::numeric_limits<U>::infinity());
			}
			temp[0] = std::complex<U>(1.,0.);
			
		} else {
			for ( size_t i = 0; i<yn.size(); i++) {
				temp[i] = djn[i]- std::complex<U>(0.,1.)*dyn[i];
			}
		}
		return temp;
		
	}
	
	
	
	/** change_x is responsible for the entire calculation of j_n and h_n.
	 *
	 * All derivative functions use the results calculated here.
	 *
	 */
	template <class T>
	void bessel<T>::change_x(T x_in) {
		
		this->x = x_in;
		this->absx = std::abs(x_in);
		this->invx = (T)(1)/x_in;
		int maxN = yn.size();
		
		//Calculate spherical bessel y functions! Start with orders 0 and 1. We want to catch the x=0 special case.
		if (absx==0){
			
			yn[0] =-std::numeric_limits<U>::infinity();
			if (maxN>0) {
				yn[1] = -std::numeric_limits<U>::infinity();
			}
		} else {
			yn[0]=-std::cos(x_in)*this->invx;
			
			if (maxN>0) {
				yn[1] = (yn[0]-std::sin(x_in))*this->invx;
			}
			
		}
		
		if (yn.size()>2) {
			for ( int i=1; i< (int) yn.size()-1; i++) {
				if ( maxN>i ) {
					//if any operations are performed on inf we get nans. These functions tend to inf so we set to ceiling.
					if ( std::abs(yn[i])!=std::abs(yn[i]) || std::abs(yn[i])==std::numeric_limits<U>::infinity()  ) {
						maxN = i;
						yn[i+1]=-std::numeric_limits<U>::infinity();
					}
					else {
						yn[i+1] = (T)( 2 * i + 1 ) * yn[i] * this->invx - yn[i-1];
					}
				}
				else {
					yn[i+1]=-std::numeric_limits<U>::infinity();
				}
				
			}
		}
		
		//Calculate spherical bessel j functions using continued fractions method! NOTE: NEED TO FIND REFERENCE FOR DUE CREDIT!
		
		int N;
		T an;
		T c1;
		T c2;
		T Q;
		
		maxN--;
		
		//Define spherical bessel j functions for order 0 and 1
		if ( absx == 0 ) {
			jn[0] = (T)1.;
		}
		else {
			jn[0] = (T)std::sin(x_in)*this->invx;
		}
		
		if (jn.size()>1) {
			if ( absx != 0 ) {
				jn[1] = (jn[0]-std::cos(x_in))*this->invx;
			}
		}
		
		if (jn.size()>2) {
			if ( absx != 0 ) {
				
				// we need to work out continued fractions to find the highest order...
				
				N = 2*( int)this->absx+10;
				
				c2 = (T)1;
				
				an=(T)(2*maxN+1)*this->invx;
				
				c1 = an;
				Q = c1/c2;
				
				an=(T)(2*maxN+3)*this->invx;
				
				c1=an-(T)1.0/c1;
				c2=an;
				
				Q=Q*c1/c2;
				
				// note this loop gets horrible when the numbers are really big... it is not worth doing for r > 1e5.
				for ( int i=2;i<N;i++) {
					an=(T)(2*(maxN+i)+1)*this->invx;
					
					c1 = an-(T)1.0/c1;
					c2 = an-(T)1.0/c2;
					
					Q = Q*c1/c2;
				}
				
				
				jn[maxN-1] = (T)1.0 / (yn[maxN-1]/Q-yn[maxN]) * this->invx * this->invx; // This uses the Wronskian to find the bessel function!
				
				jn[maxN] = jn[maxN-1] / Q; // Q = j_{n-1}/j_n
				
				for ( int i=maxN-1;i>2;i-- ) {
					jn[i-1] =  (T)(2*i+1) * this->invx * jn[i] - jn[i+1];
				}
				
			}
			
		}
		
	}
}
#endif
