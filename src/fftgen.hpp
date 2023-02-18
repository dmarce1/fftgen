/*
 * fftgen.hpp
 *
 *  Created on: Feb 4, 2023
 *      Author: dmarce1
 */

#ifndef FFTGEN_HPP_
#define FFTGEN_HPP_

#include <string>
#include <cstdio>
#include <complex>
#include <vector>
#include <unordered_map>

#define RADERS 0
#define RADIX 1
#define GOOD 2
#define COOLEY 3
#define MWEIGHT (1.5)
constexpr int NPAR = 4;
#define MAXFFT 128
#define DFFT 1
#define LOOP_N 100
#define RADER_LEN 24

struct fft_type {
	int type;
	int N1;
	int N2;
	int N3;
	int nops;
};

void print_complex_short_fft(int r, std::vector<std::string> in_, std::vector<std::string> out_);
void print_real_short_fft(int r, std::vector<std::string> in_, std::vector<std::string> out_);
void fft_real(int N);
void indent();
void deindent();
FILE* get_fp();
int get_ntabs();
void include(const char* name);
void set_file(std::string str);
std::complex<double> twiddle(int k, int N);
std::string to_str(std::complex<double> z);
int square_factor(int N);
bool is_prime(int n);
int greatest_prime_factor(int N);
long long mod_pow(long long a, long long b, long long m);
int mod_inv(int a, int m);
int generator(long long N);
const std::vector<std::complex<double>> twiddles(int N);
int index(int o, int i, int j, int N);
int index_real(int o, int i, int j, int N);
void fftw(std::vector<std::complex<double>>& x);
void fftw_real(std::vector<std::complex<double>>& xout, const std::vector<double>& xin);
std::unordered_map<int, int> prime_fac(int N);
void fft_radix(int r, int N, int o);
std::vector<int> fft_radix_bitr(int r, int N, int o, std::vector<int> indices);
std::vector<int> fft_bitreverse_indices(int N);
void gt2_fft(int N1, int N2, int o);
int gt2_fft_opcnt(int N1, int N2);
int gt3_fft_opcnt(int N1, int N2, int N3);
void gt3_fft(int N1, int N2, int N3, int o);
int fft_radix_opcnt(int r, int N);
void fft(int N, int o, bool first = false);
int fft_opcnt(int N, bool first = false);
int fft_real_inv_opcnt(int N, int o);
void fft_bitreverse(int N, std::vector<int> indices = std::vector<int>(), int o = 0);
void fft_bitreverse_real(int N, std::vector<int> indices = std::vector<int>(), int o = 0);
std::vector<int> raders_ginvq(int N);
const std::vector<int> raders_gq(int N);
const std::vector<std::complex<double>> raders_four_twiddle(int N);
void print_z(int zi, int xi, int k, int r, int N, int o);
void raders_fft(int N, int o, bool padded = false);
int raders_fft_opcnt( int N);
int print_z_opcnt(int zi, int k, int r, int N);
fft_type best_radix(int N, int o, bool first = false);
std::vector<int> gt2_fft_bitr(int N1, int N2, int o, std::vector<int> I);
std::vector<int> fft_bitr(int N, int o, std::vector<int> indices, bool first = false);

template<class ... Args>
void print(const char* fstr, Args ...args) {
	for (int n = 0; n < get_ntabs(); n++) {
		fprintf(get_fp(), "\t");
	}
	fprintf(get_fp(), fstr, std::forward<Args>(args)...);
}

inline void print(const char* str) {
	for (int n = 0; n < get_ntabs(); n++) {
		fprintf(get_fp(), "\t");
	}
	fprintf(get_fp(), "%s", str);
}

template<class ... Args>
void print_notab(const char* fstr, Args ...args) {
	fprintf(get_fp(), fstr, std::forward<Args>(args)...);
}

inline void print_notab(const char* str) {
	fprintf(get_fp(), "%s", str);
}

#endif /* FFTGEN_HPP_ */
