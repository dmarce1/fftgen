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
#define MWEIGHT (1.5)
constexpr int NPAR = 4;
#define MAXFFT 128
#define DFFT 1

struct fft_type {
	int type;
	int N1;
	int N2;
	int N3;
	int nops;
};

void indent();
void deindent();
FILE* get_fp();
int get_ntabs();
void raders_fft_real_inv(int N, int o);
void include(const char* name);
void set_file(std::string str);
std::complex<double> twiddle(int k, int N);
std::string to_str(std::complex<double> z);
int square_factor(int N);
bool is_prime(int n);
int greatest_prime_factor(int N);
int mod_pow(int a, int b, int m);
int mod_inv(int a, int m);
int generator(int N);
const std::vector<std::complex<double>> twiddles(int N);
int index(int o, int i, int j, int N);
int index_real(int o, int i, int j, int N);
int index_real_inv(int o, int i, int j, int N);
void fft_bitreverse_real_inv(int N);
std::vector<int> fft_bitreverse_indices_real_inv(int N);
std::vector<int> fft_bitr_real_inv(int N, int o, std::vector<int> indices, bool first = false);
void fftw(std::vector<std::complex<double>>& x);
void fftw_real(std::vector<std::complex<double>>& xout, const std::vector<double>& xin);
std::unordered_map<int, int> prime_fac(int N);
std::vector<int> fft_bitr_real(int N, int o, std::vector<int> indices, bool first = false);
void fft_radix(int r, int N, int o);
std::vector<int> fft_radix_bitr(int r, int N, int o, std::vector<int> indices);
std::vector<int> fft_bitreverse_indices(int N);
std::vector<int> fft_radix_bitr_real(int r, int N, int o, std::vector<int> indices);
std::vector<int> fft_radix_bitr_real(int r, int N, int o, std::vector<int> I);
std::vector<int> fft_bitreverse_indices_real(int N);
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
void raders_fft(int r, int N, int o);
int raders_fft_opcnt(int r, int N);
int print_z_opcnt(int zi, int k, int r, int N);
int fft_radix_real_opcnt(int r, int N, int o);
fft_type best_radix_real(int N, int o, bool first = false);
int fft_radix_real_opcnt(int r, int N, int o);
void fft_radix_real(int r, int N, int o, bool = false);
int fft_radix_real_opcnt(int r, int N, int o);
fft_type best_radix(int N, int o, bool first = false);
std::vector<int> gt2_fft_bitr(int N1, int N2, int o, std::vector<int> I);
std::vector<int> fft_bitr(int N, int o, std::vector<int> indices, bool first = false);
void fft_real(int N, int o, bool = false);
int fft_real_opcnt(int N, int o);
void fft_real_inv(int N, int o, bool = false);
void raders_fft_real(int N, int o);
int raders_fft_real_opcnt(int N, int o);

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
