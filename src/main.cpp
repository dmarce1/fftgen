#include <cstdio>
#include <string>
#include <complex>
#include <cmath>
#include <unordered_map>
#include <set>
#include <memory>
#include <cassert>
#include <vector>
#include <fftw3.h>

#define MAXFFT 128
#define DFFT 1

FILE* fp = nullptr;

int tabs = 0;

void indent() {
	tabs++;
}

void deindent() {
	tabs--;
}

bool null = false;

void set_null() {
	null = true;
}

void clear_null() {
	null = false;
}

constexpr int NPAR = 4;

void set_file(std::string str) {
	if (fp) {
		fclose(fp);
	}
	if (str.empty()) {
		fp = nullptr;
	} else {
		str = "./gen/" + str;
		fp = fopen(str.c_str(), "wt");
	}
}

template<class ... Args>
void print(const char* fstr, Args ...args) {
	for (int n = 0; n < tabs; n++) {
		fprintf(fp, "\t");
	}
	fprintf(fp, fstr, std::forward<Args>(args)...);
}

void print(const char* str) {
	for (int n = 0; n < tabs; n++) {
		fprintf(fp, "\t");
	}
	fprintf(fp, "%s", str);
}

std::complex<double> twiddle(int k, int N) {
	return std::polar(1.0, -2.0 * M_PI * ((k + N) % N) / N);
}

std::string to_str(std::complex<double> z) {
	char* str;
	asprintf(&str, "std::complex<double>((%24.17e), (%24.17e))", z.real(), z.imag());
	std::string rc = str;
	free(str);
	return rc;
}

int square_factor(int N) {
	int N1 = sqrt(N);
	while (N % N1 != 0) {
		N1--;
	}
	return N1;
}

bool is_prime(int n) {
	bool v = true;
	if (n == 1) {
		v = false;
	} else {
		int kmax = sqrt(n);
		for (auto i = 2; i <= kmax; i++) {
			if (n % i == 0) {
				v = false;
				break;
			}
			if (i > kmax) {
				break;
			}
		}
	}
	return v;
}

int greatest_prime_factor(int N) {
	int v;
	for (int n = 2; n <= N; n++) {
		if (N % n == 0 && is_prime(n)) {
			v = n;
		}
	}
	return v;
}

static int mod_pow(int a, int b, int m) {
	int rc = 1;
	int apow = a;
	while (b) {
		if (b & 1) {
			rc = ((rc % m) * (apow % m)) % m;
		}
		b >>= 1;
		apow = ((apow % m) * (apow % m)) % m;
	}
	return rc;
}

static int mod_inv(int a, int m) {
	return mod_pow(a, m - 2, m);
}

static int generator(int N) {
	for (int g = 2; g < 1000; g++) {
		bool I[N];
		for (int i = 0; i < N; i++) {
			I[i] = false;
		}
		bool fail = false;
		for (int m = 0; m < N - 1; m++) {
			int n = mod_pow(g, m, N);
			assert(n >= 0);
			assert(n < N);
			if (!I[n]) {
				I[n] = true;
			} else {
				fail = true;
				break;
			}
		}
		if (!fail) {
			return g;
		}
	}
	assert(false);
}

std::vector<int> raders_ginvq(int N) {
	const int g = generator(N);
	std::vector<int> ginvq;
	for (int q = 0; q < N - 1; q++) {
		ginvq.push_back(mod_inv(mod_pow(g, q, N), N));
	}
	return ginvq;
}

const std::vector<std::complex<double>> twiddles(int N) {
	std::vector<std::complex<double>> tw(N);
	for (int k = 0; k < N; k++) {
		tw[k] = std::polar(1.0, -2 * M_PI * k / N);
	}
	return tw;
}

const std::vector<int> raders_gq(int N) {
	const int g = generator(N);
	std::vector<int> gq;
	for (int q = 0; q < N - 1; q++) {
		gq.push_back(mod_pow(g, q, N));
	}
	return gq;
}

constexpr int bluestein_gpf = 2;

int index(int o, int i, int j, int N) {
	return 2 * o + 2 * i + j;
}

int bluestein_size(int N) {
	int M = 2 * N - 1;
	while (greatest_prime_factor(M) > bluestein_gpf) {
		M++;
	}
	return M;
}

static inline std::vector<std::complex<double>> blustein_series(int N);
static inline std::vector<std::complex<double>> four_blustein_series(int N);

static inline std::vector<std::complex<double>> blustein_series(int N) {
	const int M = bluestein_size(N);
	std::vector<std::complex<double>> b(M, 0);
	b[0] = 1.0;
	for (int n = 1; n < N; n++) {
		b[M - n] = b[n] = std::polar(1.0, M_PI * n * n / N);
	}
	return b;
}

void fftw(std::vector<std::complex<double>>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, fftw_complex*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);
		plans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n][0] = x[n].real();
		i[n][1] = x[n].imag();
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n].real(o[n][0]);
		x[n].imag(o[n][1]);
	}

}

static inline std::vector<std::complex<double>> four_blustein_series(int N) {
	const int M = bluestein_size(N);
	auto vec = blustein_series(N);
	fftw(vec);
	for (int m = 0; m < M; m++) {
		vec[m] /= M;
	}
	return vec;
}

const std::vector<std::complex<double>> raders_four_twiddle(int N) {
	std::vector<std::complex<double>> b(N - 1);
	const auto tws = twiddles(N);
	const auto ginvq = raders_ginvq(N);
	for (int q = 0; q < N - 1; q++) {
		b[q] = tws[ginvq[q]];
	}
	fftw(b);
	return b;
}

void fft_radix2(int N, int o);
void fft_radix3(int N, int o);
void fft_radix4(int N, int o);
void fft_radix5(int N, int o);
void fft_radix6(int N, int o);
void fft_radix(int r, int N, int o);
std::vector<int> fft_radix_bitr(int r, int N, int o, std::vector<int> indices);
std::vector<int> fft_bitreverse_indices(int N);

int fft_radix_opcnt(int r, int N);
void fft(int N, int o);
int fft_opcnt(int N);

int bluestein_fft_opcnt(int N) {
	int cnt = 0;
	const int M = bluestein_size(N);
	if (M > MAXFFT) {
		return 999999999;
	}
	cnt += 21 * N;
	cnt += 2 * fft_opcnt(M);
	return cnt;
}

void fft_bitreverse(int N, std::vector<int> indices = std::vector<int>(), int o = 0);

/*void bluestein_fft(int N, int o) {
 //	printf("Bluestein %i\n", N);
 const int M = bluestein_size(N);
 const std::vector<std::complex<double>> fourb = four_blustein_series(N);
 const std::vector<std::complex<double>> b = blustein_series(N);
 print("{\n");
 indent();
 print("// Bluestein radix - %i / %i \n", N, M);
 print("std::array<std::complex<double>, %i> y;\n", M);
 for (int n = 0; n < N; n++) {
 print("y[%i] = std::conj(x[%i]) * %s;\n", n, o + n, to_str(b[n]).c_str());
 }
 for (int n = N; n < M; n++) {
 print("y[%i] = 0;\n", n);
 }
 print("{\n");
 indent();
 print("auto* x = y.data();\n");
 fft_bitreverse(M);
 fft(M, 0);
 deindent();
 print("}\n");
 for (int k = 0; k < M; k++) {
 print("y[%i] = std::conj(y[%i]) * %s;\n", k, k, to_str(fourb[k]).c_str());
 }
 print("{\n");
 indent();
 print("auto* x = y.data();\n");
 fft_bitreverse(M);
 fft(M, 0);
 deindent();
 print("}\n");
 for (int n = 0; n < N; n++) {
 print("x[%i] = y[%i] * %s;\n", o + n, n, to_str(std::conj(b[n])).c_str());
 }
 deindent();
 print("}\n");
 }*/

void raders_fft(int r, int N, int o) {
	const int N1 = r;
	const int N2 = N / r;
	print("{\n");
	indent();
	print("// Raders radix - %i = %i x %i\n", N, N1, N2);
	printf("// Raders radix - %i = %i x %i\n", N, N1, N2);
	print("std::array<double, %i> xro, xrk0;\n", N2);
	print("std::array<double, %i> xio, xik0;\n", N2);
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	if (N2 > 1) {
		std::vector<int> I1(N);
		auto I2 = fft_bitreverse_indices(N2);
		for (int n2 = 0; n2 < N2; n2++) {
			for (int n1 = 0; n1 < N1; n1++) {
				I1[N2 * n1 + n2] = N1 * I2[n2] + n1;
			}
		}
		fft_bitreverse(N, I1, o);
		for (int n1 = 0; n1 < N1; n1++) {
			print("{\n");
			indent();
			fft(N2, o + n1 * N2);
			deindent();
			print("}\n");
		}
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				I1[N1 * n2 + n1] = N2 * n1 + n2;
			}
		}
		fft_bitreverse(N, I1, o);
	}
	const auto b = raders_four_twiddle(N1);
	const auto gq = raders_gq(N1);
	const auto ginvq = raders_ginvq(N1);
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			const int nk = (n1 * k2) % N;
			if (nk == 0) {
			} else {
				const auto W = twiddle(nk, N);
				print("tmp0 = x[%i];\n", index(o, k2 * N1 + n1, 0, N));
				print("x[%i] = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", index(o, k2 * N1 + n1, 0, N), W.real(), index(o, k2 * N1 + n1, 0, N), -W.imag(), index(o, k2 * N1 + n1, 1, N));
				print("x[%i] = std::fma((%.17e), tmp0, (%.17e) * x[%i]);\n", index(o, k2 * N1 + n1, 1, N), W.imag(), W.real(), index(o, k2 * N1 + n1, 1, N));
			}
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		print("xro[%i] = x[%i];\n", k2, index(o, k2 * N1, 0, N));
		print("xrk0[%i] = x[%i] ", k2, index(o, k2 * N1, 0, N));
		for (int k1 = 1; k1 < N1; k1++) {
			fprintf(fp, " + x[%i]", index(o, N1 * k2 + k1, 0, N));
		}
		fprintf(fp, ";\n");
		print("xio[%i] = x[%i];\n", k2, index(o, k2 * N1, 1, N));
		print("xik0[%i] = x[%i] ", k2, index(o, k2 * N1, 1, N));
		for (int k1 = 1; k1 < N1; k1++) {
			fprintf(fp, " + x[%i]", index(o, N1 * k2 + k1, 1, N));
		}
		fprintf(fp, ";\n");
	}

	std::vector<int> I1(N), I3(N);
	auto I2 = fft_bitreverse_indices(N1 - 1);
	for (int i = 0; i < I2.size(); i++) {
		I3[I2[i]] = i;
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int q = 0; q < N1 - 1; q++) {
			I1[N1 * k2 + I3[q]] = N1 * k2 + gq[q];
		}
		I1[N1 * k2 + N1 - 1] = N1 * k2;
	}
	fft_bitreverse(N, I1, o);

	for (int k2 = 0; k2 < N2; k2++) {
		print("{\n");
		indent();
		fft(N1 - 1, o + k2 * N1);
		deindent();
		print("}\n");
	}
	for (int k2 = 0; k2 < N2; k2++) {
		print("x[%i] = -x[%i];\n", index(o, k2 * N1, 0, N), index(o, k2 * N1, 0, N));
		for (int q = 1; q < N1 - 1; q++) {
			print("tmp0 = x[%i];\n", index(o, k2 * N1 + q, 0, N));
			print("x[%i] = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", index(o, k2 * N1 + q, 0, N), b[q].real(), index(o, k2 * N1 + q, 0, N), -b[q].imag(), index(o, k2 * N1 + q, 1, N));
			print("x[%i] = std::fma((%.17e), tmp0, (%.17e) * x[%i]);\n", index(o, k2 * N1 + q, 1, N), -b[q].imag(), -b[q].real(), index(o, k2 * N1 + q, 1, N));
		}
	}

	for (int k2 = 0; k2 < N2; k2++) {
		print("{\n");
		indent();
		fft_bitreverse(N1 - 1, fft_bitreverse_indices(N1 - 1), o + k2 * N1);
		fft(N1 - 1, o + k2 * N1);
		deindent();
		print("}\n");
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const auto Nm1inv = 1.0 / (N1 - 1.0);
		for (int q = 0; q < N1 - 1; q++) {
			print("x[%i] *= (%24.17e);\n", index(o, N1 * k2 + q, 0, N), Nm1inv);
			print("x[%i] = -x[%i] * (%24.17e);\n", index(o, N1 * k2 + q, 1, N), index(o, N1 * k2 + q, 1, N), Nm1inv);
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int p = N1 - 2; p >= 0; p--) {
			print("x[%i] = xro[%i] + x[%i];\n", index(o, N1 * k2 + p + 1, 0, N), k2, index(o, N1 * k2 + p, 0, N));
			print("x[%i] = xio[%i] + x[%i];\n", index(o, N1 * k2 + p + 1, 1, N), k2, index(o, N1 * k2 + p, 1, N));
		}
		print("x[%i] = xrk0[%i];\n", index(o, k2 * N1, 0, N), k2);
		print("x[%i] = xik0[%i];\n", index(o, k2 * N1, 1, N), k2);
	}
	std::vector<int> indices;
	indices.resize(N);
	for (int k2 = 0; k2 < N2; k2++) {
		indices[k2] = N1 * k2;
		for (int k1 = 0; k1 < N1 - 1; k1++) {
			indices[N2 * ginvq[k1] + k2] = N1 * k2 + k1 + 1;
		}
	}
	fft_bitreverse(N, indices, o);
	deindent();
	print("}\n");
}

constexpr double mweight = 2;

int raders_fft_opcnt(int r, int N) {
	int cnt = 0;
	const int N1 = r;
	const int N2 = N / r;
	if (N2 > 1) {
		cnt += mweight * N;
		for (int n1 = 0; n1 < N1; n1++) {
			cnt += fft_opcnt(N2);
		}
	}

	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			const int nk = (n1 * k2) % N;
			if (nk != 0) {
				cnt += 2;
			}
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int k1 = 1; k1 < N1; k1++) {
			cnt += 2;
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		cnt += mweight * (N1 - 1);
		cnt += fft_opcnt(N1 - 1);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		cnt++;
		for (int q = 1; q < N1 - 1; q++) {
			cnt += 2;
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		cnt += mweight * (N1 - 1);
		cnt += fft_opcnt(N1 - 1);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const auto Nm1inv = 1.0 / (N1 - 1.0);
		for (int q = 0; q < N1 - 1; q++) {
			cnt += 3;
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int p = N1 - 2; p >= 0; p--) {
			cnt += 2;
		}
	}
	cnt += mweight * N;
	return cnt;
}

int print_z_opcnt(int zi, int k, int r, int N) {
	if (zi * k == 0) {
		return 0;
	} else if (zi * k == N / 8 && N % 8 == 0) {
		return 4;
	} else if (zi * k == N / 4 && N % 4 == 0) {
		return 1;
	} else if (zi * k == 3 * N / 8 && N % 8 == 0) {
		return 4;
	} else if (zi * k == N / 2 && N % 2 == 0) {
		return 2;
	} else if (zi * k == 5 * N / 8 && N % 8 == 0) {
		return 4;
	} else if (zi * k == 3 * N / 4 && N % 4 == 0) {
		return 1;
	} else if (zi * k == 7 * N / 8 && N % 8 == 0) {
		return 4;
	} else {
		return 4;
	}
}

void print_z(int zi, int k, int r, int N, int o) {
	const auto W = twiddle(zi * k, N);
	const int i = k + zi * N / r;
	if (zi * k == 0) {
		print("const auto zr%i = x[%i];\n", zi, index(o, i, 0, N));
		print("const auto zi%i = x[%i];\n", zi, index(o, i, 1, N));
	} else if (zi * k == N / 8 && N % 8 == 0) {
		print("const auto zr%i = M_SQRT1_2 * (x[%i] + x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
		print("const auto zi%i = -M_SQRT1_2 * (x[%i] - x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	} else if (zi * k == N / 4 && N % 4 == 0) {
		print("const auto zr%i = x[%i];\n", zi, index(o, i, 1, N));
		print("const auto zi%i = -x[%i];\n", zi, index(o, i, 0, N));
	} else if (zi * k == 3 * N / 8 && N % 8 == 0) {
		print("const auto zr%i = -M_SQRT1_2 * (x[%i] - x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
		print("const auto zi%i = -M_SQRT1_2 * (x[%i] + x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	} else if (zi * k == N / 2 && N % 2 == 0) {
		print("const auto zr%i = -x[%i];\n", zi, index(o, i, 0, N));
		print("const auto zi%i = -x[%i];\n", zi, index(o, i, 1, N));
	} else if (zi * k == 5 * N / 8 && N % 8 == 0) {
		print("const auto zr%i = -M_SQRT1_2 * (x[%i] + x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
		print("const auto zi%i = M_SQRT1_2 * (x[%i] - x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	} else if (zi * k == 3 * N / 4 && N % 4 == 0) {
		print("const auto zr%i = -x[%i];\n", zi, index(o, i, 1, N));
		print("const auto zi%i = x[%i];\n", zi, index(o, i, 0, N));
	} else if (zi * k == 7 * N / 8 && N % 8 == 0) {
		print("const auto zr%i = M_SQRT1_2 * (x[%i] - x[%i];\n", zi, index(o, i, 0, N), index(o, i, 1, N));
		print("const auto zi%i = M_SQRT1_2 * (x[%i] + x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	} else {
		print("const auto zr%i = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", zi, W.real(), index(o, i, 0, N), -W.imag(), index(o, i, 1, N));
		print("const auto zi%i = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", zi, W.real(), index(o, i, 1, N), W.imag(), index(o, i, 0, N));
	}
}

int best_radix(int N, int o) {
	int best_cnt = 999999999;
	int best_radix = -1;
	for (int r = 2; r <= N; r++) {
		if (N % r == 0) {
			int this_cnt;
			if (r <= 6 || is_prime(r)) {
				this_cnt = fft_radix_opcnt(r, N) + mweight * N;
				if (this_cnt < best_cnt) {
					best_cnt = this_cnt;
					best_radix = r;
				}
				if (r > 6) {
					int raders_cnt = raders_fft_opcnt(r, N);
					if (raders_cnt < best_cnt) {
						best_cnt = raders_cnt;
						best_radix = r;
					}
				}
			}
		}
	}
	/*	if (greatest_prime_factor(N) > bluestein_gpf) {
	 int bluestein_cnt = bluestein_fft_opcnt(N);
	 if (bluestein_cnt < best_cnt) {
	 return -1;
	 }
	 }*/
	return best_radix;
}

std::vector<int> fft_bitr(int N, int o, std::vector<int> indices);

int fft_opcnt(int N) {
	int best_cnt = 999999999;
	int best_radix = -1;
	for (int r = 2; r <= N; r++) {
		if (N % r == 0) {
			if (r <= 6 || is_prime(r)) {
				const int this_cnt = fft_radix_opcnt(r, N) + mweight * N;
				if (this_cnt < best_cnt) {
					best_cnt = this_cnt;
					best_radix = r;
				}
				if (r > 6) {
					int raders_cnt = raders_fft_opcnt(r, N);
					if (raders_cnt < best_cnt) {
						best_cnt = raders_cnt;
						best_radix = r;
					}
				}
			}
		}
	}
	return best_cnt;
}

int fft_radix_opcnt(int r, int N) {
	int cnt = 0;
	switch (r) {
	case 1:
		return 0;
	case 2:
		if (N < 2) {
			return 0;
		}
		cnt = 0;
		if (N > 2) {
			cnt += fft_opcnt(N / 2);
			cnt += fft_opcnt(N / 2);
		}
		for (int k = 0; k < N / 2; k++) {
			for (int i = 0; i < 2; i++) {
				cnt += print_z_opcnt(i, k, 2, N);
			}
			cnt += 4;
		}
		return cnt;
	case 3:
		if (N < 3) {
			return 0;
		}
		cnt = 0;
		if (N > 3) {
			cnt += fft_opcnt(N / 3);
			cnt += fft_opcnt(N / 3);
			cnt += fft_opcnt(N / 3);
		}
		for (int k = 0; k < N / 3; k++) {
			for (int i = 0; i < 3; i++) {
				cnt += print_z_opcnt(i, k, 3, N);
			}
			cnt += 14;
		}
		return cnt;
	case 4:
		if (N < 4) {
			return 0;
		}
		cnt = 0;
		if (N > 4) {
			cnt += fft_opcnt(N / 4);
			cnt += fft_opcnt(N / 4);
			cnt += fft_opcnt(N / 4);
			cnt += fft_opcnt(N / 4);
		}
		for (int k = 0; k < N / 4; k++) {
			for (int i = 0; i < 4; i++) {
				cnt += print_z_opcnt(i, k, 4, N);
			}
			cnt += 16;
		}
		return cnt;
	case 5:
		if (N < 5) {
			return 0;
		}
		cnt = 0;
		if (N > 5) {
			cnt += fft_opcnt(N / 5);
			cnt += fft_opcnt(N / 5);
			cnt += fft_opcnt(N / 5);
			cnt += fft_opcnt(N / 5);
			cnt += fft_opcnt(N / 5);
		}
		for (int k = 0; k < N / 5; k++) {
			for (int i = 0; i < 5; i++) {
				cnt += print_z_opcnt(i, k, 5, N);
			}
			cnt += 34;
		}
		return cnt;
	case 6:
		if (N < 6) {
			return 0;
		}
		cnt = 0;
		if (N > 6) {
			cnt += fft_opcnt(N / 6);
			cnt += fft_opcnt(N / 6);
			cnt += fft_opcnt(N / 6);
			cnt += fft_opcnt(N / 6);
			cnt += fft_opcnt(N / 6);
			cnt += fft_opcnt(N / 6);
		}
		for (int k = 0; k < N / 6; k++) {
			for (int i = 0; i < 6; i++) {
				cnt += print_z_opcnt(i, k, 6, N);
			}
			cnt += 40;
		}
		return cnt;
	default:
		if (N < r) {
			return 0;
		}
		cnt = 0;
		if (N > r) {
			for (int j = 0; j < r; j++) {
				cnt += fft_opcnt(N / r);
			}
		}
		for (int k = 0; k < N / r; k++) {
			for (int i = 0; i < r; i++) {
				cnt += print_z_opcnt(i, k, r, N);
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				cnt += 4;
				for (int j = 1; j <= (r - 1) / 2; j++) {
					cnt += 4;
				}
			}
			for (int i = 1; i < r; i++) {
				cnt += 2;
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				cnt += 4;
			}
		}
		return cnt;
	};
}

std::vector<int> fft_bitreverse_indices(int N) {
	std::vector<int> I;
	std::vector<int> indices;
	I.resize(0);
	for (int n = 0; n < N; n++) {
		I.push_back(n);
	}
	indices = fft_bitr(N, 0, I);
	return indices;
}

#include <list>

void fft_bitreverse(int N, std::vector<int> indices, int o) {
	if (indices.size() == 0) {
		indices = fft_bitreverse_indices(N);
	}
	std::vector<std::list<int>> strings;
	for (int l = 0; l < 2; l++) {
		std::vector<bool> touched(N, false);
		for (int k = 0; k < N; k++) {
			if (!touched[k]) {
				std::list<int> string;
				int current = k;
				int first = current;
				int next = indices[current];
				touched[current] = true;
				if (first != next) {
					string.push_back(-1);
					string.push_back(index(o, current, l, N));
					while (first != next) {
						string.push_back(index(o, next, l, N));
						current = next;
						touched[current] = true;
						next = indices[next];
					}
				}
				strings.push_back(string);
			}
		}
	}
	while (strings.size()) {
		for (int n = 0; n < NPAR; n++) {
			if (strings.size() > n) {
				while (!strings[n].size()) {
					if (strings.size() > NPAR) {
						strings[n] = strings.back();
						strings.pop_back();
					} else {
						break;
					}
				}
				if (strings[n].size()) {
					if (strings[n].front() == -1) {
						strings[n].pop_front();
						if (strings[n].size()) {
							print("tmp%i = x[%i];\n", n, strings[n].front());
						}
					} else {
						int to = strings[n].front();
						strings[n].pop_front();
						if (strings[n].size()) {
							print("x[%i] = x[%i];\n", to, strings[n].front());
						} else {
							print("x[%i] = tmp%i;\n", to, n);
						}
					}
					while (!strings[n].size()) {
						if (strings.size() > NPAR) {
							strings[n] = strings.back();
							strings.pop_back();
						} else {
							break;
						}
					}
				}
			}
		}
		bool done = true;
		for (int n = 0; n < NPAR; n++) {
			if (strings.size() > n) {
				if (strings[n].size()) {
					done = false;
					break;
				}
			}
		}
		if (done) {
			strings.resize(0);
		}
	}
}

void fft(int N, int o) {
	if (N == 1) {
		return;
	}
	int radix = best_radix(N, o);
//	if (radix > 0) {
	fft_radix(radix, N, o);
	//} else if (radix == -1) {
	//	bluestein_fft(N, o);
	//}
}

std::vector<int> fft_bitr(int N, int o, std::vector<int> indices) {
	if (N == 1) {
		return indices;
	}
	int radix = best_radix(N, o);
	if (radix > 0) {
		return fft_radix_bitr(radix, N, o, indices);
	} else {
		return indices;
	}
}

void fft_radix(int r, int N, int o) {
	std::vector<int> L;
	std::vector<int> J;
	std::vector<int> K;
	switch (r) {
	case 1:
		return;
	case 2:
		print("{\n");
		print("// radix - 2\n");
		indent();
		for (int n = 0; n < r; n++) {
			fft(N / r, o + n * N / r);
		}
		for (int k = 0; k < N / 2; k++) {
			print("{\n");
			indent();
			for (int i = 0; i < r; i++) {
				print_z(i, k, r, N, o);
			}
			print("x[%i] = zr0 + zr1;\n", index(o, k, 0, N));
			print("x[%i] = zi0 + zi1;\n", index(o, k, 1, N));
			print("x[%i] = zr0 - zr1;\n", index(o, k + N / 2, 0, N));
			print("x[%i] = zi0 - zi1;\n", index(o, k + N / 2, 1, N));
			deindent();
			print("}\n");
		}
		deindent();
		print("}\n");
		break;
	case 3:
		print("{\n");
		print("// radix - 3\n");
		indent();
		for (int n = 0; n < r; n++) {
			fft(N / r, o + n * N / r);
		}
		for (int k = 0; k < N / 3; k++) {
			print("{\n");
			indent();
			for (int i = 0; i < r; i++) {
				print_z(i, k, r, N, o);
			}
			print("const auto tr1 = zr1 + zr2;\n");
			print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
			print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
			print("const auto ti1 = zi1 + zi2;\n");
			print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
			print("const auto ti3 = (%24.17e) * (zi1 - zi2);\n", sqrt(3) * 0.5);
			print("x[%i] = zr0 + tr1;\n", index(o, k, 0, N));
			print("x[%i] = zi0 + ti1;\n", index(o, k, 1, N));
			print("x[%i] = tr2 + ti3;\n", index(o, k + N / 3, 0, N));
			print("x[%i] = ti2 - tr3;\n", index(o, k + N / 3, 1, N));
			print("x[%i] = tr2 - ti3;\n", index(o, k + 2 * N / 3, 0, N));
			print("x[%i] = ti2 + tr3;\n", index(o, k + 2 * N / 3, 1, N));
			deindent();
			print("}\n");
		}
		deindent();
		print("}\n");
		break;
	case 4:
		print("{\n");
		print("// radix - 4\n");
		indent();
		for (int n = 0; n < r; n++) {
			fft(N / r, o + n * N / r);
		}
		for (int k = 0; k < N / 4; k++) {
			print("{\n");
			indent();
			for (int i = 0; i < r; i++) {
				print_z(i, k, r, N, o);
			}
			print("const auto tr1 = zr0 + zr2;\n");
			print("const auto ti1 = zi0 + zi2;\n");
			print("const auto tr2 = zr1 + zr3;\n");
			print("const auto ti2 = zi1 + zi3;\n");
			print("const auto tr3 = zr0 - zr2;\n");
			print("const auto ti3 = zi0 - zi2;\n");
			print("const auto tr4 = zr1 - zr3;\n");
			print("const auto ti4 = zi1 - zi3;\n");
			print("x[%i] = tr1 + tr2;\n", index(o, k, 0, N));
			print("x[%i] = ti1 + ti2;\n", index(o, k, 1, N));
			print("x[%i] = tr3 + ti4;\n", index(o, k + N / 4, 0, N));
			print("x[%i] = ti3 - tr4;\n", index(o, k + N / 4, 1, N));
			print("x[%i] = tr1 - tr2;\n", index(o, k + N / 2, 0, N));
			print("x[%i] = ti1 - ti2;\n", index(o, k + N / 2, 1, N));
			print("x[%i] = tr3 - ti4;\n", index(o, k + 3 * N / 4, 0, N));
			print("x[%i] = ti3 + tr4;\n", index(o, k + 3 * N / 4, 1, N));
			deindent();
			print("}\n");
		}
		deindent();
		print("}\n");
		break;
	case 5:
		print("{\n");
		print("// radix - 5\n");
		indent();
		for (int n = 0; n < r; n++) {
			fft(N / r, o + n * N / r);
		}
		for (int k = 0; k < N / 5; k++) {
			print("{\n");
			indent();
			for (int i = 0; i < r; i++) {
				print_z(i, k, r, N, o);
			}
			print("const auto tr1 = zr1 + zr4;\n");
			print("const auto tr2 = zr2 + zr3;\n");
			print("const auto tr3 = zr1 - zr4;\n");
			print("const auto tr4 = zr2 - zr3;\n");
			print("const auto tr5 = tr1 + tr2;\n");
			print("const auto tr6 = (%24.17e) * (tr1 - tr2);\n", sqrt(5) * 0.25);
			print("const auto tr7 = std::fma(tr5, -0.25, zr0);\n");
			print("const auto tr8 = tr7 + tr6;\n");
			print("const auto tr9 = tr7 - tr6;\n");
			print("const auto tr10 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
			print("const auto tr11 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
			print("const auto ti1 = zi1 + zi4;\n");
			print("const auto ti2 = zi2 + zi3;\n");
			print("const auto ti3 = zi1 - zi4;\n");
			print("const auto ti4 = zi2 - zi3;\n");
			print("const auto ti5 = ti1 + ti2;\n");
			print("const auto ti6 = (%24.17e) * (ti1 - ti2);\n", sqrt(5) * 0.25);
			print("const auto ti7 = std::fma(ti5, -0.25, zi0);\n");
			print("const auto ti8 = ti7 + ti6;\n");
			print("const auto ti9 = ti7 - ti6;\n");
			print("const auto ti10 = std::fma((%24.17e), ti3, (%24.17e) * ti4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
			print("const auto ti11 = std::fma((%24.17e), ti3, (%24.17e) * ti4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
			print("x[%i] = zr0 + tr5;\n", index(o, k + 0 * N / 5, 0, N));
			print("x[%i] = zi0 + ti5;\n", index(o, k + 0 * N / 5, 1, N));
			print("x[%i] = tr8 + ti10;\n", index(o, k + 1 * N / 5, 0, N));
			print("x[%i] = ti8 - tr10;\n", index(o, k + 1 * N / 5, 1, N));
			print("x[%i] = tr9 + ti11;\n", index(o, k + 2 * N / 5, 0, N));
			print("x[%i] = ti9 - tr11;\n", index(o, k + 2 * N / 5, 1, N));
			print("x[%i] = tr9 - ti11;\n", index(o, k + 3 * N / 5, 0, N));
			print("x[%i] = ti9 + tr11;\n", index(o, k + 3 * N / 5, 1, N));
			print("x[%i] = tr8 - ti10;\n", index(o, k + 4 * N / 5, 0, N));
			print("x[%i] = ti8 + tr10;\n", index(o, k + 4 * N / 5, 1, N));
			deindent();
			print("}\n");
		}
		deindent();
		print("}\n");
		break;
	case 6:
		print("{\n");
		print("// radix - 6\n");
		indent();
		for (int n = 0; n < r; n++) {
			fft(N / r, o + n * N / r);
		}
		for (int k = 0; k < N / 6; k++) {
			print("{\n");
			indent();
			for (int i = 0; i < r; i++) {
				print_z(i, k, r, N, o);
			}
			print("const auto tr1 = zr2 + zr4;\n");
			print("const auto ti1 = zi2 + zi4;\n");
			print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
			print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
			print("const auto tr3 = (%24.17e) * (zr2 - zr4);\n", sin(M_PI / 3.0));
			print("const auto ti3 = (%24.17e) * (zi2 - zi4);\n", sin(M_PI / 3.0));
			print("const auto tr4 = zr5 + zr1;\n");
			print("const auto ti4 = zi5 + zi1;\n");
			print("const auto tr5 = std::fma(tr4, -0.5, zr3);\n");
			print("const auto ti5 = std::fma(ti4, -0.5, zi3);\n");
			print("const auto tr6 = (%24.17e) * (zr5 - zr1);\n", sin(M_PI / 3.0));
			print("const auto ti6 = (%24.17e) * (zi5 - zi1);\n", sin(M_PI / 3.0));
			print("const auto tr7 = zr0 + tr1;\n");
			print("const auto ti7 = zi0 + ti1;\n");
			print("const auto tr8 = tr2 + ti3;\n");
			print("const auto ti8 = ti2 - tr3;\n");
			print("const auto tr9 = tr2 - ti3;\n");
			print("const auto ti9 = ti2 + tr3;\n");
			print("const auto tr10 = zr3 + tr4;\n");
			print("const auto ti10 = zi3 + ti4;\n");
			print("const auto tr11 = tr5 + ti6;\n");
			print("const auto ti11 = ti5 - tr6;\n");
			print("const auto tr12 = tr5 - ti6;\n");
			print("const auto ti12 = ti5 + tr6;\n");
			print("x[%i] = tr7 + tr10;\n", index(o, k + 0 * N / 6, 0, N));
			print("x[%i] = tr8 - tr11;\n", index(o, k + 1 * N / 6, 0, N));
			print("x[%i] = tr9 + tr12;\n", index(o, k + 2 * N / 6, 0, N));
			print("x[%i] = tr7 - tr10;\n", index(o, k + 3 * N / 6, 0, N));
			print("x[%i] = tr8 + tr11;\n", index(o, k + 4 * N / 6, 0, N));
			print("x[%i] = tr9 - tr12;\n", index(o, k + 5 * N / 6, 0, N));
			print("x[%i] = ti7 + ti10;\n", index(o, k + 0 * N / 6, 1, N));
			print("x[%i] = ti8 - ti11;\n", index(o, k + 1 * N / 6, 1, N));
			print("x[%i] = ti9 + ti12;\n", index(o, k + 2 * N / 6, 1, N));
			print("x[%i] = ti7 - ti10;\n", index(o, k + 3 * N / 6, 1, N));
			print("x[%i] = ti8 + ti11;\n", index(o, k + 4 * N / 6, 1, N));
			print("x[%i] = ti9 - ti12;\n", index(o, k + 5 * N / 6, 1, N));
			deindent();
			print("}\n");
		}
		deindent();
		print("}\n");
		break;
	default:
		int raders_cnt = raders_fft_opcnt(r, N);
		int this_cnt = fft_radix_opcnt(r, N);
		if (raders_cnt < this_cnt) {
			raders_fft(r, N, o);
			return;
		} else {
			print("{\n");
			print("// radix - %i\n", r);
			indent();
			for (int n = 0; n < r; n++) {
				fft(N / r, o + n * N / r);
			}
			for (int k = 0; k < N / r; k++) {
				print("{\n");
				indent();
				for (int i = 0; i < r; i++) {
					print_z(i, k, r, N, o);
				}
				for (int j = 1; j <= (r - 1) / 2; j++) {
					print("const auto txp%i = zr%i + zr%i;\n", j, j, r - j);
				}
				for (int j = 1; j <= (r - 1) / 2; j++) {
					print("const auto txm%i = zr%i - zr%i;\n", j, j, r - j);
				}
				for (int j = 1; j <= (r - 1) / 2; j++) {
					print("const auto typ%i = zi%i + zi%i;\n", j, j, r - j);
				}
				for (int j = 1; j <= (r - 1) / 2; j++) {
					print("const auto tym%i = zi%i - zi%i;\n", j, j, r - j);
				}
				for (int i = 1; i <= (r - 1) / 2; i++) {
					print("auto ap%i = zr0;\n", i);
					print("auto bp%i = zi0;\n", i);
					for (int j = 1; j <= (r - 1) / 2; j++) {
						print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
						print("bp%i = std::fma(typ%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
						if (j == 1) {
							print("double am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
							print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
						} else {
							print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
							print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
						}
					}
				}
				print("x[%i] = ", index(o, k, 0, N));
				for (int i = 0; i < r; i++) {
					fprintf(fp, " + zr%i", i);
				}
				fprintf(fp, ";\n");
				print("x[%i] = ", index(o, k, 1, N));
				for (int i = 0; i < r; i++) {
					fprintf(fp, " + zi%i", i);
				}
				fprintf(fp, ";\n");
				for (int i = 1; i <= (r - 1) / 2; i++) {
					print("x[%i] = ap%i + am%i;\n", index(o, k + i * N / r, 0, N), i, i);
					print("x[%i] = bp%i - bm%i;\n", index(o, k + i * N / r, 1, N), i, i);
					print("x[%i] = ap%i - am%i;\n", index(o, k + (r - i) * N / r, 0, N), i, i);
					print("x[%i] = bp%i + bm%i;\n", index(o, k + (r - i) * N / r, 1, N), i, i);
				}
				deindent();
				print("}\n");
			}
			deindent();
			print("}\n");
		}
	};
}

std::vector<int> fft_radix_bitr(int r, int N, int o, std::vector<int> I) {
	std::vector<int> L;
	std::vector<int> J;
	std::vector<int> K;
	switch (r) {
	case 1:
		return I;
	case 2:
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
		break;
	case 3:
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
		break;
	case 4:
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
		break;
	case 5:
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
		break;
	case 6:
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
		break;
	default:
		int raders_cnt = raders_fft_opcnt(r, N);
		int this_cnt = fft_radix_opcnt(r, N);
		if (raders_cnt < this_cnt) {
			return I;
		} else {
			for (int n = 0; n < r; n++) {
				std::vector<int> J;
				for (int k = 0; k < N / r; k++) {
					J.push_back(I[n + r * k]);
				}
				auto K = fft_bitr(N / r, o + n * N / r, J);
				L.insert(L.end(), K.begin(), K.end());
			}

		}
	};
	return L;
}

void print_fft(int N) {
	std::string fname = "fft." + std::to_string(N) + ".cpp";
	set_file(fname);
	print("#include \"fft.hpp\"\n");
	print("\nvoid fft_%i(std::complex<double>* x0) {\n", N);
	indent();
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	print("double* x = reinterpret_cast<double*>(x0);\n");
	std::vector<int> I;
	for (int n = 0; n < N; n++) {
		I.push_back(n);
	}
	fft_bitreverse(N);
	fft(N, 0);
	deindent();
	print("}\n\n");
}

int main(int argc, char **argv) {

	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print_fft(n);
	}

	set_file("fft.hpp");
	print("#include <complex>\n");
	print("#include <cmath>\n");
	print("#include <chrono>\n");
	print("#include <unordered_map>\n");
	print("#include <vector>\n\n");
	print("#include <fftw3.h>\n\n");
	print("\nvoid FFT(std::complex<double>* x, int N);\n");
	print("\nvoid FFT_inv(std::complex<double>* x, int N);\n");

	std::string header4 = "inline std::complex<double> Ix(const std::complex<double>& a) {\n"
			"\tstd::complex<double> b;\n"
			"\tb.real(-a.imag());\n"
			"\tb.imag(a.real());\n"
			"\treturn b;\n"
			"}\n"
			"";

	std::string header5 = "inline std::complex<double> No8x(const std::complex<double>& a) {\n"
			"\tstd::complex<double> c;\n"
			"\tstatic constexpr double c0 = sqrt(2) * 0.5;\n"
			"\tc.real(c0 * (a.real() + a.imag()));\n"
			"\tc.imag(c0 * (a.imag() - a.real()));\n"
			"\treturn c;\n"
			"}\n"
			"";
	std::string header41 = "inline std::complex<double> nIx(const std::complex<double>& a) {\n"
			"\tstd::complex<double> b;\n"
			"\tb.real(a.imag());\n"
			"\tb.imag(-a.real());\n"
			"\treturn b;\n"
			"}\n"
			"";

	std::string header51 = "inline std::complex<double> nNo8x(const std::complex<double>& a) {\n"
			"\tstd::complex<double> c;\n"
			"\tstatic constexpr double c0 = -sqrt(2) * 0.5;\n"
			"\tc.real(c0 * (a.real() + a.imag()));\n"
			"\tc.imag(c0 * (a.imag() - a.real()));\n"
			"\treturn c;\n"
			"}\n"
			"";

	std::string header42 = "inline std::complex<double> No8xIx(const std::complex<double>& a) {\n"
			"\tstd::complex<double> c;\n"
			"\tstatic constexpr double c0 = sqrt(2) * 0.5;\n"
			"\tc.imag(c0 * (a.real() + a.imag()));\n"
			"\tc.real(c0 * (a.real() - a.imag()));\n"
			"\treturn c;\n"
			"}\n"
			"";

	std::string header52 = "inline std::complex<double> nNo8xIx(const std::complex<double>& a) {\n"
			"\tstd::complex<double> c;\n"
			"\tstatic constexpr double c0 = -sqrt(2) * 0.5;\n"
			"\tc.imag(c0 * (a.real() + a.imag()));\n"
			"\tc.real(c0 * (a.real() - a.imag()));\n"
			"\treturn c;\n"
			"}\n"
			"";

	print("%s\n", header4.c_str());
	print("%s\n", header5.c_str());
	print("%s\n", header41.c_str());
	print("%s\n", header51.c_str());
	print("%s\n", header42.c_str());
	print("%s\n", header52.c_str());

	set_file("fft.cpp");
	print("#include \"fft.hpp\"\n\n");

	print("\n"
			"#include <chrono>\n"
			"\n"
			"class timer {\n"
			"\tstd::chrono::time_point<std::chrono::high_resolution_clock> start_time;\n"
			"\tdouble time;\n"
			"public:\n"
			"\tinline timer() {\n"
			"\t\ttime = 0.0;\n"
			"\t}\n"
			"\tinline void stop() {\n"
			"\t\tstd::chrono::time_point<std::chrono::high_resolution_clock> stop_time = std::chrono::high_resolution_clock::now();\n"
			"\t\tstd::chrono::duration<double> dur = stop_time - start_time;\n"
			"\t\ttime += dur.count();\n"
			"\t}\n"
			"\tinline void start() {\n"
			"\t\tstart_time = std::chrono::high_resolution_clock::now();\n"
			"\t}\n"
			"\tinline void reset() {\n"
			"\t\ttime = 0.0;\n"
			"\t}\n"
			"\tinline double read() {\n"
			"\t\treturn time;\n"
			"\t}\n"
			"};\n"
			"\n"
			"");

	print("\n");
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print("void fft_%i(std::complex<double>*);\n", n);
	}
	print("\n");

	print("typedef void (*func_type)(std::complex<double>*);\n");
	print("\n");
	print("const func_type fptr[] = {nullptr, nullptr, ");
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		if (!null) {
			fprintf(fp, "&fft_%i", n);
			if (n != MAXFFT) {
				fprintf(fp, ", ");
			}
		}
	}
	if (!null) {
		fprintf(fp, "};\n\n");
	}
	print("void FFT(std::complex<double>* x, int N) {\n");
	indent();
	print("(*(fptr[N]))(x);\n");
	deindent();
	print("}\n\n");
	print("void FFT_inv(std::complex<double>* x, int N) {\n");
	indent();
	print("const double Ninv = 1.0 / N;\n");
	print("for( int n = 0; n < N; n++) {\n");
	indent();
	print("x[n] = std::conj(x[n]);\n");
	deindent();
	print("}\n");
	print("(*(fptr[N]))(x);\n");
	print("for( int n = 0; n < N; n++) {\n");
	indent();
	print("x[n] = std::conj(x[n]) * Ninv;\n");
	deindent();
	print("}\n");
	deindent();
	print("}\n\n");

	std::string header2 = "\n"
			"double rand1() {\n"
			"\treturn (rand() + 0.5) / RAND_MAX;\n"
			"}\n"
			"\n"
			"";
	std::string header3 = "\n"
			"void fftw(std::vector<std::complex<double>>& x) {\n"
			"\tconst int N = x.size();\n"
			"\tstatic std::unordered_map<int, fftw_plan> plans;\n"
			"\tstatic std::unordered_map<int, fftw_complex*> in;\n"
			"\tstatic std::unordered_map<int, fftw_complex*> out;\n"
			"\tif (plans.find(N) == plans.end()) {\n"
			"\t\tin[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);\n"
			"\t\tout[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);\n"
			"\t\tplans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_ESTIMATE);\n"
			"\t}\n"
			"\tauto* i = in[N];\n"
			"\tauto* o = out[N];\n"
			"\tfor (int n = 0; n < N; n++) {\n"
			"\t\ti[n][0] = x[n].real();\n"
			"\t\ti[n][1] = x[n].imag();\n"
			"\t}\n"
			"\tfftw_execute(plans[N]);\n"
			"\tfor (int n = 0; n < N; n++) {\n"
			"\t\tx[n].real(o[n][0]);\n"
			"\t\tx[n].imag(o[n][1]);\n"
			"\t}\n"
			"\n"
			"}\n"
			"";
	print("%s\n", header3.c_str());
	print("%s\n", header2.c_str());

	print("void test() {\n");
	indent();
	print("timer tm3, tm4;\n");
	print("for( int N = 2; N <= %i; N+=%i) {\n", MAXFFT, DFFT);
	indent();

	fprintf(fp, "\t\ttimer tm1, tm2;\n"
			"\t\tdouble err;\n"
			"\t\tdouble max;\n"
			"\t\tfor (int ti = 0; ti < 256; ti++) {\n"
			"\t\t\terr = 0.0;\n"
			"\t\t\tmax = 0.0;\n"
			"\t\t\tstd::vector<std::complex<double>> X(N);\n"
			"\t\t\tstd::vector<std::complex<double>> Y(N);\n"
			"\t\t\tfor (int n = 0; n < N; n++) {\n"
			"\t\t\t\tX[n] = std::complex<double>(rand1(), rand1());\n"
			"\t\t\t}\n"
			"\t\t\tauto X0 = X;\n"
			"\t\t\tY = X;\n"
			"\t\t\ttm1.start();\n"
			"\t\t\ttm3.start();\n"
			"\t\t\tFFT(X.data(), N);\n"
			"\t\t\ttm1.stop();\n"
			"\t\t\ttm3.stop();\n"
			"\t\t\ttm2.start();\n"
			"\t\t\ttm4.start();\n"
			"\t\t\tfftw(Y);\n"
			"\t\t\ttm2.stop();\n"
			"\t\t\ttm4.stop();\n"
			"\t\t\tfor (int i = 0; i < X.size(); i++) {\n"
			"\t\t\t\tY[i] -= X[i];\n"
			"\t\t\t}\n"
			"\t\t\tfor (int n = 0; n < N; n++) {\n"
			"\t\t\t\terr += std::abs(Y[n]) * std::abs(Y[n]);\n"
			"\t\t\t\tmax = std::max(max, std::abs(X0[n]));\n"
			"\t\t\t}\n"
			"\t\t\terr = sqrt(err / N) / max;\n"
			"\t\t}\n"
			"\t\tprintf(\"%%i %%e %%e %%e %%e %%e %%e %%e\\n\", N, err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());\n"
			"");
	deindent();
	print("}\n");
	deindent();
	print("}\n\n");

	print("int main() {\n");
	indent();
	print("test();\n");
	print("printf(\"\\n\");\n");
	print("test();\n");
	deindent();
	print("}\n\n");

	set_file("Makefile");
	print("CC=g++\n");
	print("CFLAGS=-I. -Ofast -march=native\n");
	print("DEPS = fft.hpp\n");
	print("OBJ = fft.o ");
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		fprintf(fp, "fft.%i.o ", n);
	}
	fprintf(fp, "\n");

	print("\n%.o: %.cpp $(DEPS)\n");
	print("\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");

	print("ffttest: $(OBJ)\n");
	print("\t$(CC) -o $@ $^ $(CFLAGS) -lfftw3\n");
	set_file("");
	return 0;
}
