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
	print("// Raders radix - %i x %i\n", N1, N2);
//	printf("// Raders radix - %i \n", r);
	print("std::array<std::complex<double>, %i> y;\n", N);
	print("std::array<std::complex<double>, %i> xo, xk0;\n", N2);
	print("std::complex<double> tmp;\n");
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
				print("x[%i] *= %s;\n", o + k2 * N1 + n1, to_str(twiddle(nk, N)).c_str());
			}
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		print("xo[%i] = x[%i];\n", k2, o + k2 * N1);
		print("xk0[%i] = x[%i] ", k2, o + k2 * N1);
		for (int k1 = 1; k1 < N1; k1++) {
			fprintf(fp, " + x[%i]", o + (N1 * k2 + k1));
		}
		fprintf(fp, ";\n");
	}

	std::vector<int> I1(N);
	for (int k2 = 0; k2 < N2; k2++) {
		auto I2 = fft_bitreverse_indices(N1 - 1);
		auto I3 = I2;
		for (int i = 0; i < I2.size(); i++) {
			I3[I2[i]] = i;
		}
		for (int q = 0; q < N1 - 1; q++) {
//			print("y[%i] = x[%i];\n", N1 * k2 + q, o + (N1 * k2 + gq[q]));
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
		print("x[%i] = std::conj(-x[%i]);\n", o + k2 * N1, o + k2 * N1);
		for (int q = 1; q < N1 - 1; q++) {
			print("x[%i] = std::conj(x[%i] * %s);\n", o + (N1 * k2 + q), (o + N1 * k2 + q), to_str(b[q]).c_str());
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
			print("x[%i] = std::conj(x[%i]) * (%24.17e);\n", o + N1 * k2 + q, o + N1 * k2 + q, Nm1inv);
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int p = N1 - 2; p >= 0; p--) {
			print("x[%i] = xo[%i] + x[%i];\n", o + (N1 * k2 + p + 1), k2, o + N1 * k2 + p);
		}
		print("x[%i] = xk0[%i];\n", o + k2 * N1, k2);
	}
	std::vector<int> indices;
	indices.resize(N);
	for (int k2 = 0; k2 < N2; k2++) {
		indices[k2] = N1 * k2;
		for (int k1 = 0; k1 < N1 - 1; k1++) {
			indices[N2 * ginvq[k1] + k2] = N1 * k2 + k1 + 1;
		}
	}
	std::vector<int> I;
	std::vector<bool> touched(N, false);
	I.resize(0);
	for (int n = 0; n < N; n++) {
		I.push_back(n);
	}
	for (int k = 0; k < I.size(); k++) {
		if (!touched[k]) {
			int current = k;
			int first = current;
			int next = indices[current];
			touched[current] = true;
			if (first != next) {
				print("tmp = x[%i];\n", o + current);
				while (first != next) {
					print("x[%i] = x[%i];\n", o + current, o + next);
					current = next;
					touched[current] = true;
					next = indices[next];
				}
				print("x[%i] = tmp;\n", o + current);
			}
		}
	}
	deindent();
	print("}\n");
}

int raders_fft_opcnt(int r, int N) {
	int cnt = 0;
	const int N1 = r;
	const int N2 = N / r;
	if (N2 > 1) {
		cnt += N;
		cnt += N1 * fft_opcnt(N2);
		cnt += N;
	}
	cnt += 6 * (N - 1);
	for (int k2 = 0; k2 < N2; k2++) {
		cnt += 5;
		cnt += 2 * fft_opcnt(N1 - 1);
		cnt += 12 * (N1 - 1);
		cnt += 2 * N1;
	}
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
		return 6;
	}
}

void print_z(int zi, int k, int r, int N, int o) {
	const auto W = twiddle(zi * k, N);
	const int index = o + (k + zi * N / r);
	if (zi * k == 0) {
		print("const auto z%i = x[%i];\n", zi, index);
	} else if (zi * k == N / 8 && N % 8 == 0) {
		print("const auto z%i = No8x(x[%i]);\n", zi, index);
	} else if (zi * k == N / 4 && N % 4 == 0) {
		print("const auto z%i = nIx(x[%i]);\n", zi, index);
	} else if (zi * k == 3 * N / 8 && N % 8 == 0) {
		print("const auto z%i = nNo8xIx(x[%i]);\n", zi, index);
	} else if (zi * k == N / 2 && N % 2 == 0) {
		print("const auto z%i = -x[%i];\n", zi, index);
	} else if (zi * k == 5 * N / 8 && N % 8 == 0) {
		print("const auto z%i = nNo8x(x[%i]);\n", zi, index);
	} else if (zi * k == 3 * N / 4 && N % 4 == 0) {
		print("const auto z%i = Ix(x[%i]);\n", zi, index);
	} else if (zi * k == 7 * N / 8 && N % 8 == 0) {
		print("const auto z%i = No8xIx(x[%i]);\n", zi, index);
	} else {
		print("const auto z%i = %s * x[%i];\n", zi, to_str(W).c_str(), index);
	}
}

int best_radix(int N, int o) {
	int best_cnt = 999999999;
	int best_radix = -1;
	for (int r = 2; r <= N; r++) {
		if (N % r == 0) {
			int this_cnt;
			if (r <= 6 || is_prime(r)) {
				this_cnt = fft_radix_opcnt(r, N);
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
				const int this_cnt = fft_radix_opcnt(r, N);
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
			cnt += 20;
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
			cnt += 18;
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
			cnt += 48;
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
			cnt += 48;
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

			for (int j = 1; j <= (r - 1) / 2; j++) {
				cnt += 4;
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				for (int j = 1; j <= (r - 1) / 2; j++) {
					cnt += 4;
				}
			}

			for (int i = 1; i <= (r - 1) / 2; i++) {
				for (int j = 1; j <= (r - 1) / 2; j++) {
					if (j > 1) {
						cnt += 2;
					}
					cnt += 2;
				}
			}
			for (int i = 1; i < r; i++) {
				cnt++;
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				cnt += 14;
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

void fft_bitreverse(int N, std::vector<int> indices, int o) {
	if (indices.size() == 0) {
		indices = fft_bitreverse_indices(N);
	}
	std::vector<bool> touched(N, false);
	for (int k = 0; k < N; k++) {
		if (!touched[k]) {
			int current = k;
			int first = current;
			int next = indices[current];
			touched[current] = true;
			if (first != next) {
				print("tmp = x[%i];\n", o + current);
				while (first != next) {
					print("x[%i] = x[%i];\n", o + current, o + next);
					current = next;
					touched[current] = true;
					next = indices[next];
				}
				print("x[%i] = tmp;\n", o + current);
			}
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
//	} else if (radix == -1) {
//		bluestein_fft(N, o);
//	}
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

int index_lookup(const std::vector<int>& I, int i) {
	for (int j = 0; j < I.size(); j++) {
		if (i == I[j]) {
			return j;
		}
	}
	return -1;
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
			print("x[%i] = z0 + z1;\n", o + k);
			print("x[%i] = z0 - z1;\n", o + (k + N / 2));
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
			print("const auto t1 = z1 + z2;\n");
			print("const auto t2 = z0 - t1 * 0.5;\n");
			print("const auto t3 = (%24.17e) * (z1 - z2);\n", sqrt(3) * 0.5);
			print("x[%i] = z0 + t1;\n", o + k);
			print("x[%i] = t2 - Ix(t3);\n", o + (k + N / 3));
			print("x[%i] = t2 + Ix(t3);\n", o + (k + 2 * N / 3));
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
			print("const auto t1 = z0 + z2;\n");
			print("const auto t2 = z1 + z3;\n");
			print("const auto t3 = z0 - z2;\n");
			print("const auto t4 = z1 - z3;\n");
			print("x[%i] = t1 + t2;\n", o + k);
			print("x[%i] = t3 - Ix(t4);\n", o + (k + N / 4));
			print("x[%i] = t1 - t2;\n", o + (k + N / 2));
			print("x[%i] = t3 + Ix(t4);\n", o + (k + 3 * N / 4));
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
			print("const auto t1 = z1 + z4;\n");
			print("const auto t2 = z2 + z3;\n");
			print("const auto t3 = z1 - z4;\n");
			print("const auto t4 = z2 - z3;\n");
			print("const auto t5 = t1 + t2;\n");
			print("const auto t6 = (%24.17e) * (t1 - t2);\n", sqrt(5) * 0.25);
			print("const auto t7 = z0 - t5 * 0.25;\n");
			print("const auto t8 = t7 + t6;\n");
			print("const auto t9 = t7 - t6;\n");
			print("const auto t10 = (%24.17e) * t3 + (%24.17e) * t4;\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
			print("const auto t11 = (%24.17e) * t3 - (%24.17e) * t4;\n", sin(2.0 * M_PI / 10.0), sin(2.0 * M_PI / 5.0));
			print("x[%i] = z0 + t5;\n", o + (k + 0 * N / 5));
			print("x[%i] = t8 - Ix(t10);\n", o + (k + 1 * N / 5));
			print("x[%i] = t9 - Ix(t11);\n", o + (k + 2 * N / 5));
			print("x[%i] = t9 + Ix(t11);\n", o + (k + 3 * N / 5));
			print("x[%i] = t8 + Ix(t10);\n", o + (k + 4 * N / 5));
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
			print("const auto t1 = z2 + z4;\n");
			print("const auto t2 = z0 - t1 * 0.5;\n");
			print("const auto t3 = (%24.17e) * (z2 - z4);\n", sin(M_PI / 3.0));
			print("const auto t4 = z5 + z1;\n");
			print("const auto t5 = z3 - t4 * 0.5;\n");
			print("const auto t6 = (%24.17e) * (z5 - z1);\n", sin(M_PI / 3.0));
			print("const auto t7 = z0 + t1;\n");
			print("const auto t8 = t2 - Ix(t3);\n");
			print("const auto t9 = t2 + Ix(t3);\n");
			print("const auto t10 = z3 + t4;\n");
			print("const auto t11 = t5 - Ix(t6);\n");
			print("const auto t12 = t5 + Ix(t6);\n");
			print("x[%i] = t7 + t10;\n", o + (k + 0 * N / 6));
			print("x[%i] = t8 - t11;\n", o + (k + 1 * N / 6));
			print("x[%i] = t9 + t12;\n", o + (k + 2 * N / 6));
			print("x[%i] = t7 - t10;\n", o + (k + 3 * N / 6));
			print("x[%i] = t8 + t11;\n", o + (k + 4 * N / 6));
			print("x[%i] = t9 - t12;\n", o + (k + 5 * N / 6));
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
					print("const auto txp%i = z%i.real() + z%i.real();\n", j, j, r - j);
				}
				for (int j = 1; j <= (r - 1) / 2; j++) {
					print("const auto txm%i = z%i.real() - z%i.real();\n", j, j, r - j);
				}
				for (int j = 1; j <= (r - 1) / 2; j++) {
					print("const auto typ%i = z%i.imag() + z%i.imag();\n", j, j, r - j);
				}
				for (int j = 1; j <= (r - 1) / 2; j++) {
					print("const auto tym%i = z%i.imag() - z%i.imag();\n", j, j, r - j);
				}
				for (int i = 1; i <= (r - 1) / 2; i++) {
					print("const auto ap%i = z0.real()", i);
					for (int j = 1; j <= (r - 1) / 2; j++) {
						if (!null) {
							fprintf(fp, " + ");
							fprintf(fp, "txp%i * (%24.17e)", j, cos(2.0 * M_PI * j * i / r));
						}
					}
					if (!null) {
						fprintf(fp, ";\n");
					}
				}
				for (int i = 1; i <= (r - 1) / 2; i++) {
					print("const auto bp%i = z0.imag()", i);
					for (int j = 1; j <= (r - 1) / 2; j++) {
						if (!null) {
							fprintf(fp, " + ");
							fprintf(fp, "typ%i * (%24.17e)", j, cos(2.0 * M_PI * j * i / r));
						}
					}
					if (!null) {
						fprintf(fp, ";\n");
					}
				}
				for (int i = 1; i <= (r - 1) / 2; i++) {
					print("const auto am%i = ", i);
					for (int j = 1; j <= (r - 1) / 2; j++) {
						if (!null) {
							if (j > 1) {
								fprintf(fp, " + ");
							}
							fprintf(fp, "tym%i * (%24.17e)", j, sin(2.0 * M_PI * j * i / r));
						}
					}
					if (!null) {
						fprintf(fp, ";\n");
					}
				}
				for (int i = 1; i <= (r - 1) / 2; i++) {
					print("const auto bm%i = ", i);
					for (int j = 1; j <= (r - 1) / 2; j++) {
						if (!null) {
							if (j > 1) {
								fprintf(fp, " + ");
							}
							fprintf(fp, "txm%i * (%24.17e)", j, sin(2.0 * M_PI * j * i / r));
						}
					}
					if (!null) {
						fprintf(fp, ";\n");
					}
				}
				print("x[%i] = ", o + k);
				for (int i = 0; i < r; i++) {
					if (!null) {
						fprintf(fp, " + z%i", i);
					}
				}
				if (!null) {
					fprintf(fp, ";\n");
				}
				for (int i = 1; i <= (r - 1) / 2; i++) {
					print("x[%i] = (ap%i + am%i) + Ix(bp%i - bm%i);\n", o + (k + i * N / r), i, i, i, i);
					print("x[%i] = (ap%i - am%i) + Ix(bp%i + bm%i);\n", o + (k + (r - i) * N / r), i, i, i, i);
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
	print("\nvoid fft_%i(std::complex<double>* x) {\n", N);
	indent();
	print("std::complex<double> tmp;\n");
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
