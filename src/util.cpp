#include "fftgen.hpp"
#include <cassert>
#include <fftw3.h>
#include <unordered_map>
#include <vector>
#include <algorithm>

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

int greatest_fft_factor(int N) {
	int gprime = greatest_prime_factor(N);
	if (gprime > 15) {
		return gprime;
	} else {
		for (int i = 15; i >= 2; i--) {
			if (N % i == 0) {
				return i;
			}
		}
	}

}

std::vector<std::pair<int, int>> fft_factorization(int N) {
	std::unordered_map<int, int> facs;
	while (N != 1) {
		int pf = greatest_fft_factor(N);
		if (facs.find(pf) == facs.end()) {
			facs[pf] = 0;
		}
		facs[pf]++;
		N /= pf;
	}
	std::vector<std::pair<int, int>> rc(facs.begin(), facs.end());
	std::sort(rc.begin(), rc.end(), [](std::pair<int, int> a, std::pair<int, int> b ) {
		return a.first < b.first;
	});
	return rc;
}


std::vector<std::pair<int, int>> prime_factorization(int N) {
	std::unordered_map<int, int> facs;
	while (N != 1) {
		int pf = greatest_prime_factor(N);
		if (facs.find(pf) == facs.end()) {
			facs[pf] = 0;
		}
		facs[pf]++;
		N /= pf;
	}
	std::vector<std::pair<int, int>> rc(facs.begin(), facs.end());
	std::sort(rc.begin(), rc.end(), [](std::pair<int, int> a, std::pair<int, int> b ) {
		return a.first < b.first;
	});
	return rc;
}

std::string fft_factorization_string(int N) {
	const auto facs = fft_factorization(N);
	std::string str;
	for (int n = 0; n < facs.size(); n++) {
		str += std::to_string(facs[n].first);
		str += "^";
		str += std::to_string(facs[n].second);
		if (n != facs.size() - 1) {
			str += " * ";
		}
	}
	return str;
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

const std::vector<std::complex<double>> twiddles(int N) {
	std::vector<std::complex<double>> tw(N);
	for (int k = 0; k < N; k++) {
		tw[k] = std::polar(1.0, -2 * M_PI * k / N);
	}
	return tw;
}

int index(int o, int i, int j, int N) {
	return 2 * (o + i) + j;
}

int index_real(int o, int i, int j, int N) {
	if (j == 0) {
		return o + i;
	} else {
		return o + (N - i);
	}
}

int index_real_inv(int o, int i, int j, int N) {
	if (i > N / 2) {
		return index_real_inv(o, N - i, j, N);
	}
	if (j == 0) {
		return o + i;
	} else {
		return o + (N - i);
	}
}

long long mod_pow(long long a, long long b, long long m) {
	long long rc = 1;
	long long apow = a;
	while (b) {
		if (b & (long long) 1) {
			rc = ((rc % m) * (apow % m)) % m;
		}
		b >>= (long long) 1;
		apow = ((apow % m) * (apow % m)) % m;
	}
	return rc;
}

int mod_inv(int a, int m) {
	return mod_pow(a, m - 2, m);
}

int generator(long long N) {
	for (long long g = 2; g < 1000; g++) {
		std::vector<bool> I(N, false);
		bool fail = false;
		for (long long m = 0; m < N - (long long) 1; m++) {
			long long n = mod_pow(g, m, N);
			assert(n >= 0);
			assert(n < N);
			if (!I[n]) {
				if (N == 839 && g == 11) {
					printf("%i %i\n", m, n);
				}
				I[n] = true;
			} else {
				if (N == 839 && g == 11) {
					printf("%i %i * \n", m, n);
				}
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
