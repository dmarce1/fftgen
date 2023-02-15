#define MAXFFTTEST MAXFFT

void FFT_real(std::complex<double>* xout, double* xin, int N);

#include <unordered_map>
#include <set>
#include <stack>

const std::vector<std::complex<double>>& twiddles(int N) {
	static std::unordered_map<int, std::vector<std::complex<double>>>values;
	if(values.find(N) == values.end()) {
		for( int n = 0; n < N; n++) {
			values[N].push_back(std::polar(1.0, -2.0*M_PI*n/N));
		}
	}
	return values[N];
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
	static thread_local std::unordered_map<int, int> values;
	auto i = values.find(N);
	if (i == values.end()) {
		for (int g = 2;; g++) {
			std::set<int> I;
			bool fail = false;
			for (int m = 0; m < N - 2; m++) {
				int n = mod_pow(g, m, N);
				if (I.find(n) == I.end()) {
					I.insert(n);
				} else {
					fail = true;
					break;
				}
			}
			if (!fail) {
				values[N] = g;
				i = values.find(N);
				break;
			}
		}
	}
	return i->second;
}

std::vector<int> raders_ginvq(int N) {
	static thread_local std::unordered_map<int, std::vector<int>> values;
	auto i = values.find(N);
	if (i == values.end()) {
		const int g = generator(N);
		std::vector<int> ginvq;
		for (int q = 0; q < N - 1; q++) {
			ginvq.push_back(mod_pow(g, q, N));
		}
		values[N] = std::move(ginvq);
		i = values.find(N);
	}
	return i->second;

}

const std::vector<int>& raders_gq(int N) {
	static thread_local std::unordered_map<int, std::vector<int>> values;
	auto i = values.find(N);
	if (i == values.end()) {
		const int g = generator(N);
		std::vector<int> gq;
		for (int q = 0; q < N - 1; q++) {
			gq.push_back(mod_inv(mod_pow(g, q, N), N));
		}
		values[N] = std::move(gq);
		i = values.find(N);
	}
	return i->second;
}

const std::vector<std::complex<double>>& raders_four_twiddle(int N) {
	static thread_local std::unordered_map<int, std::vector<std::complex<double>>>values;
	auto i = values.find(N);
	if (i == values.end()) {
		std::vector<std::complex<double>> b(N - 1);
		const auto& tws = twiddles(N);
		const auto& ginvq = raders_ginvq(N);
		for (int q = 0; q < N - 1; q++) {
			b[q] = tws[ginvq[q]];
		}
		FFT(b.data(), N - 1);
		values[N] = std::move(b);
		i = values.find(N);
	}
	return i->second;
}

void raders_fft(std::complex<double>* x, int N) {
	static thread_local std::stack<std::vector<std::complex<double>>>astack;
	std::vector<std::complex<double>> a;
	if( astack.size() ) {
		a = std::move(astack.top());
		astack.pop();
	}
	const auto& b = raders_four_twiddle(N);
	const auto& gq = raders_gq(N);
	const auto& ginvq = raders_ginvq(N);
	a.resize(N - 1);
	for (int q = 0; q < N - 1; q++) {
		a[q] = x[gq[q]];
	}
	FFT(a.data(), N - 1);
	for (int q = 0; q < N - 1; q++) {
		a[q] = std::conj(a[q] * b[q]);
	}
	FFT(a.data(), N - 1);
	const auto Nm1inv = 1.0 / (N - 1.0);
	for (int q = 0; q < N - 1; q++) {
		a[q] = std::conj(a[q]) * Nm1inv;
	}
	const auto x0 = x[0];
	for (int k = 1; k < N; k++) {
		x[0] += x[k];
		x[k] = x0;
	}
	for (int p = 0; p < N - 1; p++) {
		x[ginvq[p]] += a[p];
	}
	astack.push(std::move(a));
}

int cooley_tukey_radix(int N) {
	static thread_local std::unordered_map<int, int> radix;
	if (radix.find(N) == radix.end()) {
		int N2 = sqrt(N);
		while (N % N2 != 0) {
			N2++;
		}
		radix[N] = N / N2;
	}
	return radix[N];
}

const std::vector<int>& permutation_indices(int N, std::vector<int> I) {
	static thread_local std::unordered_map<int, std::vector<int>> map;
	if (map.find(N) == map.end()) {
		if (N < MAXFFT) {
			const auto& J = (*(pfptr[N]))();
			std::vector<int> K(N);
			for (int i = 0; i < N; i++) {
				K[i] = I[J[i]];
			}
			map[N] = K;
		} else {
			int N1 = cooley_tukey_radix(N);
			int N2 = N / N1;
			std::vector<int> L;
			for (int n1 = 0; n1 < N1; n1++) {
				std::vector<int> J;
				for (int n2 = 0; n2 < n2; N2++) {
					J.push_back(I[N2 * n1 + n2]);
				}
				J = permutation_indices(N2, J);
				L.insert(L.end(), J.begin(), J.end());
			}
			map[N] = L;
		}
	} else {
		return map[N];
	}
}

void cooley_tukey_base(std::complex<double>* x, int N) {

}

void cooley_tukey_fft(std::complex<double>* x, int N) {
	static thread_local std::stack<std::vector<std::complex<double>>>astack;
	int N1 = cooley_tukey_radix(N);
	std::vector<std::complex<double>> y;
	if( astack.size() ) {
		y = std::move(astack.top());
		astack.pop();
	}
	y.resize(N);
	int N2 = N / N1;
	for (int n = 0; n < N; n++) {
		y[n] = x[n];
	}
	for (int n1 = 0; n1 < N1; n1++) {
		FFT(y.data() + n1, N2);
	}
	const auto& W = twiddles(N);
	for (int k2 = 1; k2 < N2; k2++) {
		for (int n1 = 1; n1 < N1; n1++) {
			y[n1 + k2 * N1] *= W[n1 * k2];
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		FFT(y.data() + k2 * N1, N1);
		for (int k1 = 0; k1 < N1; k1++) {
			x[k2 + k1 * N2] = y[k2 * N1 + k1];
		}
	}
	astack.push(std::move(y));
}

void FFT_large(std::complex<double>* x, int N) {
	if (greatest_prime_factor(N) > MAXFFT && greatest_prime_factor(N) == N) {
		raders_fft(x, N);
	} else {
		cooley_tukey_fft(x, N);
	}
}

void FFT_sym_real(double* x, int N) {
	std::vector<double> A(N / 4 + 1);
	std::vector<double> B(N / 4 + 1);
	std::vector<double> y(N / 2);
	std::vector<std::complex<double>> Y(N / 4 + 1);
	const auto index = [N](int n) {
		return (n + N) % N;
	};
	for (int n = 1; n < N / 4; n++) {
		y[n] = x[index(2 * n)] + (x[index(2 * n + 1)] - x[index(2 * n - 1)]);
		y[N / 2 - n] = x[index(2 * n)] - (x[index(2 * n + 1)] - x[index(2 * n - 1)]);
	}
	y[0] = x[0];
	y[N / 4] = x[N / 2];
	B[0] = 0.0;
	for (int n = 0; n < N / 4; n++) {
		B[0] += 2.0 * x[2 * n + 1];
	}
	FFT_real(Y.data(), y.data(), N / 2);
	for (int k = 0; k <= N / 4; k++) {
		A[k] = Y[k].real();
	}
	for (int k = 1; k <= N / 4; k++) {
		B[k] = Y[k].imag() / (2.0 * sin(2.0 * M_PI * k / N));
	}
	for (int k = 1; k <= N / 4; k++) {
		x[k] = A[k] + B[k];
		x[N / 2 - k] = A[k] - B[k];
	}
	x[0] = A[0] + B[0];
	x[N / 2] = A[0] - B[0];
}

void FFT_antsym_real(double* x, int N) {
	std::vector<double> A(N / 4 + 1);
	std::vector<double> B(N / 4 + 1);
	std::vector<double> y(N / 2);
	std::vector<std::complex<double>> Y(N / 4 + 1);
	const auto index = [N](int n) {
		return (n + N) % N;
	};
	for (int n = 1; n < N / 4; n++) {
		y[n] = x[index(2 * n)] + (x[index(2 * n + 1)] - x[index(2 * n - 1)]);
		y[N / 2 - n] = -x[index(2 * n)] + (x[index(2 * n + 1)] - x[index(2 * n - 1)]);
	}
	y[0] = x[0];
	y[N / 4] = -x[N / 2];
	B[0] = 0.0;
	for (int n = 0; n < N / 4; n++) {
		B[0] += 2.0 * x[2 * n + 1];
	}
	FFT_real(Y.data(), y.data(), N / 2);
	for (int k = 0; k <= N / 4; k++) {
		A[k] = Y[k].imag();
	}
	for (int k = 1; k <= N / 4; k++) {
		B[k] = Y[k].real() / (2.0 * sin(2.0 * M_PI * k / N));
	}
	for (int k = 1; k <= N / 4; k++) {
		x[k] = A[k] - B[k];
		x[N / 2 - k] = -(A[k] + B[k]);
	}
	x[0] = 0;
	x[N / 2] = -(A[0] + B[0]);
}

void FFT_real(std::complex<double>* xout, double* xin, int N) {
	(*(fptr_real[N]))(xin, reinterpret_cast<double*>(xout));
}

void FFT(std::complex<double>* x, int N) {
	if (N > MAXFFT) {
		FFT_large(x, N);
	} else {
		(*(fptr[N]))(reinterpret_cast<double*>(x));
	}
}

void FFT_inv(std::complex<double>* x, int N) {
	const double Ninv = 1.0 / N;
	for (int n = 1; n < N - n; n++) {
		std::swap(x[n], x[N - n]);
	}
	(*(fptr[N]))(reinterpret_cast<double*>(x));
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
void fftw_dct(std::vector<double>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, double*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (double*) malloc(sizeof(double) * N);
		plans[N] = fftw_plan_r2r_1d(N, in[N], out[N], FFTW_REDFT10, FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = x[n];
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n] = (o[n]);
	}

}

void fftw_dst(std::vector<double>& x) {
	const int N = x.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, double*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (double*) malloc(sizeof(double) * N);
		plans[N] = fftw_plan_r2r_1d(N, in[N], out[N], FFTW_RODFT10, FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = x[n];
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		x[n] = (o[n]);
	}

}

void fftw_real(std::vector<std::complex<double>>& xout, const std::vector<double>& xin) {
	const int N = xin.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));
		plans[N] = fftw_plan_dft_r2c_1d(N, in[N], out[N], FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N; n++) {
		i[n] = xin[n];
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N / 2 + 1; n++) {
		xout[n].real(o[n][0]);
		xout[n].imag(o[n][1]);
	}

}

void fftw_real_inv(const std::vector<std::complex<double>>& xout, std::vector<double>& xin) {
	const int N = xin.size();
	static std::unordered_map<int, fftw_plan> plans;
	static std::unordered_map<int, double*> in;
	static std::unordered_map<int, fftw_complex*> out;
	if (plans.find(N) == plans.end()) {
		in[N] = (double*) malloc(sizeof(double) * N);
		out[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));
		plans[N] = fftw_plan_dft_c2r_1d(N, out[N], in[N], FFTW_ESTIMATE);
	}
	auto* i = in[N];
	auto* o = out[N];
	for (int n = 0; n < N / 2 + 1; n++) {
		o[n][0] = xout[n].real();
		o[n][1] = xout[n].imag();
	}
	o[0][1] = 0;
	if (N % 2 == 0) {
		o[N / 2][1] = 0;
	}
	fftw_execute(plans[N]);
	for (int n = 0; n < N; n++) {
		xin[n] = i[n];
	}

}

void good_thomas(int N1, int N2, std::vector<std::complex<double>>& X) {
	int N = N1 * N2;
	std::vector<std::vector<std::complex<double>>>Y(N1, std::vector<std::complex<double>>(N2));
	std::vector<std::vector<std::complex<double>>>Z(N2, std::vector<std::complex<double>>(N1));
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Y[n1][n2] = X[(N1 * n2 + N2 * n1) % N];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fftw(Y[n1]);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			Z[n2][n1] = Y[n1][n2];
		}
	}
	for (int n2 = 0; n2 < N2; n2++) {
		fftw(Z[n2]);
	}
	for (int n = 0; n < N; n++) {
		int n1 = n % N1;
		int n2 = n % N2;
		X[n] = Z[n2][n1];
	}
}

void good_thomas(int N1, int N2, std::vector<std::complex<double>>& X, std::vector<double> x) {
	int N = N1 * N2;
	std::vector<std::vector<double>> y(N1, std::vector<double>(N2));
	std::vector<std::vector<std::complex<double>>>Y(N1, std::vector<std::complex<double>>(N2 / 2 + 1));
	std::vector<double> z(N1);
	std::vector<std::vector<std::complex<double>>>Z(N2 / 2 + 1, std::vector<std::complex<double>>(N1));
	std::vector<std::complex<double>> Z1(N1 / 2 + 1);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			y[n1][n2] = x[(N1 * n2 + N2 * n1) % N];
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		fftw_real(Y[n1], y[n1]);
	}
	for (int n1 = 0; n1 < N1; n1++) {
		z[n1] = Y[n1][0].real();
		for (int n2 = 1; n2 < N2 / 2 + 1; n2++) {
			Z[n2 - 1][n1] = Y[n1][n2];
		}
	}
	fftw_real(Z1, z);
	for (int n2 = 1; n2 < N2 / 2 + 1; n2++) {
		fftw(Z[n2 - 1]);
	}
	for (int n = 0; n < N / 2 + 1; n++) {
		int n1 = n % N1;
		int n2 = n % N2;
		if (n2 == 0) {
			if (n1 > N1 / 2) {
				X[n] = std::conj(Z1[N1 - n1]);
			} else {
				X[n] = Z1[n1];
			}
		} else if (n2 > N2 / 2) {
			if (n1 != 0) {
				X[n] = std::conj(Z[N2 - n2 - 1][N1 - n1]);
			} else {
				X[n] = std::conj(Z[N2 - n2 - 1][0]);
			}
		} else {
			X[n] = Z[n2 - 1][n1];
		}

	}
}

double rand1() {
	return (rand() + 0.5) / RAND_MAX;
}

void FCT(std::vector<double>& x) {
	int N = x.size();
	std::vector<double> y(N);
	std::vector<std::complex<double>> Y(N);
	int M = N / 4;
	if (N % 4 == 1) {
		for (int n = 0; n <= (N - 1) / 2; n++) {
			y[(n - M + N) % N] = x[2 * n];
		}
		for (int n = 0; n <= (N - 3) / 2; n++) {
			y[(3 * M - n + N) % N] = x[2 * n + 1];
		}
	} else if (N % 4 == 3) {
		for (int n = 0; n <= (N - 3) / 2; n++) {
			y[(n - M + N) % N] = x[2 * n + 1];
		}
		for (int n = 0; n <= (N - 1) / 2; n++) {
			y[(3 * M - n + N + 2) % N] = x[2 * n];
		}
	} else {
		for (int n = 0; n < N / 2; n++) {
			y[n] = x[2 * n];
			y[N - 1 - n] = x[2 * n + 1];
		}
	}
	FFT_real(Y.data(), y.data(), N);
	if (N % 2 == 1) {
		for (int n = 0; n <= (N - 1) / 2; n++) {
			double sgn = (n / 2) % 2 == 0 ? 1.0 : -1.0;
			if (n % 2 == 1) {
				x[n] = 2 * sgn * Y[n].imag();
			} else {
				x[n] = 2 * sgn * Y[n].real();
			}
		}
		for (int n = 1; n <= (N - 1) / 2; n++) {
			double sgn = ((n + 1) / 2) % 2 == 0 ? -1.0 : +1.0;
			if (N % 4 == 3) {
				sgn = -sgn;
			}
			if (n % 2 == 0) {
				x[N - n] = 2 * sgn * Y[n].imag();
			} else {
				x[N - n] = 2 * sgn * Y[n].real();
			}
		}
	} else {
		for (int k = 0; k < N / 2; k++) {
			x[k] = 2.0 * (Y[k].real() * cos(0.5 * M_PI * k / N) + Y[k].imag() * sin(0.5 * M_PI * k / N));
		}
		for (int k = N / 2; k < N; k++) {
			x[k] = 2.0 * (Y[N - k].real() * cos(0.5 * M_PI * k / N) - Y[N - k].imag() * sin(0.5 * M_PI * k / N));
		}
	}
}

void FST(std::vector<double>& x) {
	int N = x.size();
	for (int n = 0; n < N / 2; n++) {
		std::swap(x[n], x[N - n - 1]);
	}
	for (int n = 1; n < N; n += 2) {
		x[n] = -x[n];
	}
	FCT(x);
	for (int n = 0; n < N / 2; n++) {
		std::swap(x[n], x[N - n - 1]);
	}
	for (int n = 1; n < N; n += 2) {
		x[n] = -x[n];
	}
}

void test() {
	timer tm3, tm4;

	tm3.reset();
	tm4.reset();
	tm3.reset();
	tm4.reset();

	timer tm1, tm2;
	double err;
	double max;
	for (int N = 2; N < 64; N++) {
		for (int ti = 0; ti < 1; ti++) {
			for (int ti = 0; ti < 256; ti++) {
				err = 0.0;
				max = 0.0;
				std::vector<double> X(N);
				std::vector<double> Y(N);
				for (int n = 0; n < N; n++) {
					X[n] = double(rand1());
				}
				auto X0 = X;
				Y = X;
				tm1.start();
				tm3.start();
				FST(X);
				tm1.stop();
				tm3.stop();
				tm2.start();
				tm4.start();
				fftw_dst(Y);
				tm2.stop();
				tm4.stop();
				for (int i = 0; i < X.size(); i++) {
					Y[i] -= X[i];
				}
				for (int n = 0; n < N; n++) {
					err += std::abs(Y[n]) * std::abs(Y[n]);
					max = std::max(max, std::abs(X0[n]));
					//	printf("%i %e %e\n", n, X[n], Y[n]);
				}
				err = sqrt(err / N) / max;
			}
		}
		printf("%4i %4i %e %e %e %e %e %e %e\n", N, fft_nops[N], err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());
	}
//	return;

	printf("\ncomplex\n");
	for (int N = 2; N <= MAXFFTTEST; N++) {
		timer tm1, tm2;
		double err;
		double max;
		for (int ti = 0; ti < 256; ti++) {
			err = 0.0;
			max = 0.0;
			std::vector<std::complex<double>> X(N);
			std::vector<std::complex<double>> Y(N);
			for (int n = 0; n < N; n++) {
				X[n] = std::complex<double>(rand1(), rand1());
			}
			auto X0 = X;
			Y = X;
			tm1.start();
			tm3.start();
			FFT(X.data(), N);
			tm1.stop();
			tm3.stop();
			tm2.start();
			tm4.start();
			fftw(Y);
			tm2.stop();
			tm4.stop();
			for (int i = 0; i < X.size(); i++) {
				Y[i] -= X[i];
			}
			for (int n = 0; n < N; n++) {
				err += std::abs(Y[n]) * std::abs(Y[n]);
				max = std::max(max, std::abs(X0[n]));
				//printf("%i %e %e %e %e\n", n, X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag());
			}
			err = sqrt(err / N) / max;
		}
		printf("%4i %e %e %e %e %e %e %e\n", N, err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());
	}
	tm3.reset();
	tm4.reset();
	//return;
	printf("\nreal\n");
	for (int N = 2; N <= MAXFFTTEST; N++) {
		if (!allow_real(N)) {
			continue;
		}
		timer tm1, tm2;
		double err;
		double max;
		for (int ti = 0; ti < 256; ti++) {
			err = 0.0;
			max = 0.0;
			std::vector<std::complex<double>> X(N / 2 + 1);
			std::vector<std::complex<double>> Y(N / 2 + 1);
			std::vector<double> y(N);
			std::vector<double> x(N);
			for (int n = 0; n < N; n++) {
				x[n] = rand1();
			}
			y = x;
			tm1.start();
			tm3.start();
			FFT_real(X.data(), x.data(), N);
			tm1.stop();
			tm3.stop();
			tm2.start();
			tm4.start();
			fftw_real(Y, y);
			tm2.stop();
			tm4.stop();
			//	printf("\n");
			for (int i = 0; i < X.size(); i++) {
				Y[i] -= X[i];
			}
			for (int n = 0; n < N / 2 + 1; n++) {
				err += std::abs(Y[n]) * std::abs(Y[n]);
				//	printf("%i %16.6e %16.6e %16.6e %16.6e\n", n, X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag());
				max = std::max(max, std::abs(X[n]));
			}
			err = sqrt(err / N) / max;
		}
		printf("%4i %4i %e %e %e %e %e %e %e\n", N, fft_nops_real[N], err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());
	}
}

int main() {
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_OVERFLOW);
	feenableexcept(FE_INVALID);
	test();
	printf("\n");
	test();
	printf("EXITING\n");
}

