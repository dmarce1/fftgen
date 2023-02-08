#define MAXFFT 128

void FFT_real(std::complex<double>* xout, double* xin, int N) {
	(*(fptr_real[N]))(xin, reinterpret_cast<double*>(xout));
}

void FFT_real_inv(std::complex<double>* X, double* y, int N) {
	y[0] = X[0].real();
	for (int n = 1; n < N; n++) {
		if (n <= N - n) {
			y[n] = X[n].real() + X[n].imag();
		} else {
			y[n] = X[N - n].real() - X[N - n].imag();
		}
	}
	if (N % 2 == 0) {
		y[N / 2] = X[N / 2].real();
	}
	FFT_real(X, y, N);
	for (int n = 0; n < N; n++) {
		if (n <= N - n) {
			y[n] = X[n].real() + X[n].imag();
		} else {
			y[n] = X[N - n].real() - X[N - n].imag();
		}
	}
}

void FFT(std::complex<double>* x, int N) {
	(*(fptr[N]))(x);
}

void FFT_inv(std::complex<double>* x, int N) {
	const double Ninv = 1.0 / N;
	for (int n = 0; n < N; n++) {
		x[n] = std::conj(x[n]);
	}
	(*(fptr[N]))(x);
	for (int n = 0; n < N; n++) {
		x[n] = std::conj(x[n]) * Ninv;
	}
}

void FFT_real_sym(double* x, int N) {
	std::vector<double> A(N / 4 + 1);
	std::vector<double> B(N / 4 + 1);
	std::vector<double> y(N / 2);
	std::vector<std::complex<double>> Y(N / 4 + 1);
	const auto index = [N](int n) {
		n = (n + N) % N;
		return std::min(n, N - n);
	};
	for (int n = 0; n < N / 4; n++) {
		B[0] += x[index(2 * n + 1)];
	}
	B[0] *= 2.0;
	for (int n = 1; n < N / 4; n++) {
		y[n] = x[index(2 * n)] + (x[index(2 * n + 1)] - x[index(2 * n - 1)]);
	}
	y[0] = x[0];
	if ((N / 2) % 2 == 0) {
		y[N / 4] = x[N / 2];
	}
	for (int n = 1; n < N / 4; n++) {
		y[N / 2 - n] = x[index(2 * n)] - (x[index(2 * n + 1)] - x[index(2 * n - 1)]);
	}
	FFT_real(Y.data(), y.data(), N / 2);
	for (int k = 0; k <= N / 4; k++) {
		A[k] = Y[k].real();
	}
	for (int k = 1; k <= N / 4; k++) {
		B[k] = Y[k].imag() / (2.0 * sin(2.0 * M_PI * k / N));
	}
	for (int k = 0; k <= N / 4; k++) {
		x[k] = A[k] + B[k];
	}
	for (int k = 1; k <= N / 4; k++) {
		x[N / 2 - k] = A[k] - B[k];
	}
	if ((N / 2) % 2 == 0) {
		x[N / 2] = A[0] - B[0];
	}
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

double rand1() {
	return (rand() + 0.5) / RAND_MAX;
}

void test() {
	timer tm3, tm4;

	for (int N = 2; N <= MAXFFT; N += 1) {
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
		printf("%4i %4i %e %e %e %e %e %e %e\n", N, fft_nops[N], err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());
	}
	tm3.reset();
	tm4.reset();
	for (int N = 2; N <= MAXFFT; N += 1) {
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
			for (int n = 0; n < N / 2 + 1; n++) {
				X[n] = std::complex<double>(rand1(), rand1());
			}
			Y = X;
			tm1.start();
			tm3.start();
			FFT_real_inv(X.data(), x.data(), N);
			tm1.stop();
			tm3.stop();
			tm2.start();
			tm4.start();
			fftw_real_inv(Y, y);
			tm2.stop();
			tm4.stop();
			for (int i = 0; i < x.size(); i++) {
				y[i] -= x[i];
			}
			for (int n = 0; n < N; n++) {
				err += std::abs(y[n]) * std::abs(y[n]);
				//		printf("%i %e %e\n", n, x[n], y[n]);
				max = std::max(max, std::abs(x[n]));
			}
			err = sqrt(err / N) / max;
		}
		printf("%4i %4i %e %e %e %e %e %e %e\n", N, fft_nops_real[N], err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());
	}

	tm3.reset();
	tm4.reset();
	for (int N = 2; N <= MAXFFT; N += 1) {
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

