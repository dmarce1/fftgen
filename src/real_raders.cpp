#include "fftgen.hpp"
#include <fftw3.h>

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

std::vector<int> raders_ginvq_half(int N) {
	const int g = generator(N);
	std::vector<int> ginvq;
	for (int q = 0; q < N - 1; q++) {
		int k = mod_inv(mod_pow(g, q, N), N);
		if (k > (N - 1) / 2) {
			k = N - k;
		}
		ginvq.push_back(k);
	}
	return ginvq;
}

const std::vector<int> raders_gq_half(int N) {
	const int g = generator(N);
	std::vector<int> gq;
	for (int q = 0; q < N - 1; q++) {
		int n = mod_pow(g, q, N);
		if (n > (N - 1) / 2) {
			n = N - n;
		}
		gq.push_back(n);
	}
	return gq;
}

const std::vector<int> raders_gq_sgn(int N) {
	const int g = generator(N);
	std::vector<int> gq;
	for (int q = 0; q < N - 1; q++) {
		int n = mod_pow(g, q, N);
		if (n > (N - 1) / 2) {
			gq.push_back(-1);
		} else {
			gq.push_back(1);
		}
	}
	return gq;
}

const std::vector<std::complex<double>> raders_cos_twiddle(int N) {
	const int L = (N - 1);
	std::vector<double> b(L);
	std::vector<std::complex<double>> c(L / 2 + 1);
	const auto ginvq = raders_ginvq(N);
	for (int q = 0; q < L; q++) {
		b[q] = cos(2.0 * M_PI * ginvq[q] / N);
	}
	fftw_real(c, b);
	return c;
}

const std::vector<std::complex<double>> raders_sin_twiddle(int N) {
	const int L = (N - 1);
	std::vector<double> b(L);
	std::vector<std::complex<double>> c(L / 2 + 1);
	const auto ginvq = raders_ginvq(N);
	for (int q = 0; q < L; q++) {
		b[q] = sin(-2.0 * M_PI * ginvq[q] / N);
	}
	fftw_real(c, b);
	return c;
}

void raders_fft_real(int N, int o) {
	print("{\n");
	indent();
	print("// Raders radix - %i\n", N);
	printf("// Real Raders radix - %i\n", N);
	print("double xro, xrk0;\n");
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	const auto acos = raders_cos_twiddle(N);
	const auto bsin = raders_sin_twiddle(N);
	const auto gqh = raders_gq_half(N);
	const auto gq = raders_gq(N);
	const auto ginvq = raders_ginvq(N);
	const auto ginvqh = raders_ginvq_half(N);
	const int L = (N - 1);
	print("double xok = x[%i];\n", o);
	print("double xo = x[%i];\n", o);
	print("std::array<double, %i> z;\n", L);
	for (int q = 0; q < L; q++) {
		print("xo += x[%i];\n", o + q + 1);
	}
	for (int q = 0; q < L; q++) {
		print("z[%i] = x[%i];\n", q, o + gq[q]);
	}
	print("{\n");
	indent();
	print("auto* x = z.data();\n");
	fft_bitreverse_real(L);
	fft_real(L, 0);
	deindent();
	print("}\n");
	print("z[0] *= -1.0;\n");
	for (int q = 0; q < L / 2; q += 2) {
		if (q != 0) {
			print("tmp1 = z[%i];\n", q);
		}
		if (!(q == L / 2 - 1 && (L / 2) % 2 == 1)) {
			print("tmp0 = z[%i];\n", q + 1);
		}
		if (q != 0) {
			print("z[%i] = std::fma((%.17e), z[%i], (%.17e) * z[%i]);\n", q, acos[q].real(), q, -acos[q].imag(), L - q);
		}
		if (!(q == L / 2 - 1 && (L / 2) % 2 == 1)) {
			print("z[%i] = std::fma((%.17e), z[%i], (%.17e) * z[%i]);\n", (q + 1), bsin[(q + 1)].real(), (q + 1), -bsin[(q + 1)].imag(), L - (q + 1));
		}
		if (q != 0) {
			print("z[%i] = std::fma((%.17e), tmp1, (%.17e) * z[%i]);\n", L - q, acos[q].imag(), acos[q].real(), L - q);
		}
		if (!(q == L / 2 - 1 && (L / 2) % 2 == 1)) {
			print("z[%i] = std::fma((%.17e), tmp0, (%.17e) * z[%i]);\n", L - (q + 1), bsin[(q + 1)].imag(), bsin[(q + 1)].real(), L - (q + 1));
		}
	}
	print("z[%i] *= %.17e;\n", L / 2, (L / 2) % 2 == 0 ? acos[L / 2].real() : bsin[L / 2].real());
	print("{\n");
	indent();
	print("x[%i] = z[0];\n", o);
	for (int n = 1; n < L; n++) {
		if (n <= L - n) {
			print("x[%i] = z[%i] + z[%i];\n", o + n, n, L - n);
		} else {
			print("x[%i] = z[%i] - z[%i];\n", o + n, L - n, n);
		}
	}
	if (L % 2 == 0) {
		print("x[%i] = z[%i];\n", o + L / 2, L / 2);
	}
	fft_bitreverse_real(L, fft_bitreverse_indices_real(L), o);
	fft_real(L, o);
	print("z[0] = x[%i];\n", o);
	for (int n = 1; n < L; n++) {
		if (n <= L - n) {
			print("z[%i] = x[%i] + x[%i];\n", n, o + n, o + L - n);
		} else {
			print("z[%i] = x[%i] - x[%i];\n", n, o + L - n, o + n);
		}
	}
	if (L % 2 == 0) {
		print("z[%i] = x[%i];\n", L / 2, o + L / 2);
	}
	deindent();
	print("}\n");
	for (int q = 0; q < L; q++) {
		print("z[%i] *= %.17e;\n", q, (0.5 / L));
	}
	print("x[%i] = xo;\n", o);
	for (int p = 0; p < L / 2; p++) {
		print("x[%i] = xok + z[%i] + z[%i];\n", ginvqh[p] + o, p, p + L / 2);
		if (ginvq[p] <= L / 2) {
			print("x[%i] = z[%i] - z[%i];\n", N - ginvqh[p] + o, p, p + L / 2);
		} else {
			print("x[%i] = z[%i] - z[%i];\n", N - ginvqh[p] + o, p + L / 2, p);
		}
	}
	deindent();
	print("}\n");
}

int raders_fft_real_opcnt(int N, int o) {
	const auto acos = raders_cos_twiddle(N);
	const auto bsin = raders_sin_twiddle(N);
	const auto gq = raders_gq_half(N);
	const auto ginvq = raders_ginvq_half(N);
	const int L = N - 1;
	int cnt = 0;
	for (int q = 0; q < L; q++) {
		cnt++;
	}
	cnt += fft_real_opcnt(L, 0);
	cnt += 2;
	for (int q = 0; q < L / 2; q += 2) {
		if (q != 0) {
			cnt += 4;
		}
		cnt += 4;
	}
	cnt += fft_real_opcnt(L, 0);
	for (int q = 0; q < L; q++) {
		cnt++;
	}
	for (int n = 1; n < L; n++) {
		cnt += 2;
	}
	for (int p = 0; p < L / 2; p++) {
		cnt += 3;
	}
	return cnt;
}
