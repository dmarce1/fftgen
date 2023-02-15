
#include "fftgen.hpp"
#include <list>
#include <fftw3.h>

void fft_radix_real(int N1, int N2);
int fft_radix_real_opcnt(int N1, int N2);
std::vector<int> fft_bitreverse_indices_real(int N);

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
	print("fft_real_base_%i(x);\n", L);
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
	print("fft_real_base_%i(x);\n", L);
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


std::vector<int> factors(int N) {
	std::vector<int> f;
	for (int N1 = 2; N1 <= N; N1++) {
		if (N % N1 == 0) {
			if (is_prime(N1) || N1 == 4 || N1 == 6 || N1 == 8) {
				f.push_back(N1);
			}
		}
	}
	return f;
}

std::pair<int, int> best_radix_real(int N) {
	if (is_prime(N) && N > 23) {
		return std::make_pair(-1, 0);
	} else {
		const auto fs = factors(N);
		int best_cnt = 999999999;
		int best_radix;
		for (int n = 0; n < fs.size(); n++) {
			int cnt = fft_radix_real_opcnt(fs[n], N / fs[n]);
			if (cnt < best_cnt) {
				best_cnt = cnt;
				best_radix = fs[n];
			}
		}
		return std::make_pair(best_radix, best_cnt);
	}
}

void fft_real(int N) {
	int best = best_radix_real(N).first;
	if (best < 0) {
		raders_fft_real(N, 0);
	} else {
		fft_radix_real(best, N / best);
	}
}

void fft_radix_real(int N1, int N2) {
	int N = N1 * N2;
	print("{\n");
	indent();

	if (N2 > 1) {
		for (int n1 = 0; n1 < N1; n1++) {
			print("fft_real_base_%i(x + %i);\n", N2, N2 * n1);
		}
	}
	int ipar = 0;
	for (int k2 = 0; k2 < N2 / 2 + 1; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			const int kr = N2 * n1 + k2;
			const int ki = N2 * n1 + (N2 - k2);
			const auto W = twiddle(n1 * k2, N);
			if (k2 == 0) {
			} else if (k2 == N2 / 2 && N2 % 2 == 0) {
			} else {
				if (n1 * k2 == 0) {
				} else if (n1 * k2 == N / 2 && (N % 2 == 0)) {
					print("x[%i] = -x[%i];\n", kr, kr);
					print("x[%i] = -x[%i];\n", ki, ki);
				} else if (n1 * k2 == N / 4 && (N % 4 == 0)) {
					print("std::swap(x[%i], x[%i]);\n", kr, ki);
					print("x[%i] = -x[%i];\n", ki, ki);
				} else if (n1 * k2 == 3 * N / 4 && (N % 4 == 0)) {
					print("std::swap(x[%i], x[%i]);\n", kr, ki);
					print("x[%i] = -x[%i];\n", kr, kr);
				} else {
					print("tmp%i = x[%i];\n", ipar, kr);
					print("x[%i] = std::fma(x[%i], %.17e, x[%i] * %.17e);\n", kr, kr, W.real(), ki, -W.imag());
					print("x[%i] = std::fma(x[%i], %.17e, tmp%i * %.17e);\n", ki, ki, W.real(), ipar, W.imag());
					ipar = (ipar + 1) % NPAR;
				}
			}
		}
	}
	print("{\n");
	indent();
	for (int k1 = 0; k1 < N1; k1++) {
		const int kr = N2 * k1;
		print("const auto zr%i = x[%i];\n", k1, kr);
	}
	switch (N1) {
	case 1:
		return;
	case 2:
		print("x[%i] = zr0 + zr1;\n", 0);
		print("x[%i] = zr0 - zr1;\n", N2);
		break;
	case 3:
		print("const auto tr1 = zr1 + zr2;\n");
		print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
		print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
		print("x[%i] = zr0 + tr1;\n", 0);
		print("x[%i] = tr2;\n", N2);
		print("x[%i] = -tr3;\n", N - N2);
		break;
	case 4:
		print("const auto tr1 = zr0 + zr2;\n");
		print("const auto tr2 = zr1 + zr3;\n");
		print("const auto tr3 = zr0 - zr2;\n");
		print("const auto tr4 = zr1 - zr3;\n");
		print("x[%i] = tr1 + tr2;\n", 0);
		print("x[%i] = tr3 ;\n", N2);
		print("x[%i] = - tr4;\n", N - N2);
		print("x[%i] = tr1 - tr2;\n", 2 * N2);
		print("x[%i] = tr3;\n", N - (3 * N2));
		print("x[%i] = -tr4;\n", 3 * N2);
		break;
	case 5:
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
		print("x[%i] = zr0 + tr5;\n", 0);
		print("x[%i] = tr8;\n", N2);
		print("x[%i] = -tr10;\n", N - N2);
		print("x[%i] = tr9;\n", 2 * N2);
		print("x[%i] = -tr11;\n", N - 2 * N2);
		break;
	case 6:
		print("const auto tr1 = zr2 + zr4;\n");
		print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
		print("const auto tr3 = (%24.17e) * (zr2 - zr4);\n", sin(M_PI / 3.0));
		print("const auto tr4 = zr5 + zr1;\n");
		print("const auto tr5 = std::fma(tr4, -0.5, zr3);\n");
		print("const auto tr6 = (%24.17e) * (zr5 - zr1);\n", sin(M_PI / 3.0));
		print("const auto tr7 = zr0 + tr1;\n");
		print("const auto tr8 = tr2;\n");
		print("const auto ti8 = -tr3;\n");
		print("const auto tr9 = tr2;\n");
		print("const auto ti9 = tr3;\n");
		print("const auto tr10 = zr3 + tr4;\n");
		print("const auto tr11 = tr5;\n");
		print("const auto ti11 = -tr6;\n");
		print("const auto tr12 = tr5;\n");
		print("const auto ti12 = tr6;\n");
		print("x[%i] = tr7 + tr10;\n", 0);
		print("x[%i] = tr8 - tr11;\n", N2);
		print("x[%i] = ti8 - ti11;\n", N - N2);
		print("x[%i] = tr9 + tr12;\n", 2 * N2);
		print("x[%i] = ti9 + ti12;\n", N - 2 * N2);
		print("x[%i] = tr7 - tr10;\n", 3 * N2);

		break;
	case 7: {
		constexpr double theta = 2.0 * M_PI / 7.0;
		const double c1 = cos(theta);
		const double c2 = cos(2.0 * theta);
		const double c3 = cos(3.0 * theta);
		const double c4 = sin(theta);
		const double c5 = sin(2.0 * theta);
		const double c6 = sin(3.0 * theta);
		print("const auto tr1 = zr1 + zr6;\n");
		print("const auto tr2 = zr2 + zr5;\n");
		print("const auto tr3 = zr3 + zr4;\n");
		print("const auto tr4 = zr1 - zr6;\n");
		print("const auto tr5 = zr2 - zr5;\n");
		print("const auto tr6 = zr3 - zr4;\n");
		print("const auto tr7 = zr0 - 0.5 * tr3;\n");
		print("const auto tr8 = tr1 - tr3;\n");
		print("const auto tr9 = tr2 - tr3;\n");
		print("const auto yr0 = zr0 + tr1 + tr2 + tr3;\n");
		print("const auto yr1 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c1, c2);
		print("const auto yr2 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c2, c3);
		print("const auto yr3 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c3, c1);
		print("const auto yr4 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c6, -c4, c5);
		print("const auto yr5 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c5, -c6, -c4);
		print("const auto yr6 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c4, c5, c6);
		print("x[%i] = yr0;\n", 0);
		print("x[%i] = yr1;\n", N2);
		print("x[%i] = -yr6;\n", N - N2);
		print("x[%i] = yr2;\n", 2 * N2);
		print("x[%i] = -yr5;\n", N - 2 * N2);
		print("x[%i] = yr3;\n", 3 * N2);
		print("x[%i] = -yr4;\n", N - 3 * N2);
	}
		break;
	case 8:
		printf("N = %i - real Radix 8\n", N);
		print("const auto tr1 = zr0 + zr4;\n");
		print("const auto tr2 = zr0 - zr4;\n");
		print("const auto tr3 = zr1 + zr5;\n");
		print("const auto tr4 = zr1 - zr5;\n");
		print("const auto tr5 = zr2 + zr6;\n");
		print("const auto tr6 = zr2 - zr6;\n");
		print("const auto tr7 = zr3 + zr7;\n");
		print("const auto tr8 = zr3 - zr7;\n");
		print("const auto tr9 = tr1 + tr5;\n");
		print("const auto tr10 = tr3 + tr7;\n");
		print("const auto tr11 = (%.17e) * (tr4 - tr8);\n", 1.0 / sqrt(2));
		print("const auto tr12 = (%.17e) * (tr4 + tr8);\n", 1.0 / sqrt(2));
		print("const auto yr0 = tr9 + tr10;\n");
		print("const auto yr1 = tr2 + tr11;\n");
		print("const auto yr2 = tr1 - tr5;\n");
		print("const auto yr3 = tr2 - tr11;\n");
		print("const auto yr4 = tr9 - tr10;\n");
		print("const auto yr5 = tr12 - tr6;\n");
		print("const auto yr6 = tr3 - tr7;\n");
		print("const auto yr7 = tr12 + tr6;\n");
		print("x[%i] = yr0;\n", 0);
		print("x[%i] = yr1;\n", N2);
		print("x[%i] = -yr7;\n", N - N2);
		print("x[%i] = yr2;\n", 2 * N2);
		print("x[%i] = -yr6;\n", N - 2 * N2);
		print("x[%i] = yr3;\n", 3 * N2);
		print("x[%i] = -yr5;\n", N - 3 * N2);
		print("x[%i] = yr4;\n", N - 4 * N2);
		break;

	default:
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			print("const auto txp%i = zr%i + zr%i;\n", j, j, N1 - j);
		}
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			print("const auto txm%i = zr%i - zr%i;\n", j, j, N1 - j);
		}
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			print("auto ap%i = zr0;\n", i);
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
				if (j == 1) {
					print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
				} else {
					print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
				}
			}
		}
		print("x[%i] = ", 0);
		for (int i = 0; i < N1; i++) {
			print_notab(" + zr%i", i);
		}
		print_notab(";\n");
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			print("x[%i] = ap%i;\n", i * N2, i);
			if (N - i * N2 != i * N2) {
				print("x[%i] = -bm%i;\n", N - (i * N2), i);
			}
			print("x[%i] = ap%i;\n", N - ((N1 - i) * N2), i);
			if (N - (N1 - i) * N2 != (N1 - i) * N2) {
				print("x[%i] = -bm%i;\n", (N1 - i) * N2, i);
			}
		}
	}
	deindent();
	print("}\n");

	print("for (int k2 = 1; k2 < %i; k2++) {\n", (N2 + 1) / 2);
//		print("{\n");
	indent();
	print("auto* xp = x + k2;\n");
	print("auto* xm = x + (%i - k2);\n", N);
	for (int k1 = 0; k1 < N1; k1++) {
		const int kr = N2 * k1;
		const int ki = N2 * (N1 - k1 - 1);
		const auto W = twiddle(k1 * k1, N1);
		print("const auto zr%i = xp[%i];\n", k1, kr);
		print("const auto zi%i = xm[%i];\n", k1, -ki);
	}

	switch (N1) {
	case 1:
		return;
	case 2:
		print("xp[%i] = zr0 + zr1;\n", 0);
		print("xm[%i] = zi0 + zi1;\n", -0);
		print("xm[%i] = zr0 - zr1;\n", -N2);
		print("xp[%i] = zi1 - zi0;\n", N2);
		break;
	case 3:
		print("const auto tr1 = zr1 + zr2;\n");
		print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
		print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
		print("const auto ti1 = zi1 + zi2;\n");
		print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
		print("const auto ti3 = (%24.17e) * (zi1 - zi2);\n", sqrt(3) * 0.5);
		print("xp[%i] = zr0 + tr1;\n", 0);
		print("xm[%i] = zi0 + ti1;\n", 0);
		print("xp[%i] = tr2 + ti3;\n", N2);
		print("xm[%i] = ti2 - tr3;\n", -N2);
		print("xm[%i] = tr2 - ti3;\n", -2 * N2);
		print("xp[%i] = -ti2 - tr3;\n", 2 * N2);
		break;
	case 4:
		print("const auto tr1 = zr0 + zr2;\n");
		print("const auto ti1 = zi0 + zi2;\n");
		print("const auto tr2 = zr1 + zr3;\n");
		print("const auto ti2 = zi1 + zi3;\n");
		print("const auto tr3 = zr0 - zr2;\n");
		print("const auto ti3 = zi0 - zi2;\n");
		print("const auto tr4 = zr1 - zr3;\n");
		print("const auto ti4 = zi1 - zi3;\n");
		print("xp[%i] = tr1 + tr2;\n", 0);
		print("xm[%i] = ti1 + ti2;\n", 0);
		print("xp[%i] = tr3 + ti4;\n", N2);
		print("xm[%i] = ti3 - tr4;\n", -N2);
		print("xm[%i] = tr1 - tr2;\n", -2 * N2);
		print("xp[%i] = ti2 - ti1 ;\n", 2 * N2);
		print("xm[%i] = tr3 - ti4;\n", -3 * N2);
		print("xp[%i] = -(ti3 + tr4);\n", 3 * N2);
		break;
	case 5:
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
		print("xp[%i] = zr0 + tr5;\n", 0);
		print("xm[%i] = zi0 + ti5;\n", 0);
		print("xp[%i] = tr8 + ti10;\n", N2);
		print("xm[%i] = ti8 - tr10;\n", -N2);
		print("xp[%i] = tr9 + ti11;\n", 2 * N2);
		print("xm[%i] = ti9 - tr11;\n", -2 * N2);
		print("xm[%i] = tr9 - ti11;\n", -3 * N2);
		print("xp[%i] = -ti9 - tr11;\n", 3 * N2);
		print("xm[%i] = tr8 - ti10;\n", -4 * N2);
		print("xp[%i] = -ti8 - tr10;\n", 4 * N2);
		break;
	case 6:
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
		print("xp[%i] = tr7 + tr10;\n", 0);
		print("xm[%i] = ti7 + ti10;\n", 0);
		print("xp[%i] = tr8 - tr11;\n", N2);
		print("xm[%i] = ti8 - ti11;\n", -N2);
		print("xp[%i] = tr9 + tr12;\n", 2 * N2);
		print("xm[%i] = ti9 + ti12;\n", -2 * N2);
		print("xm[%i] = tr7 - tr10;\n", -3 * N2);
		print("xp[%i] = -ti7 + ti10;\n", 3 * N2);
		print("xm[%i] = tr8 + tr11;\n", -4 * N2);
		print("xp[%i] = -ti8 - ti11;\n", 4 * N2);
		print("xm[%i] = tr9 - tr12;\n", -5 * N2);
		print("xp[%i] = -ti9 + ti12;\n", 5 * N2);
		break;
	case 7: {
		constexpr double theta = 2.0 * M_PI / 7.0;
		const double c1 = cos(theta);
		const double c2 = cos(2.0 * theta);
		const double c3 = cos(3.0 * theta);
		const double c4 = sin(theta);
		const double c5 = sin(2.0 * theta);
		const double c6 = sin(3.0 * theta);
		print("const auto tr1 = zr1 + zr6;\n");
		print("const auto ti1 = zi1 + zi6;\n");
		print("const auto tr2 = zr2 + zr5;\n");
		print("const auto ti2 = zi2 + zi5;\n");
		print("const auto tr3 = zr3 + zr4;\n");
		print("const auto ti3 = zi3 + zi4;\n");
		print("const auto tr4 = zr1 - zr6;\n");
		print("const auto ti4 = zi1 - zi6;\n");
		print("const auto tr5 = zr2 - zr5;\n");
		print("const auto ti5 = zi2 - zi5;\n");
		print("const auto tr6 = zr3 - zr4;\n");
		print("const auto ti6 = zi3 - zi4;\n");
		print("const auto tr7 = zr0 - 0.5 * tr3;\n");
		print("const auto ti7 = zi0 - 0.5 * ti3;\n");
		print("const auto tr8 = tr1 - tr3;\n");
		print("const auto ti8 = ti1 - ti3;\n");
		print("const auto tr9 = tr2 - tr3;\n");
		print("const auto ti9 = ti2 - ti3;\n");
		print("const auto yr0 = zr0 + tr1 + tr2 + tr3;\n");
		print("const auto yi0 = zi0 + ti1 + ti2 + ti3;\n");
		print("const auto yr1 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c1, c2);
		print("const auto yi1 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c1, c2);
		print("const auto yr2 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c2, c3);
		print("const auto yi2 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c2, c3);
		print("const auto yr3 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c3, c1);
		print("const auto yi3 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c3, c1);
		print("const auto yr4 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c6, -c4, c5);
		print("const auto yi4 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c6, -c4, c5);
		print("const auto yr5 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c5, -c6, -c4);
		print("const auto yi5 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c5, -c6, -c4);
		print("const auto yr6 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c4, c5, c6);
		print("const auto yi6 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c4, c5, c6);
		print("xp[%i] = yr0;\n", 0);
		print("xm[%i] = yi0;\n", 0);
		print("xp[%i] = yr1 + yi6;\n", N2);
		print("xm[%i] = yi1 - yr6;\n", -N2);
		print("xp[%i] = yr2 + yi5;\n", 2 * N2);
		print("xm[%i] = yi2 - yr5;\n", -2 * N2);
		print("xp[%i] = yr3 + yi4;\n", 3 * N2);
		print("xm[%i] = yi3 - yr4;\n", -3 * N2);
		print("xm[%i] = yr3 - yi4;\n", -4 * N2);
		print("xp[%i] = -yi3 - yr4;\n", 4 * N2);
		print("xm[%i] = yr2 - yi5;\n", -5 * N2);
		print("xp[%i] = -yi2 - yr5;\n", 5 * N2);
		print("xm[%i] = yr1 - yi6;\n", -6 * N2);
		print("xp[%i] = -yi1 - yr6;\n", 6 * N2);
	}
		break;
	case 8:
		printf("N = %i - real Radix 8\n", N);
		print("const auto tr1 = zr0 + zr4;\n");
		print("const auto ti1 = zi0 + zi4;\n");
		print("const auto tr2 = zr0 - zr4;\n");
		print("const auto ti2 = zi0 - zi4;\n");
		print("const auto tr3 = zr1 + zr5;\n");
		print("const auto ti3 = zi1 + zi5;\n");
		print("const auto tr4 = zr1 - zr5;\n");
		print("const auto ti4 = zi1 - zi5;\n");
		print("const auto tr5 = zr2 + zr6;\n");
		print("const auto ti5 = zi2 + zi6;\n");
		print("const auto tr6 = zr2 - zr6;\n");
		print("const auto ti6 = zi2 - zi6;\n");
		print("const auto tr7 = zr3 + zr7;\n");
		print("const auto ti7 = zi3 + zi7;\n");
		print("const auto tr8 = zr3 - zr7;\n");
		print("const auto ti8 = zi3 - zi7;\n");
		print("const auto tr9 = tr1 + tr5;\n");
		print("const auto ti9 = ti1 + ti5;\n");
		print("const auto tr10 = tr3 + tr7;\n");
		print("const auto ti10 = ti3 + ti7;\n");
		print("const auto tr11 = (%.17e) * (tr4 - tr8);\n", 1.0 / sqrt(2));
		print("const auto ti11 = (%.17e) * (ti4 - ti8);\n", 1.0 / sqrt(2));
		print("const auto tr12 = (%.17e) * (tr4 + tr8);\n", 1.0 / sqrt(2));
		print("const auto ti12 = (%.17e) * (ti4 + ti8);\n", 1.0 / sqrt(2));
		print("const auto yr0 = tr9 + tr10;\n");
		print("const auto yi0 = ti9 + ti10;\n");
		print("const auto yr1 = tr2 + tr11;\n");
		print("const auto yi1 = ti2 + ti11;\n");
		print("const auto yr2 = tr1 - tr5;\n");
		print("const auto yi2 = ti1 - ti5;\n");
		print("const auto yr3 = tr2 - tr11;\n");
		print("const auto yi3 = ti2 - ti11;\n");
		print("const auto yr4 = tr9 - tr10;\n");
		print("const auto yi4 = ti9 - ti10;\n");
		print("const auto yr5 = tr12 - tr6;\n");
		print("const auto yi5 = ti12 - ti6;\n");
		print("const auto yr6 = tr3 - tr7;\n");
		print("const auto yi6 = ti3 - ti7;\n");
		print("const auto yr7 = tr12 + tr6;\n");
		print("const auto yi7 = ti12 + ti6;\n");
		print("xp[%i] = yr0;\n", 0);
		print("xm[%i] = yi0;\n", 0);
		print("xp[%i] = yr1 + yi7;\n", N2);
		print("xm[%i] = yi1 - yr7;\n", -N2);
		print("xp[%i] = yr2 + yi6;\n", 2 * N2);
		print("xm[%i] = yi2 - yr6;\n", -2 * N2);
		print("xp[%i] = yr3 + yi5;\n", 3 * N2);
		print("xm[%i] = yi3 - yr5;\n", -3 * N2);
		print("xm[%i] = yr4;\n", -4 * N2);
		print("xp[%i] = -yi4;\n", 4 * N2);
		print("xm[%i] = yr3 - yi5;\n", -5 * N2);
		print("xp[%i] = -(yi3 + yr5);\n", 5 * N2);
		print("xm[%i] = yr2 - yi6;\n", -6 * N2);
		print("xp[%i] = -(yi2 + yr6);\n", 6 * N2);
		print("xm[%i] = yr1 - yi7;\n", -7 * N2);
		print("xp[%i] = -(yi1 + yr7);\n", 7 * N2);
		break;
	default:
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			print("const auto txp%i = zr%i + zr%i;\n", j, j, N1 - j);
		}
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			print("const auto txm%i = zr%i - zr%i;\n", j, j, N1 - j);
		}
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			print("const auto typ%i = zi%i + zi%i;\n", j, j, N1 - j);
		}
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			print("const auto tym%i = zi%i - zi%i;\n", j, j, N1 - j);
		}
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			print("auto ap%i = zr0;\n", i);
			print("auto bp%i = zi0;\n", i);
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
				print("bp%i = std::fma(typ%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
				if (j == 1) {
					print("double am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
					print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
				} else {
					print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
					print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
				}
			}
		}
		print("xp[%i] = ", 0);
		for (int i = 0; i < N1; i++) {
			print_notab(" + zr%i", i);
		}
		print_notab(";\n");
		print("xm[%i] = ", 0);
		for (int i = 0; i < N1; i++) {
			print_notab(" + zi%i", i);
		}
		print_notab(";\n");
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			print("xp[%i] = ap%i + am%i;\n", i * N2, i, i);
			print("xm[%i] = bp%i - bm%i;\n", -i * N2, i, i);
			print("xm[%i] = ap%i - am%i;\n", -(N1 - i) * N2, i, i);
			print("xp[%i] = -bp%i - bm%i;\n", (N1 - i) * N2, i, i);
		}
	}
	deindent();
	print("}\n");

//	}
	if (N2 % 2 == 0) {

		int k2 = N2 / 2;
		print("{\n");
		indent();

		switch (N1) {
		case 1:
			return;
		case 2:
			for (int k1 = 0; k1 < N1; k1++) {
				const int kr = N2 * k1 + k2;
				print("const auto x%i = x[%i];\n", k1, kr);
			}
			print("x[%i] = x0;\n", k2);
			print("x[%i] = -x1;\n", N - k2);
			break;
		case 3:
			for (int k1 = 0; k1 < N1; k1++) {
				const int kr = N2 * k1 + k2;
				print("const auto x%i = x[%i];\n", k1, kr);
			}
			print("const auto t1 = x1 - x2;\n");
			print("x[%i] = x0 + 0.5 * t1;\n", k2);
			print("x[%i] = x0 - t1;\n", N2 + k2);
			print("x[%i] = (%.17e) * (x1 + x2);\n", N - k2, -sin(M_PI / 3.0));
			break;
		case 4:
			for (int k1 = 0; k1 < N1; k1++) {
				const int kr = N2 * k1 + k2;
				print("const auto x%i = x[%i];\n", k1, kr);
			}
			print("const auto t1 = (x1 - x3) * (%.17e);\n", 1.0 / sqrt(2));
			print("const auto t2 = (x1 + x3) * (%.17e);\n", 1.0 / sqrt(2));
			print("x[%i] = x0 + t1;\n", k2);
			print("x[%i] = -x2 - t2;\n", N - k2);
			print("x[%i] = x0 - t1;\n", N2 + k2);
			print("x[%i] = x2 - t2;\n", N - (N2 + k2));
			break;
		case 5:
			for (int k1 = 0; k1 < N1; k1++) {
				const int kr = N2 * k1 + k2;
				print("const auto x%i = x[%i];\n", k1, kr);
			}
			print("const auto t1 = x1 - x4;\n");
			print("const auto t2 = x1 + x4;\n");
			print("const auto t3 = x2 - x3;\n");
			print("const auto t4 = x2 + x3;\n");
			print("const auto t5 = t1 - t3;\n");
			print("const auto t6 = std::fma(0.25, t5, x0);\n");
			print("const auto t7 = (%.17e) * (t1 + t3);\n", sqrt(5) * 0.25);
			print("x[%i] = t6 + t7;\n", k2);
			print("x[%i] = std::fma((%.17e), t2, (%.17e) * t4);\n", N - k2, -sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
			print("x[%i] = t6 - t7;\n", N2 + k2);
			print("x[%i] = std::fma((%.17e), t2, (%.17e) * t4);\n", N - (N2 + k2), -sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
			print("x[%i] = x0 - t5;\n", 2 * N2 + k2);
			break;
		case 6:
			for (int k1 = 0; k1 < N1; k1++) {
				const int kr = N2 * k1 + k2;
				print("const auto x%i = x[%i];\n", k1, kr);
			}
			print("const auto t1 = (%.17e) * (x5 - x1);\n", sin(M_PI / 3.0));
			print("const auto t2 = (%.17e) * (x2 + x4);\n", sin(M_PI / 3.0));
			print("const auto t3 = x2 - x4;\n");
			print("const auto t4 = x1 + x5;\n");
			print("const auto t5 = std::fma(0.5, t3, x0);\n");
			print("const auto t6 = -std::fma(0.5, t4, x3);\n");
			print("x[%i] = t5 - t1;\n", k2);
			print("x[%i] = t6 - t2;\n", N - k2);
			print("x[%i] = x0 - t3;\n", N2 + k2);
			print("x[%i] = x3 - t4;\n", N - (N2 + k2));
			print("x[%i] = t5 + t1;\n", 2 * N2 + k2);
			print("x[%i] = t6 + t2;\n", N - (2 * N2 + k2));
			break;
		case 7: {
			for (int k1 = 0; k1 < N1; k1++) {
				const int kr = N2 * k1 + k2;
				print("const auto x%i = x[%i];\n", k1, kr);
			}
			constexpr double theta = 2.0 * M_PI / 7.0;
			const double c1 = cos(theta);
			const double c2 = cos(2.0 * theta);
			const double c3 = cos(3.0 * theta);
			const double c4 = sin(theta);
			const double c5 = sin(2.0 * theta);
			const double c6 = sin(3.0 * theta);
			print("const auto zr0 = x0;\n");
			print("const auto zr1 = x1 * (%.17e);\n", cos(1.0 * M_PI / 7.0));
			print("const auto zi1 = x1 * (%.17e);\n", -sin(1.0 * M_PI / 7.0));
			print("const auto zr2 = x2 * (%.17e);\n", cos(2.0 * M_PI / 7.0));
			print("const auto zi2 = x2 * (%.17e);\n", -sin(2.0 * M_PI / 7.0));
			print("const auto zr3 = x3 * (%.17e);\n", cos(3.0 * M_PI / 7.0));
			print("const auto zi3 = x3 * (%.17e);\n", -sin(3.0 * M_PI / 7.0));
			print("const auto zr4 = x4 * (%.17e);\n", cos(4.0 * M_PI / 7.0));
			print("const auto zi4 = x4 * (%.17e);\n", -sin(4.0 * M_PI / 7.0));
			print("const auto zr5 = x5 * (%.17e);\n", cos(5.0 * M_PI / 7.0));
			print("const auto zi5 = x5 * (%.17e);\n", -sin(5.0 * M_PI / 7.0));
			print("const auto zr6 = x6 * (%.17e);\n", cos(6.0 * M_PI / 7.0));
			print("const auto zi6 = x6 * (%.17e);\n", -sin(6.0 * M_PI / 7.0));
			print("const auto tr1 = zr1 + zr6;\n");
			print("const auto ti1 = zi1 + zi6;\n");
			print("const auto tr2 = zr2 + zr5;\n");
			print("const auto ti2 = zi2 + zi5;\n");
			print("const auto tr3 = zr3 + zr4;\n");
			print("const auto ti3 = zi3 + zi4;\n");
			print("const auto tr4 = zr1 - zr6;\n");
			print("const auto ti4 = zi1 - zi6;\n");
			print("const auto tr5 = zr2 - zr5;\n");
			print("const auto ti5 = zi2 - zi5;\n");
			print("const auto tr6 = zr3 - zr4;\n");
			print("const auto ti6 = zi3 - zi4;\n");
			print("const auto tr7 = zr0 - 0.5 * tr3;\n");
			print("const auto ti7 = -0.5 * ti3;\n");
			print("const auto tr8 = tr1 - tr3;\n");
			print("const auto ti8 = ti1 - ti3;\n");
			print("const auto tr9 = tr2 - tr3;\n");
			print("const auto ti9 = ti2 - ti3;\n");
			print("const auto yr0 = zr0 + tr1 + tr2 + tr3;\n");
			print("const auto yi0 = ti1 + ti2 + ti3;\n");
			print("const auto yr1 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c1, c2);
			print("const auto yi1 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c1, c2);
			print("const auto yr2 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c2, c3);
			print("const auto yi2 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c2, c3);
			print("const auto yr3 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c3, c1);
			print("const auto yi3 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c3, c1);
			print("const auto yr4 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c6, -c4, c5);
			print("const auto yi4 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c6, -c4, c5);
			print("const auto yr5 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c5, -c6, -c4);
			print("const auto yi5 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c5, -c6, -c4);
			print("const auto yr6 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c4, c5, c6);
			print("const auto yi6 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c4, c5, c6);
			print("x[%i] = yr0;\n", k2);
			print("x[%i] = yi0;\n", N - k2);
			print("x[%i] = yr1 + yi6;\n", N2 + k2);
			print("x[%i] = yi1 - yr6;\n", N - (N2 + k2));
			print("x[%i] = yr2 + yi5;\n", 2 * N2 + k2);
			print("x[%i] = yi2 - yr5;\n", N - (2 * N2 + k2));
			print("x[%i] = yr3 + yi4;\n", 3 * N2 + k2);
		}
			break;
		case 8:
			printf("N = %i Real - half complex\n", N);
			for (int k1 = 0; k1 < N1; k1++) {
				const int kr = N2 * k1 + k2;
				print("const auto x%i = x[%i];\n", k1, kr);
			}
			print("const auto zr0 = x0;\n");
			print("const auto zr1 = x1 * (%.17e);\n", cos(1.0 * M_PI / 8.0));
			print("const auto zi1 = x1 * (%.17e);\n", -sin(1.0 * M_PI / 8.0));
			print("const auto zr2 = x2 * (%.17e);\n", cos(2.0 * M_PI / 8.0));
			print("const auto zr3 = x3 * (%.17e);\n", cos(3.0 * M_PI / 8.0));
			print("const auto zi3 = x3 * (%.17e);\n", -sin(3.0 * M_PI / 8.0));
			print("const auto zi4 = -x4;\n");
			print("const auto zr5 = x5 * (%.17e);\n", cos(5.0 * M_PI / 8.0));
			print("const auto zi5 = x5 * (%.17e);\n", -sin(5.0 * M_PI / 8.0));
			print("const auto zr6 = x6 * (%.17e);\n", cos(6.0 * M_PI / 8.0));
			print("const auto zr7 = x7 * (%.17e);\n", cos(7.0 * M_PI / 8.0));
			print("const auto zi7 = x7 * (%.17e);\n", -sin(7.0 * M_PI / 8.0));
			print("const auto tr1 = zr0;\n");
			print("const auto ti1 = zi4;\n");
			print("const auto tr2 = zr0;\n");
			print("const auto ti2 = -zi4;\n");
			print("const auto tr3 = zr1 + zr5;\n");
			print("const auto ti3 = zi1 + zi5;\n");
			print("const auto tr4 = zr1 - zr5;\n");
			print("const auto ti4 = zi1 - zi5;\n");
			print("const auto tr5 = zr2 + zr6;\n");
			print("const auto tr6 = zr2 - zr6;\n");
			print("const auto tr7 = zr3 + zr7;\n");
			print("const auto ti7 = zi3 + zi7;\n");
			print("const auto tr8 = zr3 - zr7;\n");
			print("const auto ti8 = zi3 - zi7;\n");
			print("const auto tr9 = tr1 + tr5;\n");
			print("const auto ti9 = ti1 - tr6;\n");
			print("const auto tr10 = tr3 + tr7;\n");
			print("const auto ti10 = ti3 + ti7;\n");
			print("const auto tr11 = (%.17e) * (tr4 - tr8);\n", 1.0 / sqrt(2));
			print("const auto ti11 = (%.17e) * (ti4 - ti8);\n", 1.0 / sqrt(2));
			print("const auto tr12 = (%.17e) * (tr4 + tr8);\n", 1.0 / sqrt(2));
			print("const auto ti12 = (%.17e) * (ti4 + ti8);\n", 1.0 / sqrt(2));
			print("const auto yr0 = tr9 + tr10;\n");
			print("const auto yi0 = ti9 + ti10;\n");
			print("const auto yr1 = tr2 + tr11;\n");
			print("const auto yi1 = ti2 + ti11;\n");
			print("const auto yr2 = tr1 - tr5;\n");
			print("const auto yi2 = ti1 + tr6;\n");
			print("const auto yr3 = tr2 - tr11;\n");
			print("const auto yi3 = ti2 - ti11;\n");
//			print("const auto yr4 = tr9 - tr10;\n");
//			print("const auto yi4 = ti9 - ti10;\n");
			print("const auto yr5 = tr12 - tr6;\n");
			print("const auto yi5 = ti12 + tr5;\n");
			print("const auto yr6 = tr3 - tr7;\n");
			print("const auto yi6 = ti3 - ti7;\n");
			print("const auto yr7 = tr12 + tr6;\n");
			print("const auto yi7 = ti12 - tr5;\n");
			print("x[%i] = yr0;\n", k2);
			print("x[%i] = yi0;\n", N - k2);
			print("x[%i] = yr1 + yi7;\n", N2 + k2);
			print("x[%i] = yi1 - yr7;\n", N - (N2 + k2));
			print("x[%i] = yr2 + yi6;\n", 2 * N2 + k2);
			print("x[%i] = yi2 - yr6;\n", N - (2 * N2 + k2));
			print("x[%i] = yr3 + yi5;\n", 3 * N2 + k2);
			print("x[%i] = yi3 - yr5;\n", N - (3 * N2 + k2));
			break;
		default:
			print("const auto zr0 = x[%i];\n", k2);
			for (int n1 = 1; n1 < N1; n1++) {
				const int kr = N2 * n1 + k2;
				const auto W = twiddle(n1 * k2, N);
				print("const auto zr%i = x[%i] * (%.17e);\n", n1, kr, W.real());
				print("const auto zi%i = x[%i] * (%.17e);\n", n1, kr, W.imag());
			}
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("const auto txp%i = zr%i + zr%i;\n", j, j, N1 - j);
			}
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("const auto txm%i = zr%i - zr%i;\n", j, j, N1 - j);
			}
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("const auto typ%i = zi%i + zi%i;\n", j, j, N1 - j);
			}
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("const auto tym%i = zi%i - zi%i;\n", j, j, N1 - j);
			}
			for (int i = 1; i <= (N1 - 1) / 2; i++) {
				print("auto ap%i = zr0;\n", i);
				for (int j = 1; j <= (N1 - 1) / 2; j++) {
					print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
					if (j == 1) {
						print("double bp%i = typ%i * (%24.17e);\n", i, j, cos(2.0 * M_PI * j * i / N1));
						print("double am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
						print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
					} else {
						print("bp%i = std::fma(typ%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
						print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
						print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
					}
				}
			}
			print("x[%i] = ", k2);
			for (int i = 0; i < N1; i++) {
				print_notab(" + zr%i", i);
			}
			print_notab(";\n");
			print("x[%i] = ", N - k2);
			for (int i = 1; i < N1; i++) {
				print_notab(" + zi%i", i);
			}
			print_notab(";\n");
			for (int i = 1; i <= (N1 - 1) / 2; i++) {
				if (i == (N1 - 1) / 2) {
					print("x[%i] = ap%i + am%i;\n", k2 + i * N2, i, i);
				} else {
					print("x[%i] = ap%i + am%i;\n", k2 + i * N2, i, i);
					print("x[%i] = bp%i - bm%i;\n", N - (k2 + i * N2), i, i);
				}
			}
		}
		deindent();
		print("}\n");

	}

	deindent();
	print("}\n");

}

int fft_real_opcnt(int N) {
	int best = best_radix_real(N).first;
	return fft_radix_real_opcnt(best, N / best);
}

int fft_radix_real_opcnt(int N1, int N2) {
	int N = N1 * N2;
	int cnt = 0;
	if (N2 > 1) {
		for (int n1 = 0; n1 < N1; n1++) {
			cnt += fft_real_opcnt(N2);
		}
	}
	for (int k2 = 0; k2 < N2 / 2 + 1; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			const int kr = N2 * n1 + k2;
			const int ki = N2 * n1 + (N2 - k2);
			const auto W = twiddle(n1 * k2, N);
			if (k2 == 0) {
			} else if (k2 == N2 / 2 && N2 % 2 == 0) {
			} else {
				if (n1 * k2 == 0) {
				} else if (n1 * k2 == N / 2 && (N % 2 == 0)) {
					cnt += 2;
				} else if (n1 * k2 == N / 4 && (N % 4 == 0)) {
					cnt++;
				} else if (n1 * k2 == 3 * N / 4 && (N % 4 == 0)) {
					cnt++;
				} else {
					cnt += 4;
				}
			}
		}
	}
	switch (N1) {
	case 2:
		cnt += 2;
		break;
	case 3:
		cnt += 6;
		break;
	case 4:
		cnt += 8;
		break;
	case 5:
		cnt += 17;
		break;
	case 6:
		cnt += 17;
		break;
	case 8:
		cnt += 1;
		break;

	default:
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			cnt += 2;
		}
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				cnt += 2;
			}
		}
		for (int i = 1; i < N1; i++) {
			cnt++;
		}
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			if (N - i * N2 != i * N2) {
				cnt++;
			}
			if (N - (N1 - i) * N2 != (N1 - i) * N2) {
				cnt++;
			}
		}
	}

	switch (N1) {
	case 2:
		cnt += 4;
		break;
	case 3:
		cnt += 15;
		break;
	case 4:
		cnt += 16;
		break;
	case 5:
		cnt += 40;
		break;
	case 6:
		cnt += 43;
		break;
	case 8:
		cnt += 1;
		break;
	default:
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			cnt += 4;
		}
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				cnt += 4;
			}
		}
		for (int i = 1; i < N1; i++) {
			cnt += 2;
		}
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			cnt += 5;
		}
	}
	if (N2 % 2 == 0) {

		int k2 = N2 / 2;

		switch (N1) {
		case 2:
			cnt++;
			break;
		case 3:
			cnt += 5;
			break;
		case 4:
			cnt += 9;
			break;
		case 5:
			cnt += 15;
			break;
		case 6:
			cnt += 15;
			break;
		case 8:
			cnt += 1;
			break;

		default:
			for (int n1 = 1; n1 < N1; n1++) {
				cnt += 2;
			}
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				cnt += 4;
			}
			for (int i = 1; i <= (N1 - 1) / 2; i++) {
				for (int j = 1; j <= (N1 - 1) / 2; j++) {
					cnt += 4;
				}
			}
			for (int i = 1; i < N1; i++) {
				cnt++;
			}
			for (int i = 2; i < N1; i++) {
				cnt++;
			}
			for (int i = 1; i <= (N1 - 1) / 2; i++) {
				if (i == (N1 - 1) / 2) {
					cnt++;
				} else {
					cnt++;
					cnt++;
				}
			}
		}
	}
	return cnt;

}

std::vector<int> fft_bitr_real(int N, int o, std::vector<int> I);

std::vector<int> fft_bitreverse_indices_real(int N) {
	std::vector<int> I;
	std::vector<int> indices;
	I.resize(0);
	for (int n = 0; n < N; n++) {
		I.push_back(n);
	}
	indices = fft_bitr_real(N, 0, I);
	return indices;
}

std::vector<int> fft_bitr_real(int N, int o, std::vector<int> I) {
	if (N == 1) {
		return I;
	}
	std::vector<int> L;
	std::vector<int> J;
	std::vector<int> K;
	int r = best_radix_real(N).first;
	if (r < 0) {
		return I;
	} else {
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr_real(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
		return L;
	}
}

void fft_bitreverse_real(int N, std::vector<int> indices, int o) {
	if (indices.size() == 0) {
		indices = fft_bitreverse_indices_real(N);
//		indices = reverse_indices(indices);
	}
	std::vector<std::list<int>> strings;
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
				string.push_back(index_real(o, current, 0, N));
				while (first != next) {
					string.push_back(index_real(o, next, 0, N));
					current = next;
					touched[current] = true;
					next = indices[next];
				}
			}
			strings.push_back(string);
		}
	}
	int n = 0;
	std::vector<std::vector<std::string>> cmds(NPAR);
	while (n < strings.size()) {
		int smallest = 9999999999, smalli;
		for (int i = 0; i < NPAR; i++) {
			if (cmds[i].size() < smallest) {
				smallest = cmds[i].size();
				smalli = i;
			}
		}
		std::vector<std::string> cmd;
		char* buf;
		while (strings[n].size()) {
			if (strings[n].front() == -1) {
				strings[n].pop_front();
				if (strings[n].size()) {
					asprintf(&buf, "tmp%i = x[%i];\n", smalli, strings[n].front());
					cmd.push_back(buf);
					free(buf);
				}
			} else {
				int to = strings[n].front();
				strings[n].pop_front();
				if (strings[n].size()) {
					asprintf(&buf, "x[%i] = x[%i];\n", to, strings[n].front());
					cmd.push_back(buf);
					free(buf);
				} else {
					asprintf(&buf, "x[%i] = tmp%i;\n", to, smalli);
					cmd.push_back(buf);
					free(buf);
				}
			}
		}
		cmds[smalli].insert(cmds[smalli].end(), cmd.begin(), cmd.end());
		n++;
	}
	n = 0;
	bool done = false;
	while (!done) {
		done = true;
		for (int j = 0; j < NPAR; j++) {
			if (cmds[j].size() > n) {
				print("%s", cmds[j][n].c_str());
				done = false;
			}
		}
		n++;
	}
}
