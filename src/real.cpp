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
	print("double xok = x[%i];\n", o + L);
	print("double xo = x[%i];\n", o + L);
	print("for (int q = 0; q < %i; q++) {\n", L);
	indent();
	print("xo += x[q];\n");
	deindent();
	print("}\n");
	print("fft_real_base_%i(x);\n", L);
	if (N > LOOP_N) {
		print("static const raders_cos_twiddles<%i> acos;\n", N);
		print("static const raders_sin_twiddles<%i> bsin;\n", N);
		print("const auto b = bsin[0];\n", L / 2);
		print("tmp1 = x[1];\n");
		print("x[1] = b.real() * x[1] - b.imag() * x[%i];\n", L - 1);
		print("x[%i] = std::fma(b.imag(), tmp1, b.real() * x[%i]);\n", L - 1, L - 1);

		print("for (int q = 2; q < %i; q += 2) {\n", L / 2);
		indent();
		print("const auto a = acos[q >> 1];\n");
		print("tmp1 = x[q];\n");
		print("x[q] = a.real() * x[q] - a.imag() * x[%i - q];\n", L);
		print("x[%i - q] = std::fma(a.imag(), tmp1, a.real() * x[%i - q]);\n", L, L);
		deindent();
		print("}\n");

		print("for (int q = 3; q < %i; q += 2) {\n", (L - 1) / 2 + 1);
		indent();
		print("const auto b = bsin[q >> 1];\n");
		print("tmp1 = x[q];\n");
		print("x[q] = b.real() * x[q] - b.imag() * x[%i - q];\n", L);
		print("x[%i - q] = std::fma(b.imag(), tmp1, b.real() * x[%i - q]);\n", L, L);
		deindent();
		print("}\n");
	} else {
		for (int q = 0; q < L / 2; q += 2) {
			if (q != 0) {
				print("tmp1 = x[%i];\n", q);
			}
			if (!(q == L / 2 - 1 && (L / 2) % 2 == 1)) {
				print("tmp0 = x[%i];\n", q + 1);
			}
			if (q != 0) {
				print("x[%i] = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", q, acos[q].real(), q, -acos[q].imag(), L - q);
			}
			if (!(q == L / 2 - 1 && (L / 2) % 2 == 1)) {
				print("x[%i] = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", (q + 1), bsin[(q + 1)].real(), (q + 1), -bsin[(q + 1)].imag(), L - (q + 1));
			}
			if (q != 0) {
				print("x[%i] = std::fma((%.17e), tmp1, (%.17e) * x[%i]);\n", L - q, acos[q].imag(), acos[q].real(), L - q);
			}
			if (!(q == L / 2 - 1 && (L / 2) % 2 == 1)) {
				print("x[%i] = std::fma((%.17e), tmp0, (%.17e) * x[%i]);\n", L - (q + 1), bsin[(q + 1)].imag(), bsin[(q + 1)].real(), L - (q + 1));
			}
		}
	}
	print("x[0] *= -1.0;\n");
	print("x[%i] *= %.17e;\n", L / 2, (L / 2) % 2 == 0 ? acos[L / 2].real() : bsin[L / 2].real());
	if (N >= LOOP_N) {
		print("for (int n = 1; n < %i - n; n++) {\n", L);
		indent();
		print("tmp0 = x[n];\n");
		print("const int Lmn = %i - n;\n", L);
		print("x[n] += x[Lmn];\n", L);
		print("x[Lmn] = tmp0 - x[Lmn];\n", L, L);
		deindent();
		print("}\n");
	} else {
		for (int n = 1; n < L - n; n++) {
			print("tmp0 = x[%i];\n", n);
			print("x[%i] += x[%i];\n", n, L - n);
			print("x[%i] = tmp0 - x[%i];\n", L - n, L - n);
		}
	}
	print("fft_bitreverse_real_%i(x);\n", L);
	print("fft_real_base_%i(x);\n", L);
	if (N >= LOOP_N) {
		print("for (int n = 1; n < %i - n; n++) {\n", L);
		indent();
		print("tmp0 = x[n];\n");
		print("const int Lmn = %i - n;\n", L);
		print("x[n] += x[Lmn];\n", L);
		print("x[Lmn] = tmp0 - x[Lmn];\n", L, L);
		deindent();
		print("}\n");
	} else {
		for (int n = 1; n < L - n; n++) {
			print("tmp0 = x[%i];\n", n);
			print("x[%i] += x[%i];\n", n, L - n);
			print("x[%i] = tmp0 - x[%i];\n", L - n, L - n);
		}
	}
	print("for (int q = 0; q < %i; q++) {\n", L);
	indent();
	print("x[q] *= %.17e;\n", (0.5 / L));
	deindent();
	print("}\n");
	if (N >= LOOP_N) {
		print("for (int p = 0; p < %i; p++) {\n", L / 2);
		indent();
		print("tmp0 = x[p];\n");
		print("const int q = p + %i;\n", L / 2);
		print("x[p] += xok + x[q];\n");
		print("x[q] = tmp0 - x[q];\n");
		deindent();
		print("}\n");
		for (int p = 0; p < L / 2; p++) {
			if (ginvq[p] <= L / 2) {
			} else {
				print("x[%i] = -x[%i];\n", p + L / 2 + o, p + L / 2);
			}
		}
	} else {
		for (int p = 0; p < L / 2; p++) {
			print("tmp0 = x[%i];\n", p);
			print("x[%i] += xok + x[%i];\n", p + o, p + L / 2);
			if (ginvq[p] <= L / 2) {
				print("x[%i] = tmp0 - x[%i];\n", p + L / 2 + o, p + L / 2);
			} else {
				print("x[%i] = x[%i] - tmp0;\n", p + L / 2 + o, p + L / 2);
			}
		}
	}
	std::vector<int> I(N);
	for (int p = 0; p < L / 2; p++) {
		I[ginvqh[p]] = p;
		I[N - ginvqh[p]] = p + L / 2;
	}
	I[0] = N - 1;
	fft_bitreverse_real(N, I);
	print("x[%i] = xo;\n", o);
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
	const int rader_radix = is_prime(N) ? RADER_LEN : 1000000000;
	if (N % 7 == 0) {
		return std::make_pair(7, 0);
	}
	if (N % 9 == 0) {
		return std::make_pair(9, 0);
	}
	if (N % 6 == 0) {
		return std::make_pair(6, 0);
	}
	int twopow = 0;
	int n = N;
	while (n % 2 == 0) {
		twopow++;
		n /= 2;
	}
	if (twopow >= 3 && twopow % 3 == 0) {
		return std::make_pair(8, 0);
	}
	if (twopow >= 2 && twopow % 2 == 0) {
		return std::make_pair(4, 0);
	}
	for (int r = 2; r <= rader_radix; r++) {
		if (is_prime(r) && N % r == 0) {
			return std::make_pair(r, 0);
		}
	}
	return std::make_pair(-1, 0);
	/*


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
	 }*/
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
	if (N >= LOOP_N) {
		if (N2 > 1) {
			print("static twiddle_set<%i> twiddles;\n", N);
			print("for( int n1 = 1; n1 < %i; n1++ ) {\n", N1);
			indent();
			print("const int n1N2 = n1 * %i;\n", N2);
			print("for( int k2 = 1; k2 < %i; k2++ ) {\n", (N2 + 1) / 2);
			indent();
			print("const int i = n1N2 + k2;\n");
			print("const auto tw = twiddles[n1 * k2];\n");
			print("const int ir = %i * n1 + k2;\n", N2);
			print("const int ii = %i * n1 + (%i - k2);\n", N2, N2);
			print("tmp0 = x[ir];\n");
			print("x[ir] = x[ir] * tw.real() - x[ii] * tw.imag();\n");
			print("x[ii] = std::fma(tmp0, tw.imag(), x[ii] * tw.real());\n");
			deindent();
			print("}\n");
			deindent();
			print("}\n");
		}
	} else {
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
	}
	print("{\n");
	indent();
	std::vector<std::string> out;
	std::vector<std::string> in;
	for (int i = 0; i < N1; i++) {
		out.push_back(std::string("x[" + std::to_string((i * N2)) + "]"));
	}
	in = out;
	print_real_short_fft(N1, in, out);

	deindent();
	print("}\n");
	if ((N2 + 1) / 2 > 1) {
		if ((N2 + 1) / 2 == 2) {
			print("const int k2 = 1;\n");
			print("{\n");
		} else {
			print("for (int k2 = 1; k2 < %i; k2++) {\n", (N2 + 1) / 2);
		}
		indent();
		print("auto* xp = x + k2;\n");
		print("auto* xm = x + (%i - k2);\n", N);
		out.resize(0);
		in.resize(0);
		for (int k1 = 0; k1 < N1; k1++) {
			const int kr = N2 * k1;
			const int ki = N2 * (N1 - k1 - 1);
			const auto W = twiddle(k1 * k1, N1);
			print("const auto zr%i = xp[%i];\n", k1, kr);
			print("const auto zi%i = xm[%i];\n", k1, -ki);
			in.push_back(std::string("zr") + std::to_string(k1));
			in.push_back(std::string("zi") + std::to_string(k1));
			if (k1 < N1 - k1) {
				out.push_back(std::string("xp[") + std::to_string(kr) + "]");
				out.push_back(std::string("xm[") + std::to_string(-kr) + "]");
			} else {
				out.push_back(std::string("xm[") + std::to_string(-kr) + "]");
				out.push_back(std::string("xp[") + std::to_string(kr) + "]");
			}
		}
		print_complex_short_fft(N1, in, out);
		for (int k1 = 0; k1 < N1; k1++) {
			if (k1 >= N1 - k1) {
				const int kr = N2 * k1;
				print("xp[%i] = -xp[%i];\n", kr, kr);
			}
		}
		deindent();
		print("}\n");

	}
	if (N2 % 2 == 0) {

		int k2 = N2 / 2;
		print("{\n");
		indent();
		std::vector<std::string> out;
		std::vector<std::string> in;
		for (int i = 0; i < N1; i++) {
			out.push_back(std::string("x[" + std::to_string((i * N2 + N2 / 2)) + "]"));
		}
		in = out;
		print_skew_short_fft(N1, in, out);
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
		std::vector<int> I3(N);
		std::vector<int> I1(N);
		const auto gq = raders_gq(N);
		auto I2 = fft_bitreverse_indices_real(N - 1);
		for (int i = 0; i < I2.size(); i++) {
			I3[i] = I[gq[I2[i]]];
		}
		I3[N - 1] = I[0];
		return I3;
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
