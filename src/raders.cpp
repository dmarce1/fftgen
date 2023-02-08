#include "fftgen.hpp"



std::vector<int> raders_ginvq(int N) {
	const int g = generator(N);
	std::vector<int> ginvq;
	for (int q = 0; q < N - 1; q++) {
		ginvq.push_back(mod_inv(mod_pow(g, q, N), N));
	}
	return ginvq;
}


const std::vector<int> raders_gq(int N) {
	const int g = generator(N);
	std::vector<int> gq;
	for (int q = 0; q < N - 1; q++) {
		gq.push_back(mod_pow(g, q, N));
	}
	return gq;
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
			fft(N2, o + n1 * N2, true);
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
				print("tmp0 = x[%i];\n", index(o + N1 * k2, n1, 0, N1));
				print("x[%i] = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", index(o + N1 * k2, n1, 0, N1), W.real(), index(o + N1 * k2, n1, 0, N1), -W.imag(), index(o + N1 * k2, n1, 1, N1));
				print("x[%i] = std::fma((%.17e), tmp0, (%.17e) * x[%i]);\n", index(o + N1 * k2, n1, 1, N1), W.imag(), W.real(), index(o + N1 * k2, n1, 1, N1));
			}
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		print("xro[%i] = x[%i];\n", k2, index(o + k2 * N1, 0, 0, N2));
		print("xrk0[%i] = x[%i] ", k2, index(o + k2 * N1, 0, 0, N2));
		for (int k1 = 1; k1 < N1; k1++) {
			print_notab(" + x[%i]", index(o + N1 * k2, k1, 0, N2));
		}
		print_notab(";\n");
		print("xio[%i] = x[%i];\n", k2, index(o + k2 * N1, 0, 1, N2));
		print("xik0[%i] = x[%i] ", k2, index(o + k2 * N1, 0, 1, N2));
		for (int k1 = 1; k1 < N1; k1++) {
			print_notab(" + x[%i]", index(o + N1 * k2, k1, 1, N2));
		}
		print_notab(";\n");
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
		fft(N1 - 1, o + k2 * N1, true);
		deindent();
		print("}\n");
	}
	for (int k2 = 0; k2 < N2; k2++) {
		print("x[%i] = -x[%i];\n", index(o + k2 * N1, 0, 0, N1 - 1), index(o + k2 * N1, 0, 0, N1 - 1));
		for (int q = 1; q < N1 - 1; q++) {
			print("tmp0 = x[%i];\n", index(o + N1 * k2, q, 0, N1 - 1));
			print("x[%i] = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", index(o + N1 * k2, q, 0, N1 - 1), b[q].real(), index(o + N1 * k2, q, 0, N1 - 1), -b[q].imag(), index(o + N1 * k2, q, 1, N1 - 1));
			print("x[%i] = std::fma((%.17e), tmp0, (%.17e) * x[%i]);\n", index(o + N1 * k2, q, 1, N1 - 1), -b[q].imag(), -b[q].real(), index(o + N1 * k2, q, 1, N1 - 1));
		}
	}

	for (int k2 = 0; k2 < N2; k2++) {
		fft_bitreverse(N1 - 1, fft_bitreverse_indices(N1 - 1), o + k2 * N1);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		print("{\n");
		indent();
		fft(N1 - 1, o + k2 * N1, true);
		deindent();
		print("}\n");
	}
	for (int k2 = 0; k2 < N2; k2++) {
		const auto Nm1inv = 1.0 / (N1 - 1.0);
		for (int q = 0; q < N1 - 1; q++) {
			print("x[%i] *= (%24.17e);\n", index(o + N1 * k2, q, 0, N1 - 1), Nm1inv);
			print("x[%i] = -x[%i] * (%24.17e);\n", index(o + N1 * k2, q, 1, N1 - 1), index(o + N1 * k2, q, 1, N1 - 1), Nm1inv);
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int p = N1 - 2; p >= 0; p--) {
			print("x[%i] = xro[%i] + x[%i];\n", index(o + N1 * k2, p + 1, 0, N1 - 1), k2, index(o + N1 * k2, p, 0, N1 - 1));
			print("x[%i] = xio[%i] + x[%i];\n", index(o + N1 * k2, p + 1, 1, N1 - 1), k2, index(o + N1 * k2, p, 1, N1 - 1));
		}
		print("x[%i] = xrk0[%i];\n", index(o + k2 * N1, 0, 0, N1 - 1), k2);
		print("x[%i] = xik0[%i];\n", index(o + k2 * N1, 0, 1, N1 - 1), k2);
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



int raders_fft_opcnt(int r, int N) {
	int cnt = 0;
	const int N1 = r;
	const int N2 = N / r;
	if (N2 > 1) {
		cnt += MWEIGHT * N;
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
		cnt += MWEIGHT * (N1 - 1);
		cnt += fft_opcnt(N1 - 1);
	}
	for (int k2 = 0; k2 < N2; k2++) {
		cnt++;
		for (int q = 1; q < N1 - 1; q++) {
			cnt += 2;
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		cnt += MWEIGHT * (N1 - 1);
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
	cnt += MWEIGHT * N;
	return cnt;
}

