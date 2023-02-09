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

const std::vector<std::complex<double>> raders_four_twiddle(int N, int M) {
	std::vector<std::complex<double>> b(M, 0);
	const auto tws = twiddles(N);
	const auto ginvq = raders_ginvq(N);
	b[0] = tws[ginvq[0]];
	for (int q = 1; q < N - 1; q++) {
		b[M + q - (N - 1)] = b[q] = tws[ginvq[q]];
	}
	fftw(b);
	return b;
}
/*
void bitreverse2(int N, int o) {
	for (int n = 0; n < N; n++) {
		int j = n;
		int i = 0;
		for (int k = 0; k < ilogb(N); k++) {
			i <<= 1;
			if (j & 1) {
				i |= 1;
			}
			j >>= 1;
		}
		if (n > i) {
			print("std::swap(x[%i], x[%i]);\n", 2 * o + 2 * n, 2 * o + 2 * i);
			print("std::swap(x[%i], x[%i]);\n", 2 * o + 2 * n + 1, 2 * o + 2 * i + 1);
		}
	}
}

void fft2(int N, int o, std::vector<bool>& zero) {
	if (N == 1) {
		return;
	}
	int r = 2;
	print("{\n");
	print("// radix - %i\n", r);
	indent();
	if (N > 2) {
		std::vector<bool> zerol, zeror;
		for (int n = 0; n < N / 2; n++) {
			zerol.push_back(zero[2 * n]);
			zeror.push_back(zero[2 * n + 1]);
		}
		fft2(N / 2, o, zerol);
		fft2(N / 2, o + N / 2, zeror);
		for (int n = 0; n < N / 2; n++) {
			zero[2 * n] = zero[2 * n] && zerol[n];
			zero[2 * n + 1] = zero[2 * n + 1] && zeror[n];
		}
	}
	for (int k = 0; k < N / r; k++) {
		print("{\n");
		indent();
		if (zero[k] && zero[k + 1]) {
		} else if (!zero[k] && zero[k + 1]) {
			print_z(0, 0, k, 2, N, o);
			print("x[%i] = zr0;\n", index(o, k + N / 2, 0, N));
			print("x[%i] = zi0;\n", index(o, k + N / 2, 1, N));
			zero[k] = false;
			zero[k + 1] = false;
		} else if (zero[k] && !zero[k + 1]) {
			print_z(1, 1, k, 2, N, o);
			print("x[%i] = zr1;\n", index(o, k, 0, N));
			print("x[%i] = zi1;\n", index(o, k, 1, N));
			print("x[%i] = -zr1;\n", index(o, k + N / 2, 0, N));
			print("x[%i] = -zi1;\n", index(o, k + N / 2, 1, N));
			zero[k] = false;
			zero[k + 1] = false;
		} else {
			print_z(0, 0, k, 2, N, o);
			print_z(1, 1, k, 2, N, o);
			print("x[%i] = zr0 + zr1;\n", index(o, k, 0, N));
			print("x[%i] = zi0 + zi1;\n", index(o, k, 1, N));
			print("x[%i] = zr0 - zr1;\n", index(o, k + N / 2, 0, N));
			print("x[%i] = zi0 - zi1;\n", index(o, k + N / 2, 1, N));
			zero[k] = false;
			zero[k + 1] = false;
		}
		deindent();
		print("}\n");
	}
	deindent();
	print("}\n");
}

void raders_fft_padded(int N, int o) {
	print("{\n");
	indent();
	int m = 2 * N - 1;
	int M = m;
	while (greatest_prime_factor(M) > 2) {
		M++;
	}
	const auto b = raders_four_twiddle(N, M);
	const auto gq = raders_gq(N);
	const auto ginvq = raders_ginvq(N);
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	print("std::array<double, %i> z;\n", 2 * M);
	std::vector<bool> zeros;
	for (int n = 0; n < N - 1; n++) {
		print("z[%i] = x[%i];\n", 2 * n, index(o, gq[n], 0, N - 1));
		print("z[%i] = x[%i];\n", 2 * n + 1, index(o, gq[n], 1, N - 1));
		zeros.push_back(false);
	}
	for (int n = N - 1; n < M; n++) {
		print("z[%i] = 0;\n", 2 * n);
		print("z[%i] = 0;\n", 2 * n + 1);
		zeros.push_back(true);
	}
	print("{\n");
	indent();
	print("auto* x = z.data();\n");
	bitreverse2(M, 0);
	fft2(M, 0, zeros);
	deindent();
	print("}\n");
	for (int m = 0; m < M; m++) {
		print("tmp0 = z[%i];\n", 2 * m);
		print("z[%i] = std::fma((%.17e), z[%i], (%.17e) * z[%i]);\n", 2 * m, b[m].real(), 2 * m, -b[m].imag(), 2 * m + 1);
		print("z[%i] = std::fma((%.17e), tmp0, (%.17e) * z[%i]);\n", 2 * m + 1, -b[m].imag(), -b[m].real(), 2 * m + 1);
	}
	print("{\n");
	indent();
	print("auto* x = z.data();\n");
	bitreverse2(M, 0);
	auto z = std::vector<bool>(M, false);
	fft2(M, 0, z);
	deindent();
	print("}\n");
	for (int m = 0; m < (N - 1); m++) {
		print("z[%i] *= %.17e;\n", 2 * m, 1.0 / M);
		print("z[%i] = z[%i] * (%.17e);\n", 2 * m + 1, 2 * m + 1, -1.0 / M);
	}
	print("double xkor = x[%i];\n", index(o, 0, 0, N));
	print("double xkoi = x[%i];\n", index(o, 0, 1, N));
	for (int m = 1; m < N; m++) {
		print("x[%i] += x[%i];\n", index(o, 0, 0, N), index(o, m, 0, N));
		print("x[%i] += x[%i];\n", index(o, 0, 1, N), index(o, m, 1, N));
	}
	for (int m = 0; m < N - 1; m++) {
		print("x[%i] = xkor + z[%i];\n", index(o, ginvq[m], 0, N), 2 * m);
		print("x[%i] = xkoi + z[%i];\n", index(o, ginvq[m], 1, N), 2 * m + 1);
	}
	deindent();
	print("}\n");
}
*/
void raders_fft(int N, int o, bool padded) {
//	raders_fft_padded(N, o);
//	return;
	print("{\n");
	indent();
	const auto b = raders_four_twiddle(N);
	const auto gq = raders_gq(N);
	const auto ginvq = raders_ginvq(N);
	print("// Raders radix - %i \n", N);
	printf("// Raders radix - %i \n", N);
	print("double xro, xrk0;\n");
	print("double xio, xik0;\n");
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	print("xro = x[%i];\n", index(o, 0, 0, N));
	print("xrk0 = x[%i] ", index(o, 0, 0, N));
	for (int k1 = 1; k1 < N; k1++) {
		print_notab(" + x[%i]", index(o, k1, 0, N));
	}
	print_notab(";\n");
	print("xio = x[%i];\n", index(o, 0, 1, N));
	print("xik0 = x[%i] ", index(o, 0, 1, N));
	for (int k1 = 1; k1 < N; k1++) {
		print_notab(" + x[%i]", index(o, k1, 1, N));
	}
	print_notab(";\n");
	std::vector<int> I1(N), I3(N);
	auto I2 = fft_bitreverse_indices(N - 1);
	for (int i = 0; i < I2.size(); i++) {
		I3[I2[i]] = i;
	}
	for (int q = 0; q < N - 1; q++) {
		I1[I3[q]] = gq[q];
	}
	I1[N - 1] = 0;
	fft_bitreverse(N, I1, o);

	print("{\n");
	indent();
	fft(N - 1, o, true);
	deindent();
	print("}\n");
	print("x[%i] = -x[%i];\n", index(o, 0, 0, N - 1), index(o, 0, 0, N - 1));
	for (int q = 1; q < N - 1; q++) {
		print("tmp0 = x[%i];\n", index(o, q, 0, N - 1));
		print("x[%i] = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", index(o, q, 0, N - 1), b[q].real(), index(o, q, 0, N - 1), -b[q].imag(), index(o, q, 1, N - 1));
		print("x[%i] = std::fma((%.17e), tmp0, (%.17e) * x[%i]);\n", index(o, q, 1, N - 1), -b[q].imag(), -b[q].real(), index(o, q, 1, N - 1));
	}
	fft_bitreverse(N - 1, fft_bitreverse_indices(N - 1), o);
	print("{\n");
	indent();
	fft(N - 1, o, true);
	deindent();
	print("}\n");
	const auto Nm1inv = 1.0 / (N - 1.0);
	for (int q = 0; q < N - 1; q++) {
		print("x[%i] *= (%24.17e);\n", index(o, q, 0, N - 1), Nm1inv);
		print("x[%i] = x[%i] * (%24.17e);\n", index(o, q, 1, N - 1), index(o, q, 1, N - 1), -Nm1inv);
	}
	for (int p = N - 2; p >= 0; p--) {
		print("x[%i] = xro + x[%i];\n", index(o, p + 1, 0, N - 1), index(o, p, 0, N - 1));
		print("x[%i] = xio + x[%i];\n", index(o, p + 1, 1, N - 1), index(o, p, 1, N - 1));
	}
	print("x[%i] = xrk0;\n", index(o, 0, 0, N - 1));
	print("x[%i] = xik0;\n", index(o, 0, 1, N - 1));
	std::vector<int> indices;
	indices.resize(N);
	indices[0] = 0;
	for (int k1 = 0; k1 < N - 1; k1++) {
		indices[ginvq[k1]] = k1 + 1;
	}
	fft_bitreverse(N, indices, o);
	deindent();
	print("}\n");
}

int raders_fft_opcnt(int N) {
	int cnt = 0;
	for (int n1 = 0; n1 < N; n1++) {
		cnt += 2;
	}
	for (int k1 = 1; k1 < N; k1++) {
		cnt += 2;
	}
	cnt += MWEIGHT * (N - 1);
	cnt += fft_opcnt(N - 1);
	cnt++;
	for (int q = 1; q < N - 1; q++) {
		cnt += 2;
	}
	cnt += MWEIGHT * (N - 1);
	cnt += fft_opcnt(N - 1);
	const auto Nm1inv = 1.0 / (N - 1.0);
	for (int q = 0; q < N - 1; q++) {
		cnt += 3;
	}
	for (int p = N - 2; p >= 0; p--) {
		cnt += 2;
	}
	cnt += MWEIGHT * N;
	return cnt;
}

