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

void raders_fft(int N, int o, bool padded) {
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
	print("xro = x[%i];\n", index(o, N - 1, 0, N));
	print("xio = x[%i];\n", index(o, N - 1, 1, N));
	print("xrk0 = x[%i];\n ", index(o, N - 1, 0, N));
	print("xik0 = x[%i];\n", index(o, N - 1, 1, N));
	print("for (int k1 = 0; k1 < %i; k1++) {\n", N - 1);
	indent();
	print("const int ir = 2 * k1;\n");
	print("xrk0 += x[ir];\n");
	print("xik0 += x[ir + 1];\n");
	deindent();
	print("}\n");
	std::vector<int> I1(N), I3(N);
	print("fft_base_%i(x + %i);\n", N - 1, 2 * o);
	print("x[%i] = -x[%i];\n", index(o, 0, 0, N - 1), index(o, 0, 0, N - 1));
	if (N >= LOOP_N) {
		print("const static raders_twiddles twiddles(%i);\n", N);
		print("for( int q = 1; q < %i; q++ ) {\n", N - 1);
		indent();
		print("const auto tw = twiddles[q];\n");
		print("const auto qr = 2 * q;\n");
		print("const auto qi = 2 * q + 1;\n");
		print("tmp0 = x[qr];\n");
		print("x[qr] = tw.real() * x[qr] - tw.imag() * x[qi];\n");
		print("x[qi] = -std::fma(tw.imag(), tmp0, tw.real() * x[qi]);\n");
		deindent();
		print("}\n");
	} else {
		for (int q = 1; q < N - 1; q++) {
			print("tmp0 = x[%i];\n", index(o, q, 0, N - 1));
			print("x[%i] = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", index(o, q, 0, N - 1), b[q].real(), index(o, q, 0, N - 1), -b[q].imag(), index(o, q, 1, N - 1));
			print("x[%i] = std::fma((%.17e), tmp0, (%.17e) * x[%i]);\n", index(o, q, 1, N - 1), -b[q].imag(), -b[q].real(), index(o, q, 1, N - 1));
		}
	}
	print("fft_%i(x + %i);\n", N - 1, 2 * o);
	const auto Nm1inv = 1.0 / (N - 1.0);
	print("for (int q = 0; q < %i; q++) {\n", N - 1);
	indent();
	print("const int ir = 2 * q;\n");
	print("const int ii = ir + 1;\n");
	print("x[ir] *= (%24.17e);\n", Nm1inv);
	print("x[ii] *= (%24.17e);\n", -Nm1inv);
	deindent();
	print("}\n");
	if (N >= LOOP_N) {
		print("for (int p = %i; p >= 0; p -= 2) {\n", 2 * (N - 2));
		indent();
		print("x[p + 3] = xio + x[p + 1];\n");
		print("x[p + 2] = xro + x[p];\n");
		deindent();
		print("}\n");
	} else {
		for (int p = N - 2; p >= 0; p--) {
			print("x[%i] = xro + x[%i];\n", index(o, p + 1, 0, N - 1), index(o, p, 0, N - 1));
			print("x[%i] = xio + x[%i];\n", index(o, p + 1, 1, N - 1), index(o, p, 1, N - 1));
		}
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

