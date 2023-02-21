#include "fftgen.hpp"

void gt_fft(int N, int o);

int gt2_fft_opcnt(int N1, int N2) {
	int N = N1 * N2;
	int cnt = 0;
	cnt += N2 * fft_opcnt(N1);
	cnt += N1 * fft_opcnt(N2);
	cnt += MWEIGHT * N;
	return 2 * cnt;
}

void gt2_fft(int N1, int N2, int o) {
	printf("// good-thomas - %i = %i x %i\n", N1 * N2, N1, N2);
	print("{\n");
	indent();
	int N = N1 * N2;
	print("auto* y = x;");
	for (int n2 = 0; n2 < N2; n2++) {
		print("{\n");
		indent();
		print("auto* x = y + %i;\n", 2 * n2);
		fft_radix_dit_strided(best_radix(N1, 0).N1, N1, o, N2);
//		print("fft_base_%i(x + %i);\n", N1, 2 * (o + n2 * N1));
		deindent();
		print("}\n");
	}
	for (int n1 = 0; n1 < N1; n1++) {
		print("fft_base_%i(x + %i);\n", N2, 2 * (o + n1 * N2));
	}
	std::vector<int> I(N);
	for (int k1 = 0; k1 < N1; k1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			int k = k1 * N2 + k2;
			I[k1 * N2 + k2] = (k % N2) + (k % N1) * N2;
		}
	}
	fft_bitreverse(N, I, o);
	deindent();
	print("}\n");
}

std::vector<int> gt2_fft_bitr(int N1, int N2, int o, std::vector<int> I) {
	int N = N1 * N2;
	std::vector<int> J(N);
	std::vector<int> K(N);
	std::vector<int> L(N);
	auto P = fft_bitreverse_indices(N2);
	auto Q = fft_bitreverse_indices(N1);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			J[Q[n1] * N2 + n2] = I[(n1 * N2 + P[n2] * N1) % N];
		}
	}
	return J;
}

