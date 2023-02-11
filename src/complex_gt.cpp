#include "fftgen.hpp"

void gt_fft(int N, int o);

int gt2_fft_opcnt(int N1, int N2) {
	int N = N1 * N2;
	int cnt = 0;
	cnt += N2 * fft_opcnt(N1);
	cnt += N1 * fft_opcnt(N2);
	cnt += MWEIGHT * N;
	return cnt;
}

void gt2_fft(int N1, int N2, int o) {
	printf("// good-thomas - %i = %i x %i\n", N1 * N2, N1, N2);
	print("{\n");
	indent();
	int N = N1 * N2;
	for (int n2 = 0; n2 < N2; n2++) {
		print("fft_base_%i(x + %i);\n", N1, 2 * (o + n2 * N1));
	}
	auto J = fft_bitreverse_indices(N2);
	std::vector<int> I(N);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			I[n1 * N2 + n2] = N1 * J[n2] + n1;
		}
	}
	fft_bitreverse(N, I, o);
	std::vector<int> I(N);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
			I[n1 * N2 + n2] = (n % N2) + (n % N1) * N2;
		}
	}
	fft_bitreverse(N, I, o);
	for (int n1 = 0; n1 < N1; n1++) {
		print("fft_base_%i(x + %i);\n", N2, 2 * (o + n1 * N2));
	}
	/*for (int k1 = 0; k1 < N1; k1++) {
		for (int k2 = 0; k2 < N2; k2++) {
			int k = k1 * N2 + k2;
			I[k1 * N2 + k2] = (k % N2) + (k % N1) * N2;
		}
	}*/
	fft_bitreverse(N, I, o);
	deindent();
	print("}\n");
}

std::vector<int> gt2_fft_bitr(int N1, int N2, int o, std::vector<int> I) {
	int N = N1 * N2;
	std::vector<int> L;
	for (int n2 = 0; n2 < N2; n2++) {
		std::vector<int> J;
		for (int n1 = 0; n1 < N1; n1++) {
			J.push_back(I[(n1 * N2 + N1 * n2) % N]);
		}
		auto K = fft_bitr(N1, o + n2 * N1, J);
		L.insert(L.end(), K.begin(), K.end());
	}
	return L;
}

