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
	for (int n1 = 0; n1 < N1; n1++) {
		//	fft_bitreverse(N2, fft_bitreverse_indices(N2), o + n1 * N2);
		fft(N2, n1 * N2 + o);
	}
	print("std::array<std::array<double, %i>, %i> z;\n", 2 * N1, N2);
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			print("z[%i][%i] = x[%i];\n", n2, 2 * n1 + 0, index(N2 * n1 + o, n2, 0, N1));
			print("z[%i][%i] = x[%i];\n", n2, 2 * n1 + 1, index(N2 * n1 + o, n2, 1, N1));
		}
		print("{\n");
		indent();
		print("auto* x = z[%i].data();\n", n2);
		fft_bitreverse(N1, fft_bitreverse_indices(N1), 0);
		fft(N1, 0);
		deindent();
		print("}\n");
	}
	for (int n = 0; n < N; n++) {
		print("x[%i] = z[%i][%i];\n", index(o, n, 0, N), (n % N2), 2 * (n % N1) + 0);
		print("x[%i] = z[%i][%i];\n", index(o, n, 1, N), (n % N2), 2 * (n % N1) + 1);
	}
	deindent();
	print("}\n");
}

std::vector<int> gt2_fft_bitr(int N1, int N2, int o, std::vector<int> I) {
	int N = N1 * N2;
	std::vector<int> L;
	for (int n1 = 0; n1 < N1; n1++) {
		std::vector<int> J;
		for (int n2 = 0; n2 < N2; n2++) {
			J.push_back(I[(n1 * N2 + N1 * n2) % N]);
		}
		auto K = fft_bitr(N2, o + n1 * N2, J);
		L.insert(L.end(), K.begin(), K.end());
	}
	return L;
}

