#include "fftgen.hpp"
#include <list>



std::vector<int> fft_bitreverse_indices_real(int N) {
	std::vector<int> I;
	std::vector<int> indices;
	I.resize(0);
	for (int n = 0; n < N; n++) {
		I.push_back(n);
	}
	indices = fft_bitr_real(N, 0, I, true);
	return indices;
}

std::vector<int> fft_bitr_real(int N, int o, std::vector<int> indices, bool first) {
	if (N == 1) {
		return indices;
	}
	auto fftt = best_radix_real(N, o, first);
//	if (fftt.type == RADIX) {
	return fft_radix_bitr_real(fftt.N1, N, o, indices);
//	} else if (fftt.type == GOOD) {
//		return gt2_fft_bitr(fftt.N1, fftt.N2, o, indices);
//	} else {
//		return indices;
	//}
}


fft_type best_radix_real(int N, int o, bool first) {
	fft_type fftt;
	int best_cnt = 999999999;
	int best_radix = -1;
	for (int r = 2; r <= N; r++) {
		if (N % r == 0) {
			int this_cnt;
			if (r <= 3 || is_prime(r)) {
				this_cnt = fft_radix_real_opcnt(r, N, o);
				if (first) {
					this_cnt += N * MWEIGHT;
				}
				if (this_cnt < best_cnt) {
					best_cnt = this_cnt;
					fftt.N1 = r;
					fftt.type = RADIX;
					fftt.nops = best_cnt;
				}
			}
		}
	}
	return fftt;
}

int fft_radix_real_opcnt(int r, int N, int o) {
	int cnt = 0;
	if (N <= 1) {
		return cnt;
	}
	for (int n = 0; n < r; n++) {
		cnt += fft_real_opcnt(N / r, o + n * N / r);
	}
	for (int k = 0; k < N / r; k++) {
		if (k > N / r / 2) {
			continue;
		}
		std::vector<bool> zr_zero(r, false);
		std::vector<bool> zi_zero(r, false);
		if (k == 0 || (k == N / r / 2 && N / r % 2 == 0)) {
			for (int i = 0; i < r; i++) {
				const auto W = twiddle(i * k, N);
				const int L = N / r;
				const int M = L / 2;
				const int oo = o + i * N / r;
				const int iir = index_real(oo, k, 0, L);
				const int iii = index_real(oo, k, 1, L);
				if (i * k == 0) {
					zi_zero[i] = true;
				} else if (i * k == N / 4 && N % 4 == 0) {
					cnt++;
					zr_zero[i] = true;
				} else if (i * k == N / 2 && N % 2 == 0) {
					zi_zero[i] = true;
					cnt++;
				} else if (i * k == 3 * N / 4 && N % 4 == 0) {
					zr_zero[i] = true;
				} else {
					cnt += 2;
				}
			}
		} else {
			const int M = N / 2;
			for (int i = 0; i < r; i++) {
				const auto W = twiddle(i * k, N);
				const int L = N / r;
				const int M = L / 2;
				const int oo = o + i * N / r;
				const int iir = index_real(oo, k, 0, L);
				const int iii = index_real(oo, k, 1, L);
				if (i * k == 0) {
				} else if (i * k == N / 4 && N % 4 == 0) {
					cnt++;
				} else if (i * k == N / 2 && N % 2 == 0) {
					cnt += 2;
				} else if (i * k == 3 * N / 4 && N % 4 == 0) {
					cnt += 1;
				} else {
					cnt += 4;
				}
			}
		}

		switch (r) {
		case 2:
			cnt += int(!(zr_zero[0] || zr_zero[1]));
			if (N / 2 - k != k) {
				cnt += int(!(zr_zero[0] || zr_zero[1]));
			}
			if (k != 0) {
				cnt += int(!(zi_zero[0] || zi_zero[1]));
				if (N / 2 - k != k) {
					cnt += int(!zi_zero[0]);
				}
			}
			break;
		case 3:
			cnt += int(!(zr_zero[1] || zr_zero[2]));
			cnt += int(!(zi_zero[1] || zi_zero[2]));
			cnt += 2;
			cnt += 1 + int(!zr_zero[2]);
			cnt += 1 + int(!zi_zero[2]);
			const int k0 = k;
			const int k1 = k <= N / 2 ? k + N / 3 : N - (k + N / 3);
			const int k2 = N - (k + 2 * N / 3);
			if (k1 != k2) {
				cnt++;
			}
			if (k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
				cnt++;
			}
			if (k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
				cnt++;
			}
			if (k1 != k2) {
				if (k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					cnt++;
				}
			}
			break;
		};
	}
	return cnt;
}



std::vector<int> fft_radix_bitr_real(int r, int N, int o, std::vector<int> I) {
	std::vector<int> L;
	std::vector<int> J;
	std::vector<int> K;
	switch (r) {
	case 1:
		return I;
	case 4: {
		std::vector<int> J;
		for (int k = 0; k < N / 2; k++) {
			J.push_back(I[2 * k]);
		}
		auto K = fft_bitr_real(N / 2, o, J);
		L.insert(L.end(), K.begin(), K.end());
	}
		for (int n = 1; n < 4; n += 2) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + 4 * k]);
			}
			auto K = fft_bitr_real(N / 4, o + (2 + (n / 2)) * N / 4, J);
			L.insert(L.end(), K.begin(), K.end());
		}
		break;
	default:
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr_real(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
	};
	return L;
}

void fft_radix_real(int r, int N, int o) {
	if (N < r) {
		return;
	}
	print("{\n");
	indent();
	for (int n = 0; n < r; n++) {
		fft_real(N / r, o + n * N / r);
	}
	for (int k = 0; k < N / r; k++) {
		if (k > N / r / 2) {
			continue;
		}
		print("{\n");
		indent();
		std::vector<bool> zr_zero(r, false);
		std::vector<bool> zi_zero(r, false);
		if (k == 0 || (k == N / r / 2 && N / r % 2 == 0)) {
			for (int i = 0; i < r; i++) {
				const auto W = twiddle(i * k, N);
				const int L = N / r;
				const int M = L / 2;
				const int oo = o + i * N / r;
				const int iir = index_real(oo, k, 0, L);
				const int iii = index_real(oo, k, 1, L);
				if (i * k == 0) {
					print("const auto zr%i = x[%i];\n", i, iir);
					print("constexpr double zi%i = 0;\n", i);
					zi_zero[i] = true;
				} else if (i * k == N / 8 && N % 8 == 0) {
					print("const auto zr%i = M_SQRT1_2 * (x[%i]);\n", i, iir);
					print("const auto zi%i = -M_SQRT1_2 * (x[%i]);\n", i, iir);
				} else if (i * k == N / 4 && N % 4 == 0) {
					print("constexpr double zr%i = 0;\n", i);
					print("const auto zi%i = -x[%i];\n", i, iir);
					zr_zero[i] = true;
				} else if (i * k == 3 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = -M_SQRT1_2 * x[%i];\n", i, iir);
					print("const auto zi%i = -M_SQRT1_2 * x[%i];\n", i, iir);
				} else if (i * k == N / 2 && N % 2 == 0) {
					print("const auto zr%i = -x[%i];\n", i, iir);
					print("constexpr double zi%i = 0;\n", i);
					zi_zero[i] = true;
				} else if (i * k == 5 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = -M_SQRT1_2 * x[%i];\n", i, iir);
					print("const auto zi%i = M_SQRT1_2 * x[%i];\n", i, iir);
				} else if (i * k == 3 * N / 4 && N % 4 == 0) {
					print("constexpr double zr%i = 0;\n", i);
					print("const auto zi%i = x[%i];\n", i, iir);
					zr_zero[i] = true;
				} else if (i * k == 7 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = M_SQRT1_2 * (x[%i]);\n", i, iir);
					print("const auto zi%i = M_SQRT1_2 * (x[%i]);\n", i, iir);
				} else {
					print("const auto zr%i = (%.17e) * x[%i];\n", i, W.real(), iir);
					print("const auto zi%i = (%.17e) * x[%i];\n", i, W.imag(), iir);
				}
			}
		} else {
			const int M = N / 2;
			for (int i = 0; i < r; i++) {
				const auto W = twiddle(i * k, N);
				const int L = N / r;
				const int M = L / 2;
				const int oo = o + i * N / r;
				const int iir = index_real(oo, k, 0, L);
				const int iii = index_real(oo, k, 1, L);
				if (i * k == 0) {
					print("const auto zr%i = x[%i];\n", i, iir);
					print("const auto zi%i = x[%i];\n", i, iii);
				} else if (i * k == N / 8 && N % 8 == 0) {
					print("const auto zr%i = M_SQRT1_2 * (x[%i] + x[%i]);\n", i, iir, iii);
					print("const auto zi%i = -M_SQRT1_2 * (x[%i] - x[%i]);\n", i, iir, iii);
				} else if (i * k == N / 4 && N % 4 == 0) {
					print("const auto zr%i = x[%i];\n", i, iii);
					print("const auto zi%i = -x[%i];\n", i, iir);
				} else if (i * k == 3 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = -M_SQRT1_2 * (x[%i] - x[%i]);\n", i, iir, iii);
					print("const auto zi%i = -M_SQRT1_2 * (x[%i] + x[%i]);\n", i, iir, iii);
				} else if (i * k == N / 2 && N % 2 == 0) {
					print("const auto zr%i = -x[%i];\n", i, iir);
					print("const auto zi%i = -x[%i];\n", i, iii);
				} else if (i * k == 5 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = -M_SQRT1_2 * (x[%i] + x[%i]);\n", i, iir, iii);
					print("const auto zi%i = M_SQRT1_2 * (x[%i] - x[%i]);\n", i, iir, iii);
				} else if (i * k == 3 * N / 4 && N % 4 == 0) {
					print("const auto zr%i = -x[%i];\n", i, iii);
					print("const auto zi%i = x[%i];\n", i, iir);
				} else if (i * k == 7 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = M_SQRT1_2 * (x[%i] - x[%i]);\n", i, iir, iii);
					print("const auto zi%i = M_SQRT1_2 * (x[%i] + x[%i]);\n", i, iir, iii);
				} else {
					print("const auto zr%i = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", i, W.real(), iir, -W.imag(), iii);
					print("const auto zi%i = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", i, W.real(), iii, W.imag(), iir);
				}
			}
		}

		switch (r) {
		case 2:
			print("// Radix - 2\n");
			print("x[%i] = zr0 + zr1;\n", index_real(o, k, 0, N));
			if (N / 2 - k != k) {
				print("x[%i] = zr0 - zr1;\n", index_real(o, N / 2 - k, 0, N));
			}
			if (k != 0) {
				print("x[%i] = zi0 + zi1;\n", index_real(o, k, 1, N));
				if (N / 2 - k != k) {
					print("x[%i] = zi1 - zi0;\n", index_real(o, N / 2 - k, 1, N));
				}
			}
			break;
		case 3:
			print("// Radix - 3\n");
			print("const auto tr1 = zr1 + zr2;\n");
			print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
			print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
			print("const auto ti1 = zi1 + zi2;\n");
			print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
			print("const auto ti3 = (%24.17e) * (zi1 - zi2);\n", sqrt(3) * 0.5);
			{
				const int k0 = k;
				const int k1 = k <= N / 2 ? k + N / 3 : N - (k + N / 3);
				const int k2 = N - (k + 2 * N / 3);
				print("x[%i] = zr0 + tr1;\n", index_real(o, k0, 0, N));
				print("x[%i] = tr2 + ti3;\n", index_real(o, k1, 0, N));
				if (k1 != k2) {
					print("x[%i] = tr2 - ti3;\n", index_real(o, k2, 0, N));
				}
				if (k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					print("x[%i] = zi0 + ti1;\n", index_real(o, k0, 1, N));
				}
				if (k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					print("x[%i] = ti2 - tr3;\n", index_real(o, k1, 1, N));
				}
				if (k1 != k2) {
					if (k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
						print("x[%i] = -(ti2 + tr3);\n", index_real(o, k2, 1, N));
					}
				}
			}
			break;
		case 5:
			// Radix -- 5
			print("// Radix - 5\n");
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
			{
				int k0 = k;
				int k1 = k + N / 5;
				int k2 = k + 2 * N / 5;
				int k3 = N - (k + 3 * N / 5);
				int k4 = N - (k + 4 * N / 5);
				if (k > N / 2) {
					k2 = N - k2;
				}
				print("x[%i] = zr0 + tr5;\n", index_real(o, k0, 0, N));
				print("x[%i] = tr8 + ti10;\n", index_real(o, k1, 0, N));
				print("x[%i] = tr9 + ti11;\n", index_real(o, k2, 0, N));
				print("x[%i] = tr9 - ti11;\n", index_real(o, k3, 0, N));
				print("x[%i] = tr8 - ti10;\n", index_real(o, k4, 0, N));
				if (k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					print("x[%i] = zi0 + ti5;\n", index_real(o, k0, 1, N));
				}
				if (k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					print("x[%i] = ti8 - tr10;\n", index_real(o, k1, 1, N));
				}
				if (k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					print("x[%i] = ti9 - tr11;\n", index_real(o, k2, 1, N));
				}
				if (k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					print("x[%i] = -(ti9 + tr11);\n", index_real(o, k3, 1, N));
				}
				if (k4 != 0 && !(k4 == N / 2 && N % 2 == 0)) {
					print("x[%i] = -(ti8 + tr10);\n", index_real(o, k4, 1, N));
				}
			}
			break;
		default:
			print("// radix - %i\n", r);
			for (int j = 1; j <= (r - 1) / 2; j++) {
				print("const auto txp%i = zr%i + zr%i;\n", j, j, r - j);
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				print("const auto txm%i = zr%i - zr%i;\n", j, j, r - j);
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				print("const auto typ%i = zi%i + zi%i;\n", j, j, r - j);
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				print("const auto tym%i = zi%i - zi%i;\n", j, j, r - j);
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				print("auto ap%i = zr0;\n", i);
				print("auto bp%i = zi0;\n", i);
				for (int j = 1; j <= (r - 1) / 2; j++) {
					print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
					print("bp%i = std::fma(typ%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
					if (j == 1) {
						print("auto am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
						print("auto bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
					} else {
						print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
						print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
					}
				}
			}
			print("x[%i] = ", index_real(o, k, 0, N));
			for (int i = 0; i < r; i++) {
				print_notab(" + zr%i", i);
			}
			print_notab(";\n");
			for (int i = 1; i <= (r - 1) / 2; i++) {
				const int ka = k + i * N / r;
				const int rk0 = index_real(o, ka, 0, N);
				print("x[%i] = ap%i + am%i; // 1 %i\n", rk0, i, i, i);
			}
			if (!(k == N / r / 2 && N % 2 == 0)) {
				for (int i = 1; i <= (r - 1) / 2; i++) {
					const int ka = k + i * N / r;
					const int kb = N - (k + (r - i) * N / r);
					if (ka != kb) {
						const int rk1 = index_real(o, kb, 0, N);
						print("x[%i] = ap%i - am%i; // 2 %i\n", rk1, i, i, i);
					}
				}
			}
			if (k != 0 && !(k == N / 2 && N % 2 == 0)) {
				print("x[%i] = ", index_real(o, k, 1, N));
				for (int i = 0; i < r; i++) {
					print_notab(" + zi%i", i);
				}
				print_notab(";\n");
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				const int ka = k + i * N / r;
				const int ik0 = index_real(o, ka, 1, N);
				if (ka != 0 && !(ka == N / 2 && N % 2 == 0)) {
					print("x[%i] = bp%i - bm%i; // 3 %i\n", ik0, i, i, i);
				}
			}
			if (!(k == N / r / 2 && N % 2 == 0)) {
				for (int i = 1; i <= (r - 1) / 2; i++) {
					const int ka = k + i * N / r;
					const int kb = N - (k + (r - i) * N / r);
					const int ik1 = index_real(o, kb, 1, N);
					if (ka != kb) {
						if (kb != 0 && !(kb == N / 2 && N % 2 == 0)) {
							print("x[%i] = -(bp%i + bm%i); // 6 %i\n", ik1, i, i, i);
						}
					}
				}
			}
		}
		deindent();
		print("}\n");
	}

	deindent();
	print("}\n");

}

void fft_bitreverse_real(int N, std::vector<int> indices, int o) {
	if (indices.size() == 0) {
		indices = fft_bitreverse_indices_real(N);
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


