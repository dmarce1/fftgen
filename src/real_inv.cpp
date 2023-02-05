#include "fftgen.hpp"

int fft_radix_real_inv_opcnt(int r, int N, int o);
void fft_radix_real_inv(int r, int N, int o);
std::vector<int> fft_radix_bitr_real_inv(int r, int N, int o, std::vector<int> I);

fft_type best_radix_real_inv(int N, int o, bool first = false) {
	fft_type fftt;
	int best_cnt = 999999999;
	int best_radix = -1;
	for (int r = 2; r <= N; r++) {
		if (N % r == 0) {
			int this_cnt;
			if (r <= 6 || is_prime(r)) {
				this_cnt = fft_radix_real_inv_opcnt(r, N, o);
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

void fft_real_inv(int N, int o) {
	if (N > 1) {
		fft_type fftt = best_radix_real_inv(N, o);
		fft_radix_real_inv(fftt.N1, N, o);
	}
}

int fft_real_inv_opcnt(int N, int o) {
	fft_type fftt = best_radix_real_inv(N, o);
	return fft_radix_real_inv_opcnt(fftt.N1, N, o);
}

int fft_radix_real_inv_opcnt(int r, int N, int o) {
	int cnt = 0;
	if (N <= 1) {
		return cnt;
	}
	for (int n = 0; n < r; n++) {
		cnt += fft_real_inv_opcnt(N / r, o + n * N / r);
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

std::vector<int> fft_bitreverse_indices_real_inv(int N) {
	std::vector<int> I;
	std::vector<int> indices;
	I.resize(0);
	for (int n = 0; n < N; n++) {
		I.push_back(n);
	}
	indices = fft_bitr_real_inv(N, 0, I, true);
	return indices;
}

std::vector<int> fft_bitr_real_inv(int N, int o, std::vector<int> indices, bool first) {
	if (N == 1) {
		return indices;
	}
	auto fftt = best_radix_real_inv(N, o, first);
//	if (fftt.type == RADIX) {
	return fft_radix_bitr_real_inv(fftt.N1, N, o, indices);
//	} else if (fftt.type == GOOD) {
//		return gt2_fft_bitr(fftt.N1, fftt.N2, o, indices);
//	} else {
//		return indices;
	//}
}

std::vector<int> fft_radix_bitr_real_inv(int r, int N, int o, std::vector<int> I) {
	std::vector<int> L;
	std::vector<int> J;
	std::vector<int> K;
	switch (r) {
	case 1:
		return I;
		/*case 4: {
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
		 break;*/
	default:
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr_real_inv(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
	};
	return L;
}

void fft_radix_real_inv(int r, int N, int o) {

	print("{\n");
	print("// radix - %i\n", r);
	indent();
	for (int n = 0; n < N / r; n++) {
		if (n > N / r / 2) {
			continue;
		}
		print("{\n");
		indent();
		for (int i = 0; i < r; i++) {
			const int nn = n + i * N / r;
			print("const auto& zr%i = x[%i];\n", i, index_real_inv(o, nn, 0, N));
			if (nn != 0 && !(nn == N / 2 && N % 2 == 0)) {
				if (nn > N / 2) {
					print("const auto& zi%i = x[%i];\n", i, index_real_inv(o, nn, 1, N));
				} else {
					print("const auto zi%i = -x[%i];\n", i, index_real_inv(o, nn, 1, N));
				}
			} else {
				print("constexpr double zi%i = 0;\n", i);
			}
		}
		switch (r) {
		case 1:
			abort();
			return;
		case 2:
			print("const auto qr0 = zr0 + zr1;\n");
			print("const auto qi0 = zi0 + zi1;\n");
			print("const auto qr1 = zr0 - zr1;\n");
			print("const auto qi1 = zi0 - zi1;\n");
			break;
		case 3:
			print("const auto tr1 = zr1 + zr2;\n");
			print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
			print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
			print("const auto ti1 = zi1 + zi2;\n");
			print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
			print("const auto ti3 = (%24.17e) * (zi1 - zi2);\n", sqrt(3) * 0.5);
			print("const auto qr0 = zr0 + tr1;\n");
			print("const auto qi0 = zi0 + ti1;\n");
			print("const auto qr1 = tr2 + ti3;\n");
			print("const auto qi1 = ti2 - tr3;\n");
			print("const auto qr2 = tr2 - ti3;\n");
			print("const auto qi2 = ti2 + tr3;\n");
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
			print("const auto qr0 = tr1 + tr2;\n");
			print("const auto qi0 = ti1 + ti2;\n");
			print("const auto qr1 = tr3 + ti4;\n");
			print("const auto qi1 = ti3 - tr4;\n");
			print("const auto qr2 = tr1 - tr2;\n");
			print("const auto qi2 = ti1 - ti2;\n");
			print("const auto qr3 = tr3 - ti4;\n");
			print("const auto qi3 = ti3 + tr4;\n");
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
			print("const auto qr0 = zr0 + tr5;\n");
			print("const auto qi0 = zi0 + ti5;\n");
			print("const auto qr1 = tr8 + ti10;\n");
			print("const auto qi1 = ti8 - tr10;\n");
			print("const auto qr2 = tr9 + ti11;\n");
			print("const auto qi2 = ti9 - tr11;\n");
			print("const auto qr3 = tr9 - ti11;\n");
			print("const auto qi3 = ti9 + tr11;\n");
			print("const auto qr4 = tr8 - ti10;\n");
			print("const auto qi4 = ti8 + tr10;\n");
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
			print("const auto qr0 = tr7 + tr10;\n");
			print("const auto qr1 = tr8 - tr11;\n");
			print("const auto qr2 = tr9 + tr12;\n");
			print("const auto qr3 = tr7 - tr10;\n");
			print("const auto qr4 = tr8 + tr11;\n");
			print("const auto qr5 = tr9 - tr12;\n");
			print("const auto qi0 = ti7 + ti10;\n");
			print("const auto qi1 = ti8 - ti11;\n");
			print("const auto qi2 = ti9 + ti12;\n");
			print("const auto qi3 = ti7 - ti10;\n");
			print("const auto qi4 = ti8 + ti11;\n");
			print("const auto qi5 = ti9 - ti12;\n");
			break;

		default:
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
						print("double am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
						print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
					} else {
						print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
						print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
					}
				}
			}
			print("const auto qr0 = ", index_real_inv(o, n, 0, N / r));
			for (int i = 0; i < r; i++) {
				print_notab(" + zr%i", i);
			}
			print_notab(";\n");
			print("const auto qi0 = ", index_real_inv(o, n, 1, N / r));
			for (int i = 0; i < r; i++) {
				print_notab(" + zi%i", i);
			}
			print_notab(";\n");
			for (int i = 1; i <= (r - 1) / 2; i++) {
				print("const auto qr%i = ap%i + am%i;\n", i, i, i);
				print("const auto qi%i = bp%i - bm%i;\n", i, i, i);
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				print("const auto qr%i = ap%i - am%i;\n", r - i, i, i);
				print("const auto qi%i = bp%i + bm%i;\n", r - i, i, i);
			}
		}
		for (int q = 0; q < r; q++) {
			const auto W = twiddle(q * n, N);
			print("x[%i] = std::fma((%.17e), qr%i, (%.17e) * qi%i);\n", index_real_inv(o + q * N / r, n, 0, N / r), W.real(), q, -W.imag(), q);
			if (n != 0 && !(n == N / 2 / r && (N / r) % 2 == 0)) {
				print("x[%i] = -std::fma((%.17e), qr%i, (%.17e) * qi%i);\n", index_real_inv(o + q * N / r, n, 1, N / r), W.imag(), q, W.real(), q);
			}
		}

		deindent();
		print("}\n");
	}
	if (N > r) {
		for (int n = 0; n < r; n++) {
			fft_real_inv(N / r, o + n * N / r);
		}
	}
	deindent();
	print("}\n");
}

